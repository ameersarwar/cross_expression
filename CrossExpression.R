# documentation

# install required packages if needed
if (!"Rfast" %in% rownames(installed.packages())){install.packages("Rfast", dependencies = TRUE, quiet = TRUE)}
if (!"RANN" %in% rownames(installed.packages())){install.packages("RANN", dependencies = TRUE, quiet = TRUE)}
if (!"Matrix" %in% rownames(installed.packages())){install.packages("Matrix", dependencies = TRUE, quiet = TRUE)}
if (!"tidyverse" %in% rownames(installed.packages())){install.packages("tidyverse", dependencies = TRUE, quiet = TRUE)}

# load required packages
library(Rfast); library(RANN); library(Matrix); library(tidyverse)

#' Computes Pearson's correlation between pairs of columns.
#' If one matrix is provided, the output is the pairwise correlations between its columns.
#' If two matrices are provided, the output is the pairwise correlations between their columns.
#'
#' @param matrix1 The first input matrix.
#' @param matrix2 The second input matrix (optional).
#'
#' @return Returns a correlation matrix.
#'
correlation <- function(matrix1, matrix2 = NULL){
  
  # correlation with columns of one matrix
  if (is.null(matrix2)){
    output = (t(scale(matrix1)) %*% scale(matrix1)) / (nrow(matrix1)-1)
  } else { # correlation with columns of two matrices
    output = (t(scale(matrix1)) %*% scale(matrix2)) / (nrow(matrix1)-1)
  }
  
  # return output
  return(output)
}

#' Calculates the number of elements common between columns of two matrices.
#' This function performs a simple dot product when binarize = FALSE.
#'
#' @param matrix1 The first input matrix.
#' @param matrix2 The second input matrix.
#' @param sparse If TRUE, coerces matrices to sparse format for efficiency.
#' @param binarize If TRUE, converts matrices to binary format to count the number of common elements.
#'
#' @return Returns a column by column matrix, where entries represent the number of common elements.
#'
get_cooccurrence_stats <- function(matrix1, matrix2, sparse = TRUE, binarize = TRUE){
  
  # convert to sparse matrices
  if (sparse){
    matrix1 = Matrix(data = Matrix::as.matrix(matrix1), sparse = TRUE)
    matrix2 = Matrix(data = Matrix::as.matrix(matrix2), sparse = TRUE)
  }
  
  # binarize for count-based calculation
  if (binarize){
    matrix1[matrix1 > 0] = 1
    matrix2[matrix2 > 0] = 1
  }
  
  # number of shared elements
  cooccur = t(matrix1) %*% matrix2
  cells   = t(matrix1) %*% matrix1
  
  # convert to regular matrix
  cooccur = base::as.matrix(cooccur)
  cells   = base::as.matrix(cells)
  
  # return output
  output  = list(cooccur,cells); names(output) = c("cooccur","cells")
  return(output)
}

#' Scales a vector between 0 and 1, including for negative values.
#'
#' @param vector The input vector to be scaled.
#'
#' @return Returns a vector scaled between 0 and 1.
#' 
scale_01 <- function(vector){
  vec       = as.vector(vector)
  min_value = min(vec, na.rm = TRUE)
  output    = (vec - min_value) / (max(vec, na.rm = TRUE) - min_value)
  return(output)
}

#' This function takes x and y coordinates and rotates them counterclockwise
#' by the specified number of degrees, and mean centers the points.
#'
#' @param x The x coordinates.
#' @param y The y coordinates.
#' @param n_degrees The degree of rotation in counterclockwise direction.
#' @param center If TRUE, mean centers the points.
#' @param scale If TRUE, standardizes the points.
#' @param flip_x Reflects the points across y = 0 line.
#' @param flip_y Reflects the points across x = 0 line.
#'
#' @return Returns a data frame containing the new x and y coordinates.
#'
rotate_coordinates <- function(x, y, n_degrees = 0, center = FALSE, scale = FALSE, flip_x = FALSE, flip_y = FALSE) {
  
  # convert degrees to radians
  theta <- n_degrees * (pi / 180)
  
  # apply rotation matrix to x and y coordinates
  x_rotated <- (cos(theta) * x) - (sin(theta) * y)
  y_rotated <- (sin(theta) * x) + (cos(theta) * y)
  
  # center and/or scales the coordinates
  x_rotated <- scale(x = x_rotated, center = center, scale = scale)
  y_rotated <- scale(x = y_rotated, center = center, scale = scale)

  # flip x or y coordinates
  if (flip_x){
    x_rotated <- -1 * x_rotated
  }
  if (flip_y){
    y_rotated <- -1 * y_rotated
  }
  
  # return output
  output <- data.frame(pos_x = x_rotated, pos_y = y_rotated)
  return(output)
}

#' Computes cross-expression and co-expression p-values between all gene pairs.
#'
#' @param data A cells by genes expression matrix.
#' @param locations A cells by coordinates (x-y or higher dimensions) matrix.
#' @param neighbor The n-th nearest neighbor for calculating cross-expression.
#' @param alpha_cross The significance threshold for cross-expression.
#' @param alpha_co The significance threshold for co-expression.
#' @param output_matrix If TRUE, outputs the cross-expression p-value matrix.
#'
#' @return Returns a list containing gene pairs with co-expression and cross-expression p-values
#'         before and after false discovery rate (FDR) correction.
#'

cross_expression <- function(data, locations, neighbor = 1, alpha_cross = 0.05, alpha_co = 0, output_matrix = FALSE){
  
  # convert data to sparse matrix and binarize
  data = Matrix(data = Matrix::as.matrix(data), sparse = TRUE)
  data[data > 0] = 1
  
  # find nth nearest neighbors and create neighbors x genes matrix
  distances = RANN::nn2(locations, locations, k = neighbor + 1, searchtype = "standard", treetype = "kd")
  distances = distances$nn.idx[, neighbor + 1]
  data_temp = data[distances,]
  
  # number of shared genes between cell-neighbor pairs (cross-expression) and within the same cells (co-expression)
  hyper_matrices = get_cooccurrence_stats(data, data_temp, sparse = FALSE, binarize = FALSE)
  cell_total     = nrow(data)
  
  # re-order matrices in ascending order by gene name
  names_seq = colnames(data)
  
  data      = data[,order(colnames(data))]
  data_temp = data_temp[,order(colnames(data_temp))]
  
  # co-expression p-values
  hyper_matrices$cells = hyper_matrices$cells[order(rownames(hyper_matrices$cells)), order(colnames(hyper_matrices$cells))]
  
  m = diag(hyper_matrices$cells)
  k = rep(m, each  = nrow(hyper_matrices$cells))
  m = rep(m, times = nrow(hyper_matrices$cells))
  n = cell_total - m
  
  q = as.vector(hyper_matrices$cells) - 1
  names(q) = colnames(hyper_matrices$cells)
  
  pval = phyper(q = q, m = m, n = n, k = k, lower.tail = FALSE)
  names(pval) = NULL
  pval = p.adjust(pval, method = "BH")
  
  co_expression = matrix(pval, nrow = nrow(hyper_matrices$cells), byrow = TRUE, dimnames = dimnames(hyper_matrices$cells))
  name_order    = match(names_seq, colnames(co_expression))
  co_expression = co_expression[name_order, name_order]
  
  # cross-expression p-values
  hyper_matrices$cooccur = hyper_matrices$cooccur[order(rownames(hyper_matrices$cooccur)), order(colnames(hyper_matrices$cooccur))]
  
  q = t(data * (1 - data_temp)) %*% ((1 - data) * data_temp)
  
  coexp = as.vector(hyper_matrices$cells)
  k = diag(hyper_matrices$cells)
  k = rep(k, each = nrow(hyper_matrices$cells))
  k = k - coexp
  
  neig = t(data_temp) %*% data_temp
  m = diag(neig)
  m = rep(m, times = nrow(hyper_matrices$cells))
  m = m - as.vector(neig)
  n = cell_total - m
  
  q = t(q)
  q = as.vector(q) - 1
  names(q) = rep(rownames(hyper_matrices$cells), times = nrow(hyper_matrices$cells))
  pval = phyper(q = q, m = m, n = n, k = k, lower.tail = FALSE)
  pvalx= pval # w/o FDR correction
  pval = p.adjust(pval, method = "BH")
  
  cross_expression = matrix(pval, nrow = nrow(hyper_matrices$cells), byrow = TRUE, dimnames = dimnames(hyper_matrices$cells))
  cross_expression = cross_expression[name_order, name_order]
  
  cross_x = matrix(pvalx, nrow = nrow(hyper_matrices$cells), byrow = TRUE, dimnames = dimnames(hyper_matrices$cells))
  cross_x = cross_x[name_order, name_order]
  
  # cross-expression is significant if A-B or B-A direction is significant
  cross_expression = Pmin(cross_expression, t(cross_expression))
  dimnames(cross_expression) = dimnames(co_expression)
  
  cross_x = Pmin(cross_x, t(cross_x))
  dimnames(cross_x) = dimnames(co_expression)
  
  # vectorized versions of co-expression and cross-expression p-values
  co_expr    = upper_tri(co_expression)
  cross_expr = upper_tri(cross_expression)
  cross_x    = upper_tri(cross_x)
  
  # co-expressing gene pairs are considered to not cross-express
  cross_expr[co_expr <= alpha_co] = 1
  cross_x[co_expr    <= alpha_co] = 1
  
  # assemble output with gene pairs and corresponding p-values for co-expression and cross-expression
  output = data.frame(which(upper.tri(cross_expression), arr.ind = TRUE))
  output = data.frame(output,
                      gene1 = colnames(cross_expression)[output$row],
                      gene2 = colnames(cross_expression)[output$col],
                      co_pvalue    = co_expr,
                      cross_pvalue = cross_x,
                      cross_padj   = cross_expr,
                      cross_sig    = as.integer(cross_expr <= alpha_cross))
  
  # assemble output as a gene-gene matrix of cross-expression p-values
  if (output_matrix){
    
    # output matrices with and w/o FDR correction
    mat_out = vector(mode = "list", length = 3)
    names(mat_out) = c("Cross_without_FDR","Cross_with_FDR","Co_with_FDR")
    
    # cross without FDR
    output_mat = matrix(0, nrow = nrow(co_expression), ncol = ncol(co_expression))
    dimnames(output_mat) = dimnames(co_expression)
    output_mat[upper.tri(output_mat)] = output$cross_pvalue
    output_mat = output_mat + t(output_mat)
    diag(output_mat) = 1
    
    mat_out[[1]] = output_mat
    
    # cross with FDR
    output_mat = matrix(0, nrow = nrow(co_expression), ncol = ncol(co_expression))
    dimnames(output_mat) = dimnames(co_expression)
    output_mat[upper.tri(output_mat)] = output$cross_padj
    output_mat = output_mat + t(output_mat)
    diag(output_mat) = 1
    
    mat_out[[2]] = output_mat
    
    # cross with FDR
    output_mat = matrix(0, nrow = nrow(co_expression), ncol = ncol(co_expression))
    dimnames(output_mat) = dimnames(co_expression)
    output_mat[upper.tri(output_mat)] = output$co_pvalue
    output_mat = output_mat + t(output_mat)
    diag(output_mat) = 1
    
    mat_out[[3]] = output_mat
    
    # return output
    return(mat_out)
  }
  
  # return output
  return(output)
}

#' Computes gene-gene correlations between cross-expressing cell-neighbor pairs.
#' Cell and neighbor masks are used to consider mutually exclusive expression per gene pair.
#'
#' @param data A cells by genes expression matrix.
#' @param locations A cells by coordinates (x-y or higher dimensions) matrix.
#' @param neighbor The n-th nearest neighbor for computing correlations.
#' @param output_matrix If TRUE, outputs the cross-expression correlation matrix.
#' 
#' @return Returns a gene list with cross-expression correlation for each gene pair.
#'
cross_expression_correlation <- function(data, locations, neighbor = 1, output_matrix = FALSE){
  
  # convert data to matrix
  data = as.matrix(data)
  
  # find nth nearest neighbors and create neighbors x genes matrix
  neighbor  = neighbor + 1
  distances = RANN::nn2(locations, locations, k = neighbor, searchtype = "standard", treetype = "kd")
  distances = distances$nn.idx[,neighbor]
  data_temp = data[distances,]
  
  # create masks for mutually exclusive gene expression
  mask_data = data
  mask_data[mask_data > 0] = 1
  
  mask_data_temp = data_temp
  mask_data_temp[mask_data_temp > 0] = 1
  
  X = mask_data * (1 - mask_data_temp)
  Y = (1 - mask_data) * mask_data_temp
  
  # keep mutually exclusive gene pairs
  X = X * data
  Y = Y * data_temp
  
  # gene pair correlations between cells and their neighbors
  corr = correlation(X, Y)
  corr = (corr + t(corr)) / 2
  
  # return correlations if matrix form requested
  if (output_matrix){
    return(corr)  
  }
  
  # return correlations as an edge list
  ids  = which(upper.tri(corr), arr.ind = TRUE)
  corr = data.frame(gene1 = rownames(corr)[ids[,1]], gene2 = colnames(corr)[ids[,2]], correlation = upper_tri(corr))
  
  if (!output_matrix){
    return(corr)
  }
}

#' Smooths cells' gene expression by averaging its expression by the nearest neighbors.
#' Optionally computes genes by genes Pearson's correlation matrix between cells by genes
#' and neighbors by genes matrices.
#'
#' @param data A cells by genes expression matrix.
#' @param locations A cells by coordinates (x-y or higher dimensions) matrix.
#' @param neighbors_smoothed Numbers of neighbors used for smoothing (0 is the cell itself; 1 is the nearest neighbor).
#' @param corr If TRUE, computes genes by genes correlation matrix between regions.
#' 
#' @return Returns a smoothed gene expression matrix. If corr = TRUE, additionally returns a gene-gene correlation matrix.
#'
smooth_cells <- function(data, locations, neighbors_smoothed = 1, corr = TRUE){
  
  # convolution kernel (0 is cell itself; 1 is first neighbor)
  neighbors_smoothed = neighbors_smoothed + 2
  
  # include nearest n cells (incl. cell itself) and exclude others
  distances = RANN::nn2(locations, locations, k = neighbors_smoothed, searchtype = "standard", treetype = "kd")
  distances = distances$nn.idx
  neighbor_regions = distances[,ncol(distances)]
  distances = distances[,ncol(distances) - 1]
  
  cellbycell = Matrix(0, nrow = nrow(locations), ncol = nrow(locations), doDiag = FALSE)
  
  # insert 1's for nearest neighbors (incl. cell itself) but 0's otherwise
  ids = cbind(rep(1:nrow(locations), neighbors_smoothed),
              as.vector(distances))
  cellbycell[ids] = 1
  
  # smooth gene expression
  data = Matrix(data = Matrix::as.matrix(data), sparse = TRUE)
  data_smooth = cellbycell %*% data
  data_smooth = data_smooth / neighbors_smoothed
  
  # gene x gene correlation between smoothed cells x genes matrix and smoothed neighbors x genes matrix
  # 'neighbors' are regions at n + 1 cells
  if (corr){
    corr = data_smooth[neighbor_regions,]
    corr = correlation(matrix1 = data_smooth, matrix2 = corr)
    
    # output smoothed expression matrix and correlation matrix
    data_smooth = base::as.matrix(data_smooth)
    output = list(data_smooth, corr); names(output) = c("smooth_expression", "smooth_correlation")
    
    # return output
    return(output)
  }
  
  # return output
  data_smooth = base::as.matrix(data_smooth)
  return(data_smooth)
}

#' Calculates bullseye statistics for ALL gene pairs.
#' Counts the number of cells with gene B in those with gene A
#' And in their neighbors
#' Neighbor scores are cumulative and normalized by window size.
#'
#' @param data A cells by genes expression matrix.
#' @param locations A cells by coordinates (x-y or higher dimensions) matrix.
#' @param window_sizes Vector of window sizes.
#' @param ratio_to_co If TRUE, results are normalized by co-expressing cells.
#' @param log_2 If TRUE, results are transformed by log2.
#' @param output_matrix If TRUE, results are provided in matrix form.
#'
#' @return Returns a dataframe with each row representing a gene pair and each column representing
#'         the bullseye score for that pair in each window size.
#'         output_matrix = TRUE presents a matrix for each window size.
#'
bullseye_scores <- function(data, locations, window_sizes = 1:10, ratio_to_co = FALSE, log_2 = FALSE, output_matrix = FALSE){
  
  # convert to sparse matrix and binarize
  data = Matrix(data = Matrix::as.matrix(data), sparse = TRUE)
  data[data > 0] = 1
  
  # find neighbors of each cell per window
  distances = RANN::nn2(locations, locations, k = max(window_sizes) + 1, searchtype = "standard", treetype = "kd")
  distances = distances$nn.idx
  
  # storage across window sizes
  outcome = matrix(data = 0, ncol = length(c(0,window_sizes)), nrow = choose(ncol(data),2) * 2)
  colnames(outcome) = c("Cell", str_c("Neig = ", window_sizes))
  
  # compute scores
  for (i in 2:ncol(outcome)){
    hits = t(data) %*% data[distances[,i],]
    outcome[,i] = c(hits[upper.tri(hits)], t(hits)[upper.tri(hits)])
  }
  
  # cumulative sums and normalize by window size
  outcome = t(apply(outcome, 1, cumsum))
  outcome = sweep(outcome, 2, c(1,window_sizes), FUN = "/")
  
  # add co-expression
  cells   = t(data) %*% data[distances[,1],]
  outcome[,1] = c(cells[upper.tri(cells)], t(cells)[upper.tri(cells)])
  
  # if ratio_to_co = TRUE, transform results as ratio
  if (ratio_to_co){
    outcome = sweep(outcome + 1, 1, outcome[,"Cell"] + 1, FUN = "/")
  }
  
  # if log_2 = TRUE, transform results as log2
  if (log_2){
    if (ratio_to_co) {outcome = log2(outcome)}
    if (!ratio_to_co){outcome = log2(outcome + 1)}
  }
  
  # append gene names
  upper_indices   = which(upper.tri(cells), arr.ind = TRUE)
  row = rownames(cells)[upper_indices[, 1]]
  col = rownames(cells)[upper_indices[, 2]]
  
  row1 = c(row, col)
  col1 = c(col, row)
  
  outcome = data.frame(gene1 = row1, gene2 = col1, outcome)
  colnames(outcome)[4:ncol(outcome)] = str_c("Neig = ", window_sizes)
  
  # if output_matrix = TRUE, then return outcome as matrix
  if (output_matrix){
    
    # storage list
    outcome_matrix = vector(mode = "list", length = length(c(1,window_sizes)))
    names(outcome_matrix) = colnames(outcome)[3:ncol(outcome)]
    
    # indices
    reconstructed_data = matrix(data = 0, nrow = ncol(data), ncol = ncol(data))
    upper_triangle_indices = which(upper.tri(reconstructed_data), arr.ind = TRUE)
    lower_triangle_indices = which(upper.tri(t(reconstructed_data)), arr.ind = TRUE)
    
    # create matrices
    for (i in 3:ncol(outcome)){
      
      reconstructed_data[upper_triangle_indices] = outcome[,i][1:length(upper_triangle_indices[,1])]
      reconstructed_data = t(reconstructed_data)
      reconstructed_data[lower_triangle_indices] = outcome[,i][(length(upper_triangle_indices[,1]) + 1):nrow(outcome)]
      reconstructed_data = t(reconstructed_data)
      
      rownames(reconstructed_data) = colnames(data); colnames(reconstructed_data) = colnames(data)
      outcome_matrix[[i-2]] = reconstructed_data
    }
    
    # return outcome matrix
    return(outcome_matrix)
  }
  
  # return outcome
  return(outcome)
}

#' Outputs a circular bullseye plot for a gene pair.
#' The central circle is gene B in cells expressing gene A.
#' Rings indicate neighbors with gene B, where the first ring is the first neighbor.
#'
#' @param scores Bullseye score as a vector.
#'
#' @return Returns a circular bullseye plot.
#'
bullseye_plot <- function(scores){
  
  # initialization
  vec_values   = as.numeric(scores)
  window_sizes = 1:length(vec_values)
  radius_ratio = 1
  
  # initial radius, ratio for subsequent circles, and no. of circles
  initial_radius = 1/length(window_sizes)
  n_circles      = length(window_sizes)
  
  # calculate the radii
  radii    = vector(mode = "numeric", length = length(window_sizes))
  radii[1] = initial_radius
  for (i in 2:n_circles) {radii[i] = radii[i-1] + initial_radius * radius_ratio^(i-1)}
  
  # dataframe for coordinates and hues of each circle
  df_list = list()
  for (i in 1:length(radii)) {
    theta = seq(0, 2 * pi, length.out = 100)
    if (i > 1) {inner_radius = radii[i - 1]} else {inner_radius = 0}
    outer_radius = radii[i]
    
    df_circle = data.frame(x = c(cos(theta) * inner_radius, rev(cos(theta) * outer_radius)),
                           y = c(sin(theta) * inner_radius, rev(sin(theta) * outer_radius)),
                           id = rep(i, length(theta) * 2), vec_value = rep(vec_values[i], length(theta) * 2))
    df_list[[i]] = df_circle
  }
  
  df = bind_rows(df_list); colnames(df) = c("x","y","window","hue")
  
  # create bullseye plot
  p = ggplot(df, aes(x = x, y = y, fill = hue)) +
    geom_polygon(aes(group = window)) + coord_fixed(ratio = 1) +
    scale_fill_gradient(low = "lightblue", high = "#08306B") +
    labs(x = "", y = "", fill = "") +
    theme_classic() +
    theme(axis.line = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          legend.text = element_text(size = 10))
  
  # return plot
  return(p)
}

#' Determines whether the supplied genes show spatial enrichment in cross-expression.
#' Spatial enrichment can be interpreted as delineating anatomical boundaries.
#'
#' @param data A cells by genes expression matrix.
#' @param locations A cells by coordinates (x-y or higher dimensions) matrix.
#' @param gene1 Name of gene1.
#' @param gene2 Name of gene2.
#' @param neighbor The n-th nearest neighbor.
#' @param max_pairs Specify maximum number of cell pairs to consider. Lower number increases computational efficiency.
#' 
#' @return Returns a p-value and distance distributions between cross-expressing cells and cross-expressing and random cells.
#'
spatial_enrichment <- function(data, locations, gene1, gene2, neighbor = 1, max_pairs = 20000){
  
  # subset genes and convert data to matrix
  data = data[,c(gene1,gene2)]
  data = as.matrix(data)
  data[data > 0] = 1
  
  # find nth nearest neighbors and create neighbors x genes matrix
  neighbor  = neighbor + 1
  distances = RANN::nn2(locations, locations, k = neighbor, searchtype = "standard", treetype = "kd")
  distances = distances$nn.idx[,neighbor]
  data_temp = data[distances,]
  
  # create masks for mutually exclusive gene expression
  mask_data = data
  mask_data_temp = data_temp
  
  X = mask_data * (1 - mask_data_temp)
  Y = (1 - mask_data) * mask_data_temp
  
  # keep mutually exclusive gene pairs
  X = X * data
  Y = Y * data_temp
  
  # average using locations of neighboring cells and append cross-expressing pairs
  locations = (locations + locations[distances,]) / 2
  locations = data.frame(locations, cross = as.numeric((X[,1] > 0 & Y[,2] > 0) | (X[,2] > 0 & Y[,1] > 0)))
  
  # subset by cross-expressing and random cells
  locations1 = locations[locations$cross == 1,]; locations1 = locations1[,c(1,2)]
  locations1 = locations1[sample(1:nrow(locations1), size = min(c(max_pairs, nrow(locations1))), replace = FALSE),]
  
  locations0 = locations[locations$cross == 0,]; locations0 = locations0[,c(1,2)]
  locations0 = locations0[sample(1:nrow(locations0), size = nrow(locations1), replace = FALSE),]
  
  # compute distances
  dist1 = Dist(locations1)
  dist1 = upper_tri(dist1)
  
  dist0 = rbind(locations1, locations0)
  dist0 = Dist(dist0)
  dist0 = dist0[1:nrow(locations1),(nrow(locations0)+1):ncol(dist0)]
  dist0 = as.numeric(dist0)
  
  # compute p-value and return p-value plus two distance distributions
  pval = wilcox.test(x = scale(dist1), y = scale(dist0), alternative = "less")
  pval = pval$p.value
  
  # make plot comparing distances
  target = dist1
  null   = dist0
  
  pp = data.frame(vals = c(target, null), type = rep(c("Cross-expressing", "Random"), times = c(length(target), length(null))))
  pp$type = factor(pp$type, levels = c("Random","Cross-expressing"))
  pp = ggplot(pp) + aes(x = vals, fill = type, y = after_stat(scaled)) +
    geom_density(alpha = 0.8) +
    labs(x = "Distance to cells", y = "Density", fill = "") + theme_classic() +
    guides(fill = guide_legend(reverse = TRUE)) +
    scale_fill_manual(values = c("Cross-expressing" = "lightblue", "Random" = "gray88")) +
    theme(legend.position = "top",
          axis.text = element_text(size = 10),
          axis.title = element_text(size = 12),
          legend.text = element_text(size = 12))
    
  # return p-value, two distance distributions, and plot
  pval = list(pvalue = pval, target = dist1, null = dist0, plot = pp)
  return(pval)
}

#' Plots gene expression and cross-expression on tissue by coloring cells.
#'
#' @param data A cells by genes expression matrix.
#' @param locations A cells by coordinates (x-y or higher dimensions) matrix.
#' @param gene1 Name of gene 1.
#' @param gene2 Name of gene 2.
#' @param cross_expression If TRUE, only cross-expressing cell pairs are shown.
#' @param neighbor The nearest neighbor (for cross-expression).
#' @param point_size Point size on the scatter plot.
#' @param scale_bar Length of the scale bad in microns.
#'
#' @return Returns a plot with cells shown as points and color indicating the genes it expresses.
#'
tissue_expression_plot <- function(data, locations, gene1, gene2, cross_expression = TRUE, neighbor = 1, point_size = 0, scale_bar = 0){
  
  # convert to sparse matrix and binarize
  data = data[,c(gene1,gene2)]
  gene_names = colnames(data)
  colnames(data) = c("gene1","gene2")
  data = Matrix(data = Matrix::as.matrix(data), sparse = TRUE)
  data[data > 0] = 1
  
  colnames(locations) = c("x","y")
  
  # scale bar
  x_start = as.numeric(quantile(locations$x, probs = 0.95))
  x_end   = x_start + scale_bar
  
  y_start = as.numeric(quantile(locations$y, probs = 0))
  y_end   = y_start
  
  # plot without cross-expression
  if (!cross_expression){
    
    # type of cell by gene expression
    type = vector(mode = "character", length = nrow(data))
    type[data[,1] == 1] = "Gene1"
    type[data[,2] == 1] = "Gene2"
    type[(data[,1] == 1) & (data[,2] == 1)] = "Both"
    type[(data[,1] == 0) & (data[,2] == 0)] = "Neither"
    
    # plot cells and color by type
    locations = data.frame(locations, type)
    locations$type = factor(locations$type, levels = c("Neither", "Both", "Gene2", "Gene1"))
    locations = locations[order(locations$type), ]
    
    p = ggplot(locations) + aes(x = x, y = y, color = type) +
      geom_point(size = point_size) +
      scale_color_manual(values = c("Neither" = "gray88", "Both" = "chartreuse3", "Gene2" = "deepskyblue4", "Gene1" = "brown3"),
                         labels = c("Neither" = "Neither", "Both" = "Both", "Gene2" = gene_names[2], "Gene1" = gene_names[1])) +
      labs(x = "x coordinates", y = "y coordinates", color = "") +
      guides(color = guide_legend(override.aes = list(size = 2), reverse = TRUE)) +
      annotate("segment", x = x_start, xend = x_end, y = y_start, yend = y_end, linewidth = 1, color = "black") +
      theme_classic() +
      theme(axis.text = element_text(size = 10),
            axis.title = element_text(size = 12),
            legend.text = element_text(size = 10))
    
    # return plot
    return(p)
  }
  
  # genes' cell-neighbor pairs
  neighbor  = neighbor + 1
  distances = RANN::nn2(locations, locations, k = neighbor, searchtype = "standard", treetype = "kd")
  distances = distances$nn.idx[,neighbor]
  df_neig   = data[distances,]
  
  cells     = rownames(data)
  neigs     = rownames(data)[distances]
  
  # cell-neighbor pairs to light up
  pair1 = data[,1] > 0 & data[,2] == 0 & df_neig[,2] > 0 & df_neig[,1] == 0
  pair2 = data[,2] > 0 & data[,1] == 0 & df_neig[,1] > 0 & df_neig[,2] == 0
  
  # cells lit up with genes 1 or 2
  gene1 = unique(c(cells[pair1], neigs[pair2]))
  gene2 = unique(c(cells[pair2], neigs[pair1]))
  
  # vector for lighting up cells with genes 1 or 2
  light = vector(mode = "character", length = nrow(data))
  light[rownames(data) %in% gene1] = "Gene1"
  light[rownames(data) %in% gene2] = "Gene2"
  light[light == ""] = "Neither"
  
  # process for plotting
  meta = data.frame(locations, light)
  meta$light = factor(meta$light)
  
  meta_neither = meta[meta$light == "Neither", ]
  meta_others  = meta[meta$light != "Neither", ]
  meta_reordered = rbind(meta_neither, meta_others)
  
  # plot cross-expression
  p = ggplot(meta_reordered) + aes(x = x, y = y, color = light) +
    geom_point(size = point_size) +
    scale_color_manual(values = c("Neither" = "gray88", "Both" = "chartreuse3", "Gene2" = "deepskyblue4", "Gene1" = "brown3"),
                       labels = c("Neither" = "Neither", "Both" = "Both", "Gene2" = gene_names[2], "Gene1" = gene_names[1])) +
    labs(x = "x coordinates", y = "y coordinates", color = "") +
    guides(color = guide_legend(override.aes = list(size = 2))) +
    annotate("segment", x = x_start, xend = x_end, y = y_start, yend = y_end, linewidth = 1, color = "black") +
    theme_classic() +
    theme(axis.text = element_text(size = 10),
          axis.title = element_text(size = 12),
          legend.text = element_text(size = 10))
  
  # return plot
  return(p)
}
