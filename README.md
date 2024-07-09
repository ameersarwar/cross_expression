# Cross-expression
Spatial transcriptomic technologies measure gene expression in individual cells and record their physical locations, allowing us to ask how cells influence one another within the tissue. To address this question, we introduce cross-expression, a novel framework for understanding coordinated gene expression between spatially adjacent cells. Whereas co-expression measures the degree to which the expression of two genes is coordinated within the same cells, cross-expression quantifies how their expression is associated across neighboring cells. Since two genes can trivially cross-express if they are co-expressed in neighboring cells, we define cross-expression mutually exclusively, where the target cell expresses gene A but not gene B and its neighbor expresses gene B but not gene A, thus revealing genuine coordination as opposed to mere co-localization. Here, we provide an efficient R package to perform cross-expression analysis, which makes pairwise comparisons between all genes and indicates which pairs are significantly cross-expressed across the tissue. Rather than targetting previously known genes, our framework uses the high-throughput of spatial transcriptomics data to discover spatial gene expression programs in tissue.

![Main Figures](https://github.com/ameersarwar/cross_expression/assets/174621170/39552a8a-d29f-4a14-8949-05a6c3d0f01e)

## Part 1 - Setup
We will conduct cross-expression analysis on a dataset collected using BARseq (barcoded anatomy resolved by sequencing). As an example, here we provide one complete slice sectioned sagittally from the left hemisphere of an adult mouse brain. This data is in the `example_data` directory above and contains the gene expression matrix `expression.csv` and cell coordinates matrix `metadata.csv`.

### Download package and example data
To work through the demo, the first step is to download (in the working directory of your choice) the git repository by running the following commands in your terminal:
```{sh}
git clone --depth 1 --filter=blob:none --sparse https://github.com/ameersarwar/cross_expression
cd cross_expression
git sparse-checkout set example_data
```
### Initialize R
In your terminal, run the following command:
```{r}
touch demo_cross_expression.R
```
This will create a new `R` file in which you can paste the upcoming commands and play around with various functions.

**!! Important !!**

From now on, we will work exclusively in the `demo_cross_expression.R` file.
1. Open `demo_cross_expression.R`, e.g., in `RStudio`.
2. Use `setwd("/path/to/parent/directory/cross_expression/")` to set the working directory. Make sure to change `/path/to/parent/directory/` to the parent directory of `cross_expression`.

For example, in my case the full path is `setwd("/Users/AmeerSarwar/Downloads/cross_expression/")`.

Having set the working directory, run the following command from within `demo_cross_expression.R`
```{r}
source("CrossExpression.R")
```
This will load all the relevant functions into your global environment. (We will eventually make `CrossExpression.R` into a package that can be loaded simply as a library.)

### Load data
Load the matrices `expression.csv` and `metadata.csv` from the `example_data` folder by running:
```{r}
# gene expression matrix
data = read.csv(file = "example_data/expression.csv")
rownames(data) = data$X; data = data[,2:ncol(data)]

# cell coordinates matrix
locations = read.csv(file = "example_data/metadata.csv")
rownames(locations) = locations$X; locations = locations[,2:ncol(locations)]
```
Check that the datasets are loaded correctly by running:
```{r}
data[1:5,1:5]
head(locations)
```
You should see the following:

<img width="354" alt="Screenshot 2024-07-04 at 2 32 29 AM" src="https://github.com/ameersarwar/cross_expression/assets/174621170/e29da77f-e57e-4b34-9bf9-41e95b3e0adf">

The `data` is a cells by genes expression matrix with 94,100 cells assayed across 133 genes, and `locations` contains the x (`pos_x`) and y (`pos_y`) coordinates (centroids) for each of the 94,100 cells.

## Part 2 - Core functions
Before analyzing our data, let us briefly look at it by running:
```{r}
ggplot(locations) + aes(x = pos_x, y = pos_y) + geom_point(size = 0) + theme_classic()
```
This outputs the image shown below, where each dot is a cell plotted using its x and y coordinates.

<img width="1153" alt="Screenshot 2024-07-05 at 12 32 00 AM" src="https://github.com/ameersarwar/cross_expression/assets/174621170/7cd1017e-3ef4-47d2-990f-349da44fd1ab">

### Cross-expression across all gene pairs
We will now perform cross-expression analysis, which tells us whether a gene is preferentially expressed in cells whose neighbors express another gene. The two main inputs are the gene expression matrix `data` and cell coordinates matrix `locations`. Run the function `cross_expression` and view its (default) output:
```{r}
cross = cross_expression(data = data, locations = locations, neighbor = 1, alpha_cross = 0.05, alpha_co = 0, output_matrix = FALSE)
head(cross)
```

The output is given as a dataframe:

<img width="676" alt="Screenshot 2024-07-05 at 12 44 51 AM" src="https://github.com/ameersarwar/cross_expression/assets/174621170/90b23fb7-987f-4ec2-8a6b-5240dddbc95f">

The function compares each gene pair and reports the p-values of cross-expression before (`cross_pvalue`) and after (`cross_padj`) Benjamini-Hochberg false discovery rate (FDR) multiple test correction. It also reports whether the p-values (after FDR) are significant (`cross_sig`) at `alpha ≤ 0.05` as well as the p-values of these genes' co-expression (`co_pvalue`) after FDR correction. You can play with the other parameters to test how the output changes.

A critical feature of `cross` is `cross_sig`, so let us only keep the gene pairs with significant cross-expression.
```{r}
cross = cross[as.logical(cross$cross_sig),]
head(cross)
```
Inspect the output:

<img width="717" alt="Screenshot 2024-07-05 at 1 04 13 AM" src="https://github.com/ameersarwar/cross_expression/assets/174621170/eb205442-3840-47c7-9dbd-a7ae587ba9a6">

The gene pairs in `cross` now only include those showing statistically significant cross-expression across our tissue slice. In total, you should have `52` pairs (out of `8778` possible pairs in the `133` gene panel).

**!! Important !!**

Gene expression is spatially variable, e.g., a gene may be expressed in some slices but not in others. Since cross-expression captures coordinated gene expression, it may also be spatially variable, e.g., a gene pair may cross-express in some slices but not in others. Accordingly, understanding the input gene expression and cell coordinates is essential when interpreting the results of the `cross_expression` algorithm.

**!! Important !!**

Here, we have assumed that a gene is expressed within a cell if its count is non-zero. However, one may think that gene expression is stochastic and apply a threshold to distinguish background expression from genuine transcription. Consequently, the user may input the thresholded matrix, and the `cross_expression` function will count the non-zero entries as instances of gene expression.

### Cross-expressing cells on tissue
We now have statistical evidence that, for the genes listed in `cross`, the expression of one gene in a cell predicts the expression of another gene in the neighboring cell. But spatial transcriptomics allows us to see gene expression in space.

To this end, let us color the cells based on the expression of `Tafa1` and `Col19a1`, the first gene pair in `cross`.

Run the `tissue_expression_plot` function to do so:
```{r}
tissue_expression_plot(data = data, locations = locations, gene1 = "Tafa1", gene2 = "Col19a1", cross_expression = FALSE, neighbor = 1, point_size = 0, scale_bar = 0)
```
This shows the following image, where the cells are colored by the expression, co-expression, or non-expression of these genes:

<img width="1038" alt="Screenshot 2024-07-05 at 1 25 33 AM" src="https://github.com/ameersarwar/cross_expression/assets/174621170/abee8094-b8f3-4dba-a9ba-91feb60bd768">

However, it is still difficult to distinguish the cross-expressing cell-neighbor pairs from individual cells expressing each gene. To only color cross-expressing cells, call the `tissue_expression_plot` function while setting `cross_expression = TRUE`:
```{r}
tissue_expression_plot(data = data, locations = locations, gene1 = "Tafa1", gene2 = "Col19a1", cross_expression = TRUE, neighbor = 1, point_size = 0, scale_bar = 0)
```
This produces the following image, which clearly shows that these two genes are preferentially expressed in neighboring cells:

<img width="1039" alt="Screenshot 2024-07-05 at 2 02 43 AM" src="https://github.com/ameersarwar/cross_expression/assets/174621170/eccf7b7d-84ac-4d9c-968e-e1f5e75531e8">

This allows one to view any gene pair of interest or customize the plot using `R`'s `ggplot` `library`, etc.

### Spatial enrichment of cross-expression
A salient feature of the last two images is that the cross-expressing cells are located towards the top of the slice (cortical brain regions) even though the individual genes, especially `Tafa1`, are expressed fairly broadly across the tissue. This raises the question, "are cross-expressing cells spatially enriched"?

We can test for the hypothesis that the average distance between cross-expressing cells is smaller than that between cross-expressing and randomly selected cells.

We do the test by running `spatial_enrichment`:
```{r}
enrich = spatial_enrichment(data = data, locations = locations, gene1 = "Tafa1", gene2 = "Col19a1", neighbor = 1, max_pairs = 20000)
enrich$pvalue
```
The p-value from `enrich$pvalue` is smaller than `alpha = 0.05`, suggesting that cross-expression patterns between `Tafa1` and `Col19a1` are spatially enriched, confirming our observation that most such cell-neighbor pairs are towards the top of the tissue.

**!! Important !!**

The `spatial_enrichment` algorithm is stochastic in the sense that we do not use all of the non-cross-expressing cells when computing distances. Instead, we use `max_pairs = 20000` as default. Using a larger number makes the p-value more accurate but reduces computational efficiency.

We can view the distances using:
```{r}
enrich$plot
```

This shows the following image:

<img width="969" alt="Screenshot 2024-07-05 at 2 32 47 AM" src="https://github.com/ameersarwar/cross_expression/assets/174621170/f9d77d36-e7be-40ce-92bd-e58d680a7789">

The `spatial_enrichment` function also contains the two distance distributions, which can be obtained via `enrich$target` and `enrich$null` for further analysis.

### Cross-expression correlation
Cross-expression tells us whether the expression of one gene in a cell predicts the expression of another gene in its neighbor. The `cross_expression` algorithm quantifies this in terms of probabilities, where p-values lower than `alpha` indicate significant spatial coordination given the genes' counts in the population at large.

However, this formalism does not provide a continuous metric of the strength of the spatial relationship. Specifically, it does not tell us whether cells with high expression of a given gene are neighbors of cells with similarly high (or low) expression of another gene.

To this end, we compute Pearson's correlation between genes across cells and neighbors. Whereas in co-expression the correlation is computed between gene expression vectors obtained from the same cells, in cross-expression the correlation is computed between gene expression vectors obtained from the cells and their neighbors. Like before, the cell-neighbor pairs must show mutually exclusive expression.

Find the correlations between gene pairs using `cross_expression_correlation`:
```{r}
corr = cross_expression_correlation(data = data, locations = locations, neighbor = 1, output_matrix = FALSE)
head(corr)
```
The output is as follows:

<img width="286" alt="Screenshot 2024-07-08 at 10 23 18 PM" src="https://github.com/ameersarwar/cross_expression/assets/174621170/5f28d3bd-848a-4530-b636-b5e858163f22">

The `cross_expression_correlation` function can be used in conjunction with `cross_expression`, e.g., by considering genes with significant cross-expression, or in isolation, e.g., comparing cross-expression correlations between nearby tissue sections.

## Part 3 - Auxiliary functions
We offer a few functions that augment the core analysis provided above. While these could be used independently, they will often be used as steps within longer analysis pipelines.

### Bulleseye scores - relative enrichment
In addition to statistical significance (`cross_expression`) and the strength of association (`cross_expression_correlation`), an important idea is `effect size`, which compares the number of neighbors expressing a gene to the number of target cells co-expressing it alongside the cognate gene.

The `effect size` - also called `bullseye scores` (see below) - are computed using:
```{r}
bull = bullseye_scores(data = data, locations = locations, window_sizes = 1:5, ratio_to_co = FALSE, log_2 = FALSE, output_matrix = FALSE)
head(bull)
```
This generates the output:

<img width="636" alt="Screenshot 2024-07-08 at 11 28 49 PM" src="https://github.com/ameersarwar/cross_expression/assets/174621170/f9959852-b604-48c5-8c34-dbd0246199f5">

We can see the scores for the target cell (`Cell`) and the neighbors `1:5` as specified in `window_sizes = 1:5`. The `window_sizes` input can be non-continuous and in any order, e.g., `window_sizes = c(100, 2:5, 12)`

The bullseye scores can be used for subsequent analysis. For example, one can present the scores as ratio of the neighbor to target scores (`ratio_to_co = TRUE` in the `bullseye_scores` function) and compare these across different genes or tissues, etc.

**!! Important !!**

Cross-expression is conceptualized without considering direction. Specifically, the p-values in `cross_expression` or correlations in `cross_expression_correlation` consider the cases where gene A or gene B is in the target cell. The outputs are the FDR-corrected lower p-values and average correlation, respectively.

Like the two cross-expression functions, the `bullseye_scores` consider both cases above. Unlike these functions, however, the `bullseye_scores` present both directions as outputs, making the output size here twice as large as in the two cross-expression functions. Thus, the user can average the directional information or take the minimum, etc., if a unidirectional output is desired for `bullseye_scores`.

### Bullseye plot

An effective way to visualize the effect size is in the form of bullseye plots. These require a vector of input scores, where the first score corresponds to the target cell and subsequent scores represent `window_sizes`. Therefore, the most straightforward way of making bullseye plots is to extract the scores from `bull` (which is storing the output of the `bullseye_scores` function).

Here, we make a bullseye plot for an example gene pair `Galntl6` (central cell) and `Nrg1` (neighbors):
```{r}
# extract scores for relevant genes
score_vector = as.numeric(bull[bull$gene1 %in% "Galntl6" & bull$gene2 %in% "Nrg1", 3:ncol(bull)])

# plot bullseye using scores
bullseye_plot(scores = score_vector)
```
This creates the following plot:

<img width="653" alt="Screenshot 2024-07-08 at 11 16 40 PM" src="https://github.com/ameersarwar/cross_expression/assets/174621170/ac0c53cc-cbd7-4b0e-9725-e1ce228139af">

The plot - called the `bullseye plot` - shows few cells co-expressing `Nrg1` with `Galntl6` (center) and many nearest neighbors (first ring) of `Galntl6`-positive cells expressing `Nrg1`. Importantly, the number of distant neighbors (subsequent rings) expressing `Nrg1` is about the same as the number of cells co-expressing both genes.

### Smooth gene expression
Our analysis considers cross-expression as a relationship between individual cells. However, groups of cells may form a spatially contiguous niche and cross-expression may obtain between niches.

We facilitate this analysis by smoothing (convolving) genes' expression in cells with that in their neighbors (kernel), which is a user-defined parameter.

Gene expression can be smoothed using:
```{r}
smth = smooth_cells(data = data, locations = locations, neighbors_smoothed = 5, corr = TRUE)
```
The `neighbors_smoothed = 0` performs no smoothing and any number above specifies the number of neighbors to use as the smoothing kernel. For example, `neighbors_smoothed = 5` smooths the gene expression using `5` nearest neigbhors. The `corr = TRUE` specifies that correlations between the 5-neighbor niches should be computed.

### Rotate tissue
During spatial transcriptomic data collection, tissues are often unaligned with respect to the glass slide or to each other. While sophisticated methods exist to non-linearly align them in the same coordinate plane, one is often interested in simply rotating, translating, or reflecting (linear transformations) the tissue to aid biological interpretation. For example, the cortical brain regions are usually presented towards the top of the slide.

```{r}
# rotate coordinates
locations2 = rotate_coordinates(x = locations$pos_x, y = locations$pos_y, n_degrees = 20, center = TRUE, scale = TRUE, flip_x = TRUE, flip_y = TRUE)

# view result
ggplot(locations2) + aes(x = pos_x, y = pos_y) + geom_point(size = 0) + theme_classic()
```

<img width="1061" alt="Screenshot 2024-07-09 at 12 40 55 AM" src="https://github.com/ameersarwar/cross_expression/assets/174621170/ed25f762-c505-4865-b7b5-e79e1c79834f">

The output in `locations2` rotates the xy coordinates `20` degrees counter-clockwise, z-score normalize them (`center = TRUE` and `scale = TRUE`), and reflect in the x and y axes (`flip_x = TRUE` and `flip_y = TRUE`).

In practice, the original tissue may have coordinates in a position resembling those in `location2`, thereby needing rotation in light with biological knowledge about the tissue and common practices regarding its presentation.

## Citation
If you find the cross-expression framework or the accompanying functions useful, please cite the following manuscript:

```{sh}
Sarwar A, Rue M, French L, Cross H, Chen X, & Gillis J.
Cross-expression analysis reveals patterns of coordinated gene expression in spatial transcriptomics
Journal Name (20xx). DOI-info.
```
