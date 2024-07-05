# Cross-expression
Spatial transcriptomic technologies measure gene expression in individual cells and record their physical locations, allowing us to ask how cells influence one another within the tissue. To address this issue, we introduce cross-expression, a novel framework to understand coordinated gene expression between cells. Whereas co-expression measures the degree to which the expression of two genes is coordinated within the same cells, cross-expression quantifies how their expression is associated across neighboring cells. Since two genes can trivially cross-express if they are co-expressed in neighboring cells, we define cross-expression mutually exclusively, where the target cell expresses gene A but not gene B and its neighbor expresses gene B but not gene A, thus revealing genuine coordination as opposed to simple co-localization. Here, we provide an efficient R package to perform cross-expression analysis, which makes pairwise comparisons between all genes and outputs a gene-by-gene p-value matrix indicating which pairs are significantly cross-expressed across the tissue. Rather than targetting previously known genes, our framework therefore uses the high-throughput of spatial transcriptomics data to discover spatial gene expression programs in tissue.

![Main Figures](https://github.com/ameersarwar/cross_expression/assets/174621170/39552a8a-d29f-4a14-8949-05a6c3d0f01e)

## Part 1 - Demo setup
We will conduct cross-expression analysis on a dataset collected using BARseq. As an example, here we provide one complete slice sectioned sagittally from the left hemisphere of an adult mouse brain. This data is in the `example_data` directory above and contains the gene expression matrix `expression.csv` and cell coordinates matrix `metadata.csv`.

### Download package and example data
To work through the demo, the first step is to download the git repository by running the following commands in your terminal/shell `sh`:
```{sh}
git clone --depth 1 --filter=blob:none --sparse https://github.com/ameersarwar/cross_expression
cd cross_expression
git sparse-checkout set example_data
```
### Initialize R and load package plus data
In your terminal, run the following command:
```{r}
touch demo_cross_expression.R
```
This will create a new `R` file in which you can paste the upcoming commands and play around with various functions.

**!!Important!!**

From now on, we will work exclusively in the `demo_cross_expression.R` file. Open it and set the working directory to the `cross_expression` folder, making sure to change `/path/to/parent/directory/` to the parent directory of `cross_expression`.
```{r}
setwd("/path/to/parent/directory/cross_expression/")
```
Having set the working directory, run the following command from within `demo_cross_expression.R`
```{r}
source("CrossExpression.R")
```
This will load all the relevant functions into your global environment. (We will eventually make `CrossExpression.R` into a package that can be loaded simply as a library.)

Now, load the `example_data` matrices `expression.csv` and `metadata.csv` by running:
```{r}
# gene expression matrix
data = read.csv(file = "example_data/expression.csv")
rownames(data) = data$X; data = data[,2:ncol(data)]

# cell coordinates matrix
locations = read.csv(file = "example_data/metadata.csv")
rownames(locations) = locations$X; locations = locations[,2:ncol(locations)]
```
Check that the datasets loaded correctly by running:
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
cross = cross_expression(data = data, locations = locations)
head(cross)
```
The output is given as a dataframe:

<img width="676" alt="Screenshot 2024-07-05 at 12 44 51 AM" src="https://github.com/ameersarwar/cross_expression/assets/174621170/90b23fb7-987f-4ec2-8a6b-5240dddbc95f">

The function compares each gene pair and reports the p-values of cross-expression before (`cross_pvalue`) and after (`cross_padj`) Benjamini-Hochberg false discovery rate (FDR) multiple test correction. It also reports whether the p-values are significant (`cross_sig`) at `alpha ≤ 0.05` (after FDR) as well as the p-values of these genes' co-expression (`co_pvalue`) after FDR correction. You can play with the other parameters to test how the output changes.

The most important feature of `cross` is `cross_sig`, so let us only keep the gene pairs with significant cross-expression.
```{r}
cross = cross[as.logical(cross$cross_sig),]
head(cross)
```
Inspect the output:

<img width="717" alt="Screenshot 2024-07-05 at 1 04 13 AM" src="https://github.com/ameersarwar/cross_expression/assets/174621170/eb205442-3840-47c7-9dbd-a7ae587ba9a6">

The gene pairs in `cross` now show statistically significant cross-expression across our tissue slice. In total, you should have `52` pairs (out of `8778` possible non-directional pairs that can be formed using `133` genes in our panel). Since gene expression is regional, cross-expression is also regional, so understanding the input gene expression and cell coordinates is essential when interpreting the results of the `cross_expression` algorithm.

### Cross-expressing cells on tissue
We now have statistical evidence that, for the genes listed in `cross`, the expression of one gene in a cell predicts the expression of another gene in the neighboring cell. But spatial transcriptomics allows to see see gene expression in space.

To this end, let us color the cells based on the expression of `Tafa1` and `Col19a1`, the first gene pair in `cross`. Run the `tissue_expression_plot` function to do so:
```{r}
tissue_expression_plot(data = data, locations = locations, gene1 = "Tafa1", gene2 = "Col19a1", cross_expression = FALSE)
```
This shows the following image, where the cells are colored by the expression, co-expression, or non-expression of these genes:

<img width="1038" alt="Screenshot 2024-07-05 at 1 25 33 AM" src="https://github.com/ameersarwar/cross_expression/assets/174621170/abee8094-b8f3-4dba-a9ba-91feb60bd768">

However, it is still difficult to distinguish the cross-expressing cell-neighbor pairs from individual cells expressing each gene. To only color cross-expressing cells, call the `tissue_expression_plot` function while setting `cross_expression = TRUE`:
```{r}
tissue_expression_plot(data = data, locations = locations, gene1 = "Tafa1", gene2 = "Col19a1", cross_expression = TRUE)
```
This produces the following image, which clearly shows that these two genes are preferentially expressed in neighboring cells:

<img width="1039" alt="Screenshot 2024-07-05 at 2 02 43 AM" src="https://github.com/ameersarwar/cross_expression/assets/174621170/eccf7b7d-84ac-4d9c-968e-e1f5e75531e8">

One can view gene pairs of interest or customize the plot using `R`'s `ggplot` `library`, etc.

### Spatial enrichment of cross-expression
A salient feature of the last two images is that the cross-expressing cells are located towards the top of the slice (cortical brain regions) even though the individual genes, especially `Tafa1`, are expressed fairly broadly. This raises the question, "are cross-expressing cells spatially enriched"?

We can test for the hypothesis that the average distance between cross-expressing cells is smaller than that between cross-expressing and randomly selected cells. In other words, cross-expressing cells are spatially enriched. We test thus by running `spatial_enrichment`:
```{r}
enrich = spatial_enrichment(data = data, locations = locations, gene1 = "Tafa1", gene2 = "Col19a1")
enrich$pvalue
```
The p-value `5.839144e-06` from `enrich$pvalue` is smaller than `alpha = 0.05`, suggesting that cross-expression patterns between `Tafa1` and `Col19a1` are spatially enriched, confirming our observation that most such cell-neighbor pairs are towards the top of the tissue.

We can view the distances using:
```{r}
enrich$plot
```

This shows the following image:

<img width="969" alt="Screenshot 2024-07-05 at 2 32 47 AM" src="https://github.com/ameersarwar/cross_expression/assets/174621170/f9d77d36-e7be-40ce-92bd-e58d680a7789">

The `spatial_enrichment` function also contains the two distance distributions, which can be obtained via `enrich$target` and `enrich$null` for further analysis.


cross_expression_correlation()



## Part 3 - Auxiliary functions



bullseye_scores()
bullseye_plot()
rotate_coordinates()
smooth_cells()
