# Cross-expression
## Introduction
Spatial transcriptomic technologies measure gene expression in individual cells while recording their physical locations, allowing us to ask how cells influence one another within the tissue. To address this question, we introduce cross-expression, a novel conceptual and statistical framework to understand coordinated gene expression between cells. Whereas co-expression measures the degree to which the expression of two genes is coordinated within the same cells, cross-expression quantifies how their expression is associated across neighboring cells. Since two genes can trivially cross-express if they are co-expressed in neighboring cells, we define cross-expression mutually exclusively, where the target cell expresses gene A but not gene B and its neighbor expresses gene B but not gene A, thus revealing genuine coordination as opposed to simple co-localization. Here, we provide an efficient R package to perform cross-expression analysis, which makes pairwise comparisons between all genes and outputs a gene-by-gene p-value matrix indicating which pairs are significantly cross-expressed across the tissue. Rather than targetting previously known genes, our framework therefore uses the high-throughput of spatial transcriptomics data to discover spatial gene expression programs.
![Main Figures](https://github.com/ameersarwar/cross_expression/assets/174621170/39552a8a-d29f-4a14-8949-05a6c3d0f01e)
## Cross-expression analysis
We will conduct the cross-expression analysis on a data collected using BARseq. Here, we provide one complete slice sectioned sagittally from the left hemisphere of an adult mouse brain. The data is in the `example_data` directory above and contains the gene expression matrix `expression.csv` and cell coordinates matrix `metadata.csv`.
### Download package and example data
To work through the demo, download the git repository by running the following commands in your terminal/shell `sh`:
```{sh}
git clone --depth 1 --filter=blob:none --sparse https://github.com/ameersarwar/cross_expression
cd cross_expression
git sparse-checkout set example_data
```
### Initialize R and load package plus data
In your terminal, run the folling command:
```{r}
touch demo_cross_expression.R
```
This will create an empty `R` file in which you can paste the upcoming commands and change the parameters to test various functions.

From now on, we will work exclusively in the `demo_cross_expression.R` file. Open it and set the working directory to `cross_expression`
```{r}
setwd("/path/to/root/directory/cross_expression/")
```
Run the package file `CrossExpression.R` from within `demo_cross_expression.R`
```{r}
source("CrossExpression.R")
```
This will load all the relevant functions into your global environment. (We will eventually make `CrossExpression.R` into a package that can be loaded simply as a library.)
Load the two `example_data` by running:
```{r}
# gene expression matrix
data = read.csv(file = "example_data/expression.csv")
rownames(data) = data$X; data = data[,2:ncol(data)]

# cell coordinates matrix
locations = read.csv(file = "example_data/metadata.csv")
rownames(locations) = locations$X; locations = locations[,2:ncol(locations)]
```
