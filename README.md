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
### Cross-expression demo
Create an `R` file by running the command:
```{r}
touch demo_cross_expression.R
```
Now, run the `CrossExpression.R` script in terminal/shell:
