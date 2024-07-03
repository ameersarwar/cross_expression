# Cross-expression
## Introduction
Spatial transcriptomic technologies measure gene expression in individual cells while recording their physical locations in tissue. This allows us to ask how cells influence one another within the tissue. To address this question, we introduce cross-expression, a novel conceptual and statistical framework to understand coordinated gene expression between cells. Whereas co-expression measures the degree to which the expression of two genes is coordinated within the same cells, cross-expression quantifies how their expression is associated across neighboring cells. Here, we provide an efficient R package to perform cross-expression analysis, which makes pairwise comparisons between all genes and outputs a gene-by-gene p-value matrix indicating which pairs are significantly cross-expressed. Rather than targetting previously known genes, the cross-expression framework is unbiased, allowing us to discover spatial gene expression programs, thus leveraging the high-throughput spatial transcriptomics data to the fullest extent.
