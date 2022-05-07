# R---Microarray-Study-
A project that analyzes the gene expression dataset.



There are three groups in the main dataset: A=Sepsis B=Control C=SIRS. 

The main dataset contains the difference in gene expression between each group was measured. Gene expressions show how an organism responds to different conditions.

A binary significance value was assigned (1 significant, 0 not significant). 

Here I looked for genes with significant difference in Sepsis vs SIRS that do not have a significant difference in Control vs SIRS.

To reveal the specific genes that are upregulated or downregulated due to infection, I selected from siglist a group of genes with the following conditions: a fold change in gene expression between Sepsis vs Control, Sepsis vs SIRS greater than 1 and a significant difference between Control vs SIRS. 

These genes are the ones that are downregulated or upregulated due to infection

A second set of genes differentially expressed in Control vs SIRS was selected from siglist. Genes with a fold change in gene expression in Control vs SIRS greater than 2 were selected


Actions:
* Read siglist
* Extract the genes that have sigAC1=1 & sigBC=0
* Construct Venn Diagrams
* Export Files
* Fold change AB/AC>2 extraction and HEATMAP
* Fold change AB/AC>1 extraction and HEATMAP
* Fold change BC>2 extraction and HEATMAP 


The dataset is not uploaded as I do not have the rights to publicly distribute the dataset.
