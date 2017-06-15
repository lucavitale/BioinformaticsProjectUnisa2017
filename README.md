# BioinformaticsProjectUnisa2017
A project developed for the bioinformatics course at the University of Salerno 2016/2017. 

The goal of the project was to develop a classifier, based on pathways, to identify subclass of patients affected by tumors.

The proposed methodology is divided into four steps: 
  (i) Dimensionality reduction: since the gene expression data is high dimensional the DFP algorithm was used to identify the discriminant genes
  (ii) Pathways identification: the enrichment analysis is used to identify those biological pathways (KEGG or Reactome) that are statistically represented by the genes identified in the previous step
  (iii) Classification of the patients based on SVM: for each pathway a linear SVM is trained in cross validation. The most relevant genes are also identified. 
  (iv) Graphical representation: a graph is construct to represent the pathways relationships
  
The code was developed in R and Python.
