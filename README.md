##Transcriptomics data analysis

This repository contains tools to analysis the microarray data from the "System Biology in Practice" project. This course is part of the Bioninformatics and System Biology master program from the VU and UvA. 

Part of the code comes from the course material provided by Martijs Jonker.

#### Instructions:
To run these documents you need to place them in the same folder than "Expression.RMA.txt" and "design.txt". You also need to install the packages detailed later. When packages are installed, open the files in Rstudio and click in "Run document". 

Tools:

+ **pvalue-interactive.Rmd**: interactive document to peform t-test between genes of 2 groups of samples and visualize the results. Select one or more genes and get the expression levels in a barplot and t-test adjusted P-values.
+ **PCA.Rmd:** interactive document to perform PCA on all the samples or a subset. 


####Packages

You need some packages to run these tools.

+ limma from bioconductor
+ reshape2
+ dplyr
+ ggplot2
+ knitr
+ rmarkdown
+ shiny

*install packages using "install.packages("package")"*