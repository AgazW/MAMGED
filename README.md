# MAMGED (Meta-analyis for Microarray Gene Expression Data)
[MAMGED](http://mamged.ibab.ac.in/ "mamged online and offline version") Is a tool written in R to perform meta-analysis  of  microarray  gene  expression  data  to
1. Quantify  absolute  expression  calls.
2. Quantify differential expression calls. 

The current version of this tool is able to handle raw/processed data from three different platforms:
1. Affymetrix
2. Codelink
3. Illumina

The process is to quantify the consistency of absolute expression calls (transcribed or not) made across experiments. The tool can be used to analyze the data from multiple studies across different types of microarray platforms to derive a consensus expression status of genes in a specific tissue or cell-type and a condition of interest. In addition, the tool extends the application of the method to provide a score for the consistency of differential expression pattern for each gene. MAMGED supports single study analysis, multiple study analysis and cross platform study analysis. The tool can be run online/offline. For online version user need not to install any package ([R](https://cran.r-project.org/) and [Rstudio](https://www.rstudio.com/products/RStudio/#Desktop)) to use the tool. However, for using offline version, [R](https://cran.r-project.org/) and [Rstudio](https://www.rstudio.com/products/RStudio/#Desktop) need to be installed. More details can be found on [MAMGED Wiki](https://github.com/AgazW/MAMGED/wiki) and [MAMGED main page](http://mamged.ibab.ac.in/).

___________________________________________________________________________________________________________________________
___________________________________________________________________________________________________________________________

h6(If you are having any issue with the online/offline version of the tool, please feel free to contact agazhussain@mangaloreuniversity.ac.in)
