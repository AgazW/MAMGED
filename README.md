# MAMGED (Meta-analyis for Microarray Gene Expression Data)
## Introduction

A tremendous amount of transcriptomic data has been generated, mostly via microarray experiments, for various species, tissues, cell-types, and conditions. Such data have accumulated in public repositories and have been often used by other researchers with various biological research goals such as identification of potential biomarkers, understanding the molecular mechanisms of diseases, the discovery of disease subtypes and predicting patient survival rate. Meta-analysis approach can aid in making a better use of existing gene expression data, particularly by increasing the overall sample size and countering the variations in results across studies. A promising meta-analysis method was reported a few years earlier to quantify the consistency of absolute expression calls (transcribed or not) made across experiments. However, biologists have not been able to make use of this method as there is no software/tool to apply the method to the data sets of interest. There was also a major limitation with that method in terms of deriving a consensus of differential levels of transcription associated with the conditions of interest. Hence, we developed MAMGED (Meta-Analysis of Microarray Gene Expression Data), a user-friendly tool to do a meta-analysis of microarray gene expression data across studies. This software, available off-line and online at [mamged.ibab.ac.in](http://mamged.ibab.ac.in/ "mamged online and offline version"), can be used to analyze the combined data from multiple studies, even across different types of microarray platforms, to derive a consensus expression status of genes in a specific tissue or cell-type and a condition of interest. In addition, the tool extends the application of the method to give a 'reliability score' for the consistency of differential levels of expression for each gene, and thus improves the original method. MAMGED is expected to aid biologists in making a better use of multiple datasets corresponding to a condition of interest by providing a hierarchical list of genes prepared based on the reproducibility of expression status and/or levels. The strength of association of genes with conditions indicated by the tool could be used to short-list potential biomarkers using existing datasets and open new ways for exploring molecular mechanisms of specific diseases.

[MAMGED](http://mamged.ibab.ac.in/ "mamged online and offline version") Is a tool written in [R](https://cran.r-project.org/) using [Rstudio](https://www.rstudio.com/products/RStudio/#Desktop) to perform meta-analysis of microarray  gene  expression  data   by
1. Quantifying  absolute  expression  calls.
2. Quantifying differential expression calls. 

The tool can be run online/offline. For online version, user need not to install any package ([R](https://cran.r-project.org/) and [Rstudio](https://www.rstudio.com/products/RStudio/#Desktop)) to use the tool. However, for using offline version, [R](https://cran.r-project.org/) and [Rstudio](https://www.rstudio.com/products/RStudio/#Desktop) need to be installed. More details can be found on [MAMGED Wiki](https://github.com/AgazW/MAMGED/wiki) and [MAMGED main page](http://mamged.ibab.ac.in/).

The current version of this tool is able to handle raw/processed data from three different platforms:
1. Affymetrix
* HC_G110 Affymetrix Human Cancer Array (GPL 74)
* Mu11KsubA Affymetrix Murine 11K SubA Array (GPL75)
* Hu6800 Affymetrix Human Full Length HuGeneFL Array (GPL80)
* HG_U95A Affymetrix Human Genome U95A Array (GPL 91)
* HG_U95B- Affymetrix Human Genome U95B Array (GPL92)
* HG_U95C- Affymetrix Human Genome U95C Array(GPL93)
* HG_U95D- Affymetrix Human Genome U95D Array (GPL94)
* HG_U95E- Affymetrix Human Genome U95E Array(GPL95)
* HG-U133B Affymetrix Human Genome U133B Array (GPL96)
* Hu35KsubA Affymetrix Human 35K SubA Array (GPL98)
* Hu35KsubB Affymetrix Human 35K SubB Array (GPL99)
* Hu35KsubC Affymetrix Human 35K SubC Array (GPL100)
* Hu35KsubD Affymetrix Human 35K SubD Array (GPL101)
* HG-Focus Affymetrix Human HG-Focus Target Array (GPL 201)
* HG-U133_Plus_2 (GPL 570)
* [HG-U133A_2] Affymetrix Human Genome U133A 2.0 Array (GPL571)
* [HT_HG-U133A] Affymetrix HT Human Genome U133A Array (GPL3921)
* HG_U95Av2 Affymetrix Human Genome U95 Version 2 Array (GPL 8300)
* [HG-U219] Affymetrix Human Genome U219 Array (GPL 13667)
* HG-U133_Plus_2 Affymetrix Gene Chip Human Genome HG-U133 Plus 2 Array Brainarray Version 12 (GPL 22321)
* HGU133Plus2_Hs_REFSEQ_12.1.0 (GPL10558)
* Affymetrix Human Genome U133 Plus 2.0 HGU133Plus2_Hs_REFSEQ_12.1.0] (GPL10881)
2. Codelink
* CodeLink Human Whole Genome Bioarray (GPL2895)
* CodeLink Human Whole Genome Bioarray (GPL6248)
* CodeLink Human Whole Genome Bioarray (GPL8547)
* CodeLink UniSet Human 20K I Bioarray (GPL2891)
* CodeLink UniSet Human 20K I Bioarray (GPL4044)
3. Illumina
* Illumina human-6 v2.0 expression beadchip (GPL6102)
* Illumina HumanHT-12 V3.0 expression beadchip (GPL694
* Illumina HumanHT-12 V4.0 expression beadchip (GPL10558)
* Illumina humanRef-8 v2.0 expression beadchip (GPL6104)
* Illumina HumanRef-8 v3.0 expression beadchip (GPL6883)
* Illumina HumanWG-6 v3.0 expression beadchip (GPL6884)


In addition to support for many more platforms from 1, 2 and 3, support for data from Agilent Platform will be available soon.


## PREREQUISITES

To use MAMGED offline version, the R 3.4 and Rstudio 1.1.383 or greater version has to be installed in the local machine. MAMGED offline requests the following R Packages that can be intalled using the `install.packages("packagename")` and then loaded by `library("packagename")`: 

NOTE: Replace `packagename` with the actual package name.

* shiny
* shinyjs
* simpleaffy
* affy
* oligo
* markdown
* plyr
* data.table
* dplyr
* codelink
* lumi
* h20kcod.db
* h10kcod.db
* hwgcod.db
* limma
* tibble
* dtplyr
* tools
* affyio
* gtools

After installing all the packages, shiny package must be loaded before running the MAMGED application. For  loading  shiny  R  package, enter `library(shiny)` in to the  Rstudio  consol.  All  other packages will be loaded automatically  once  the  application  is  started. 

## Running OF MAMGED

Once you are done with all the installation and loading process, put the tool folder in the current working directory. The working directory can also be set at the top of the `Rstudio session -> Set Working Directory -> Choose Directory`. After setting the path, check the current working  directory  by `getwd()` to  make  sure  that  correct  working  directory  is  set.  To  run  the program enter `runApp('applicationname)` i.e. the name of tool folder into R, open in browser and work. 

Note: Replace `applicationname` with the actual application name.
If  shiny  package  is  not  installed  and  an  attempt  is  made  to  run  the  application  by runApp('name  of  your  application’),  R  will  throw  an  error Error:  could  not  find  function "runApp". To  avoid such error,  make  sure  shiny  package  is  loaded  before  running  the application. 

## REPORTING ERRORS, FEATURE REQUEST, AND HELP
MAMGED is still a BETA version. If you are having any issue with the online/offline version of the MAMGED tool, please raise at https://github.com/AgazW/MAMGED/issues.

Feature requests and pull request are appreciated at https://github.com/AgazW/MAMGED/pulls and will be attended as soon as possible.

Send additional enquiries to mamged@ibab.ac.in
_____________________________________________________________________________________________________________________________

##### If you are having any issue with the online/offline version of the tool, please feel free to contact agazhussain@mangaloreuniversity.ac.in
