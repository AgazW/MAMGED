# MAMGED (Meta-analyis for Microarray Gene Expression Data)
## INTRODUCTION

A tremendous amount of transcriptomic data has been generated, mostly via microarray experiments, for various species, tissues, cell-types, and conditions. Such data have accumulated in public repositories and have been often used by other researchers with various biological research goals such as identification of potential biomarkers, understanding the molecular mechanisms of diseases, the discovery of disease subtypes and predicting patient survival rate. Meta-analysis approach can aid in making a better use of existing gene expression data, particularly by increasing the overall sample size and countering the variations in results across studies. A promising meta-analysis method was reported a few years earlier to quantify the consistency of absolute expression calls (transcribed or not) made across experiments. However, biologists have not been able to make use of this method as there is no software/tool to apply the method to the data sets of interest. There was also a major limitation with that method in terms of deriving a consensus of differential levels of transcription associated with the conditions of interest. Hence, we developed MAMGED (Meta-Analysis of Microarray Gene Expression Data), a user-friendly tool to do a meta-analysis of microarray gene expression data across studies. This software, available off-line and online at [mamged.ibab.ac.in](http://mamged.ibab.ac.in/ "mamged online and offline version"), can be used to analyze the combined data from multiple studies, even across different types of microarray platforms, to derive a consensus expression status of genes in a specific tissue or cell-type and a condition of interest. In addition, the tool extends the application of the method to give a 'reliability score' for the consistency of differential levels of expression for each gene, and thus improves the original method. MAMGED is expected to aid biologists in making a better use of multiple datasets corresponding to a condition of interest by providing a hierarchical list of genes prepared based on the reproducibility of expression status and/or levels. The strength of association of genes with conditions indicated by the tool could be used to short-list potential biomarkers using existing datasets and open new ways for exploring molecular mechanisms of specific diseases.

Three types of analysis can be performed for meta-analysis of absolute expression calls and differential expression calls:
1. Individual study analysis (Data samples belonging to one study).
2. Multi-study analysis (Data samples belonging to different studies by single platform).
3. Cross-platform study analysis (Data samples belonging to different studies and different platforms).

## Data that can be anayzed

The current version of the MAMGED is able to handle microarray gene expression raw/processed data from three different platforms:

1. Affymetrix
   1. Human
      - Affymetrix HT HG-U133+ PM Array Plate       
      - Affymetrix HT Human Genome U133 Array Plate Set(hthgu133a)     
      - Affymetrix HT Human Genome U133 Array Plate Set(hthgu133b)
      - Affymetrix hugene10(hugene10stprobeset) 
      - Affymetrix hugene10(hugene10sttranscriptcluster)       
      - Affymetrix hugene11(hugene11stprobeset)       
      - Affymetrix hugene11(hugene11sttranscriptcluster)       
      - Affymetrix hugene21(hugene21stprobeset)   
      - Affymetrix hugene21(hugene21sttranscriptcluster)    
      - Affymetrix HuGeneFL Genome Array(hu6800)       
      - Affymetrix Human Cancer G110 Array(hcg110)      
      - Affymetrix Human Genome 219 Plate(hgu219)     
      - Affymetrix Human Genome Focus Array(hgfocus)      
      - Affymetrix Human Genome HU35K Set(hu35ksuba)      
      - Affymetrix Human Genome HU35K Set(hu35ksubb)      
      - Affymetrix Human Genome HU35K Set(hu35ksubc)      
      - Affymetrix Human Genome HU35K Set(hu35ksubd)       
      - Affymetrix Human Genome U133 Plus 2.0 Array(hgu133plus2)     
      - Affymetrix Human Genome U133 Set(hgu133a)     
      - Affymetrix Human Genome U133 Set(hgu133b)      
      - Affymetrix Human Genome U133A 2.0 Array(hgu133a2)      
      - Affymetrix Human Genome U95 Set(hgu95a)      
      - Affymetrix Human Genome U95 Set(hgu95av2)      
      - Affymetrix Human Genome U95 Set(hgu95b)       
      - Affymetrix Human Genome U95 Set(hgu95c)      
      - Affymetrix Human Genome U95 Set(hgu95d)       
      - Affymetrix Human Genome U95 Set(hgu95e)     
      - [PrimeView] Affymetrix Human Gene Expression Array

   2. Mouse   
      - [HT_MG-430_PM] Affymetrix HT MG-430 PM Array Plate 
      - Affymetrix Mouse Expression Set 430(moe430a)                       
      - Affymetrix Mouse Expression Set 430(moe430b)
      - Affymetrix Mouse Genome 430A 2.0 Array(mouse430a2)
      - Affymetrix Mouse(mgug4121a)
      - Affymetrix Murine 11K Set(mu11ksuba)
      - Affymetrix Murine 11K Set(mu11ksubb)
      - Affymetrix Murine Genome U74v2(mgu74a)
      - Affymetrix Murine Genome U74v2(mgu74av2) 
      - Affymetrix Murine Genome U74v2(mgu74b)
      - Affymetrix Murine Genome U74v2(mgu74bv2)
      - Affymetrix Murine Genome U74v2(mgu74c)
      - Affymetrix Murine Genome U74v2(mgu74cv2)

      
   3. Rat   
      - Affymetrix Rat Expression Set 230(rae230a)
      - Affymetrix Rat Expression Set 230(rae230b)            
      - Affymetrix Rat Genome 230 2.0 Array(rat2302)
      - Affymetrix Rat Genome U34 Set(rgu34a)
      - Affymetrix Rat Genome U34 Set(rgu34b)
      - Affymetrix Rat Genome U34 Set(rgu34c)                      

   4. Arabidopsis
      - Affymetrix Arabidopsis ATH1 Genome Array(ath1121501) 
   5. Drosophila
      - Affymetrix Drosophila Genome 2.0 Array(drosophila2)
      - Affymetrix Drosophila Genome Array(drosgenome1)
   6. Yeast
      - Affymetrix Yeast Genome 2.0 Array(yeast2)
      - Affymetrix Yeast Genome S98 Array(ygs98)
   7. Zebrafish
      - Affymetrix zebrafish(zebrafish) 
   8. Celegans
      - Affymetrix celegans(celegans)    
      
2. Agilent
   1. Human
      - Agilent-012391 Whole Human Genome Oligo Microarray G4112A
      - Agilent-027478 Homo sapiens 420K HOTSPOT3.1_SG
      - Agilent-028004 SurePrint G3 Human GE 8x60K Microarray
      - Agilent-028005 SurePrint G3 Mouse GE 8x60K Microarray
      - Agilent-039494 SurePrint G3 Human GE v2 8x60K Microarray 039381
      - Agilent Human Genome,  Whole(hgug4112a)
      - Chips that use Agilent design number 026652(HsAgilentDesign026652)
   2. Mouse
      - Agilent Mouse Genome,  Whole(mgug4122a)
      - Chips that use Agilent design number 026655(MmAgilentDesign026655) 
   3. Rat
      - Chips that use Agilent design number 028282(RnAgilentDesign028282) 
   4. Celegans
      - Agilent-020186 C. elegans (V2) Gene Expression Microarray          
   5. Arabidopsis
      - Agilent-021169 Arabidopsis 4 Oligo Microarray (V4)
   6. E-coli
      - Agilent-029412 E. coli K-12 MG1655
      
3. Illumina
   1. Human
      - HumanHT12v3(illuminaHumanv3)                                       
      - HumanHT12v4(illuminaHumanv4)
      - HumanWGDASLv4(illuminaHumanWGDASLv4)
      - HumanHT12WGDASLv3
      - HumanWG6v2(illuminaHumanv2) 
   2. Mouse
      - MouseWG6v1p1(illuminaMousev1p1)
      - MouseWG6v2(illuminaMousev2)                                       


4. Codelink
   1. Human
      - UniSet Human 20k I Bioarray(h20kcod)
      - UniSet Human I Bioarray (~10 000 human genes)(h10kcod)             
      - Human Whole Genome Bioarray (~55 000 human genes)(hwgcod)          





Support for many more platforms from 1, 2 and 3, and Agilent Platform will be available soon.

## STRUCTURE OF MAMGED
The online version of the MAMGED consists of a web interface to a set of modules pipelined for analyses and a database holding GPL files. This allows the user to submit the data and set all the necessary parameters. Once the job is submitted, the user is provided with a job id for future reference. The submitted job is put in the queue and is processed on first come first serve basis. After the submitted job is analysed, the user can download the data by providing the job id (no online visualization of the loaded data and results is provided).  Due to large space requirements, the submitted job and the results are cleared after two weeks.  There is also a standalone version of the tool for remote work, which is downloadable from the [MAMGED](http://mamged.ibab.ac.in/v01/) under stand "Stand-Alone Version". In this version, the analysis process starts immediately after submitting the job (without putting the job in the queue) and the user is able to view the data and final results before downloading it for further use.

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
* affyio
* gtools


## SAMPLE DATA

The sample data can be downloaded directly for the [MAMGED](http://mamged.ibab.ac.in/v01/) under "Sample Data" section (which is recommended because the data contatins target files also, that are needed for differential expression analysis). The Affymetrix data sample [GSE45016](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE45016) was used to perform meta-analysis of absolute expression calls. A total of 78 samples from [GSE29721](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE29721), [GSE55092](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE55092), and [GSE6764](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE6764)  were used to perform meta-analysis of differential expression calls. 

Also prostate cancer gene expression samples from a study [GSE45016](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE45016), and annotation file [HG-U133_Plus_2] Affymetrix Human Genome U133 Plus 2.0 Array (GPL570) was used to validate individual study analysis). Data samples from three different studies of prostate cancer [GSE3325](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE3325), [GSE6369](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE6369), and [GSE9666](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE9666) were used in multi-study analysis; annotation was performed using GPL570. For cross-platform analysis, prostate cancer data samples of three different Affymetrix platforms [GSE7930](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE7930), [GSE9633](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE9633), and [GSE9666](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE9666),  using GPL96, GPL571 and GPL570 respectively were downloaded and used. 

## RUNNING OF MAMGED

[MAMGED](http://mamged.ibab.ac.in/ "mamged online and offline version") is written in [R](https://cran.r-project.org/) using [Rstudio](https://www.rstudio.com/products/RStudio/#Desktop) to perform meta-analysis of microarray  gene  expression  data   by
1. Quantifying  absolute  expression  calls.
2. Quantifying differential expression calls. 

The tool can be run online/offline. More details about usage and data can be found on [MAMGED](http://mamged.ibab.ac.in/).

#### Online Version
For online version, user need not to install any package. The user can directly access the [MAMGED](http://mamged.ibab.ac.in/) tool for meta-analysis. Once a job is submited, the user will be provided with a job id for future reference. The user will be intimated via emial about finishing of the job or any issue in the submited job. If the job is successfully processed, the results can be dowloaded by providing the job id.

#### Offline Version

To access MAMGED offline, the user needs to download the MAMGED package from the [MAMGED](http://mamged.ibab.ac.in/v01/) and install [R](https://cran.r-project.org/) and [Rstudio](https://www.rstudio.com/products/RStudio/#Desktop) on local computer. All other dependency R packages mentioned above need to be installed. Once all the packages are installed, the user needs to load the `shiny` R package by typing `library(shiny)` into the `Rstudio` console before running the tool. All the other packages will be automatically loaded once the application starts. As a final step to run the tool offline, tool folder needs to be put in the current working directory. The working directory can also be set at the top of the `Rstudio session -> Set Working Directory -> Choose Directory`. After setting the directory (path), check the current working  directory  by `getwd()` to  make  sure  that  correct  working  directory  is  set.  To  run  the program enter `runApp("applicationname")` i.e. the name of tool folder into R, open in browser and work.  In this version no job id is provided and user will be able to see and download the results once the analysis is finished. 

Note: Replace `applicationname` with the actual application name.
If  shiny  package  is  not  installed  and  an  attempt  is  made  to  run  the  application  by ``runApp("applicationname")`  ,  R  will  throw an error.

```diff
- Error:  could  not  find  function "runApp" 
``` 
To  avoid such error,  make  sure  shiny  package  is  loaded  before  running  the application. 

## REPORTING BUGS, REQUESTING FEATURES, AND HELP
MAMGED is still a BETA version. If you are having any issue with the online/offline version of the MAMGED tool, please raise at https://github.com/AgazW/MAMGED/issues.

Feature requests and pull request are appreciated at https://github.com/AgazW/MAMGED/pulls and will be attended as soon as possible.

For more help, reference manual is available under "Help" section of the [MAMGED](http://mamged.ibab.ac.in/v01/). Help about getting started with the stand alone version can be found under "Stand-Alone Version" section of the [MAMGED](http://mamged.ibab.ac.in/v01/).

Send additional enquiries to mamged@ibab.ac.in
_____________________________________________________________________________________________________________________________

