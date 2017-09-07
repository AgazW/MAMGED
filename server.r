 
# This application Is to perform Meta-Analysis of gene expression data.
# In the initial phase we deal with Affymetrix, Codelink and Illumina data to get 
# the genes strongly associated with the disease condition. The genes selected 
# here are going to analysed according to their functions and pathway analysis in order to 
# find the potential cancer biomarkers. 
# In addition, the tool also perform differential expression of microarray gene expression data.
 

cat("\014")  


# Packages requried
#library(simpleaffy)
library(affy)
library(shiny)
library(markdown)
library(plyr)
library(data.table)
library(dplyr)
library(codelink)
library(lumi)
library(h20kcod.db)
library(h10kcod.db)
library(hwgcod.db)
library(limma)
library(shinyjs)
library(tibble)
#library(dtplyr)
#library(DBI)
#library(RMySQL)

#library(gcrma)
#library(hgu133a.db)
#library(hgu133acdf)
#library(hgu133plus2.db)
#library(hgu133plus2cdf)
#library(hgu133a2frmavecs)
#library(hgu133b.db)
#library(hgu133bcdf)
#library(hgu219.db)
#library(hgu219cdf)
#library(hgu95a.db)
#library(hgu95acdf)
#library(org.Hs.eg.db)
#library(hgu133a2.db)
#library(hgu133a2cdf)
#library(hgu95av2.db)
#library(hgu95av2cdf)
#library(hgu133plus2cdf)
#library(affyPLM)
##library(makecdfenv)
#library(parallel)
#library(base)
#library(S4Vectors)
#library(IRanges)
#library(stats4)
#library(BiocInstaller)
#library(Biobase)
#library(BiocGenerics)
#library(BiocParallel)
#library(biomaRt)
#library(Biostrings)
#library(preprocessCore)
#library(affyio)
#library(zlibbioc)
#library(graphics)
#library(grDevices)
#library(methods)
#library(genefilter)
#library(stats)
#library(AnnotationDbi)
#library(utils)

 

## Defining the size of file to be accepted. Here it can accept any size.
options(shiny.maxRequestSize= -1) 


## Function for affymetric processed data
affydata <- list()
affymetrix <- function(df, files, df1, affy_range){
  if(length(df) < 2){
    stop("Less than two files provided, please select 2 or more files to perform
         meta-analysis")
  }
  if(!is.null(affy_range)){
  vec <- as.numeric(unlist(strsplit(as.character(affy_range), split = "," )))
  print(vec)
  }else{
    vec <- length(df)
  }
  ## Through an error if length of range is greater than the
  ## length of annotation files.
  if(length(vec) > length(df1)) 
    stop("Range argument not a multiple of length of annotation")
  
  if(length(vec) < length(df1))
    stop("length of range is not equal to length of annotation ")
 
  if(any(vec== 0))
    stop("Range element is zero")
  
  if(sum(vec) != length(df))
    stop("length of range is not equal to length of input files")
  
  cat("testing more than one gpls ....\n")
  
  ## Performing annotation of input files with input gpl file/files
  ## The behavior of annotation is different if the sum of range is < length of input files.
  ## Eg. if vec = c(3,4), and length of input files is 10. If there are two input annotation files then
  ## the last 4 files will be annotated by 2nd annotation file and the remaining files will be annotated 
  ## by the first annotation file. if vec= (11,4) i,e the first vec element is greater than
  ## the length of input files, than all the files will be annotated by first annotation file.
  ## If vec = c(7,5), i,e the sum of vec is greater than the length of input files,
  ## then the first 7 files will be annotated by first annotation file and the remaning files
  ## will be annotated by the second file. 
  affy_combined <- do.call(c, 
                           Map(function(x,y) lapply(x, function(dat) 
                             merge(y, dat, by = c("ID"), all.y = TRUE)), 
                             split(df, rep(seq_along(vec), vec)),
                             df1))
  print(head(affy_combined[[1]]))
  
  for(i in 1:length(affy_combined)){
    incProgress(0.1, detail = paste("File", i)) 
    cat("Annotation done...",i,"\n")
    dataAnot <- na.omit(affy_combined[[i]][, -c(1,2)])
    #dataAnot <- dataAnot[, -c(1,2)]
    pvalue <- grepl("PVALUE", colnames(dataAnot))
    if(any(pvalue)){
      dt <- data.table(dataAnot)
      setkey(dt, PVALUE) ; indx <- dt[,.I[1L], by = SYMBOL]$V1 
      dataWoDuplic <- dt[indx]
    }else{
      cat("Trying to remove the duplicates\n")
      dataWoDuplic <- setDT(dataAnot)[order(factor(CALL, levels=c('P', 'A', 'M'))),
                                      .SD[1L], by = SYMBOL]
      #print(head(dataWoDuplic, 20))
    }
    affydata[[length(affydata)+1]] <- dataWoDuplic
  }
  affyfinal <- lapply(seq_along(affydata), function(i) setnames(affydata[[i]],
                          2:ncol(affydata[[i]]), paste0(names(affydata[[i]])[-1],i)))
  affyfinal <- Reduce(function(x, y) merge(x, y, by = c("SYMBOL"), all=TRUE), affyfinal)
  return(affyfinal)
}

###############################  End  ############################################## 


## Function to process Affymetrix raw  data.
processed <- list()
processed_data <- function(genAffyData){
  if(length(genAffyData)<2){
    stop("Less than two files provided, please select 2 or more files to perform
         meta-analysis")
  }
  for(i in 1: length(genAffyData))
  {
    incProgress(0.1, detail = paste("File", i))
    print(genAffyData[, i])
    cat("Length:", length(genAffyData), "\n")
    eset.mas5 = mas5(genAffyData[, i])
    
    ## getting the expression matrix (probesets/genes in rows, chips in columns).
    exprSet.nologs = exprs(eset.mas5)
    
    
    ## At this time let's log-transform the expression values to get a more normal distribution. 
    ## We have to remember we've done this when we calculate ratios. Logarithms can use any
    ## base, but base 2 is easiest when transforming ratios, since transformed 2-fold
    ## ratios up or down will be +1 or -1. As a result, we'll do all logs with base
    ## 2 to keep thing simplest.
    
    
    ## While we're doing Affymetrix-specific preprocessing, let's calculate
    ##an Absent/Present call for each probeset.
    # Run the Affy A/P call algorithm on the CEL files we processed above
    
    data.mas5calls = mas5calls(genAffyData[, i])
    
    data.mas5calls.calls = exprs(data.mas5calls)
     
    pvalue <- assayData(data.mas5calls)[["se.exprs"]]
    
    data.full <- cbind(exprSet.nologs, data.mas5calls.calls, pvalue)
    
    data.processed <- data.frame(data.full)
    
    #data.processed <- add_rownames(data.processed, "VALUE")
    data.processed <- rownames_to_column(data.processed, var = "VALUE")
    
    setnames(data.processed,c("ID","INTENSITY","CALL","PVALUE")) 
    
    processed[[length(processed)+1]] <- data.processed
    
    #cat("Affymetrix processing done...\n")
  }
  return(processed)
}

####################################  End  ###############################################


## Function for raw codelink data
codelink_processed <- list()

codelink_data <- function(codelinkpath){
  
  if(length(codelinkpath) < 2){
    stop("Less than two files provided, please select 2 or more files to perform
         meta-analysis")
  }
  for(i in 1:length(codelinkpath))
  {
    incProgress(0.1, detail = paste("File", i))  
    
    codset = readCodelinkSet(filename = codelinkpath[[i]])
    
    #print(head(codset))
    
    features <- readCodelink(codelinkpath[[i]])  ## An addition to get the ids
    
    ids <- features$id
    print(head(ids))
    
    if(all(!is.na(ids))){
      
      codset = codCorrect(codset, method = "half", offset = 0)
      
      codset = codNormalize(codset, method = "loess", weights = getWeight(codset),
                            loess.method = "fast")
      
      exprs <- exprs(codset)
      print(head(exprs))
      
      snr  <- getSNR(codset)
      print(head(snr))
      
      data <- data.frame(cbind(ids, exprs, snr), stringsAsFactors = FALSE)
      
      setnames(data, c("ID", "SIGNALINTENSITY", "SNR"))
      
      codelink_processed[[length(codelink_processed)+1]] <- data
      
      cat("Codelink processing done...\n") 
    }
  }
  if(length(codelink_processed) < 2){
    stop("Either none or less than two studies has id column present.
          Please provide data with ids present and try again.")
  }
  return(codelink_processed)
}


## Function for codelink processed data.
codelinkdata <- list()

codelink <- function(df, files, df1, codelink_range){
  if(length(df)<2){
    stop("Less than two files provided, please select 2 or more files to perform
         meta-analysis")
  }
  
  if(!is.null(codelink_range)){
    vec <- as.numeric(unlist(strsplit(as.character(codelink_range), split = "," )))
    print(vec)
  }else{
    vec <- length(df)
  }
  
  if(length(vec) > length(df1)) 
    stop("Range argument not a multiple of length of annotation")
  
  if(length(vec) < length(df1))
    stop("length of range is not equal to length of annotation ")
  
  if(any(vec== 0))
    stop("Range element is zero")
  print(length(df))
  
  ## this option is disapled here because some files are not valid 
  ## the invalid files are skiped if any
  #if(sum(vec) != length(df)) 
    #stop("length of range is not equal to length of input files")
  
  codelink_combined <- do.call(c, 
                           Map(function(x,y) lapply(x, function(dat) 
                             merge(y, dat, by = c("ID"), all.y = TRUE)), 
                             split(df, rep(seq_along(vec), vec)),
                             df1))
  print(head(codelink_combined[[1]]))
  
  for(i in 1:length(codelink_combined)) {
    incProgress(0.1, detail = paste("File", i))
    
    dataAnot <- na.omit(codelink_combined[[i]][, -c(1, 2)])
    
    #dataAnot <- dataAnot[, -c(1, 2)]    
    #print(head(dataAnot),10)
    
    dataAnot[with(dataAnot, order(SNR, decreasing = TRUE)),]
    
    dt <- data.table(dataAnot)
    
    #print(head(dt),10)
    
    ## Removing duplicate genes, set SNR as key column to remove the duplicates
    ## here max need to be kept, so use order() to sort the data frame and then
    ##removing the duplicates.
     
    indx <- dt[,.I[1L], by = SYMBOL]$V1 
    
    dataWoDuplic <- dt[indx]    
    #print(head(dataWoDuplic))
    
    codelinkdata[[length(codelinkdata)+1]] <- dataWoDuplic
  }
  ##Naming the columns uniquely and merging the list of data frames
  cdlnkdata <- lapply(seq_along(codelinkdata), 
                      function(i) setnames(codelinkdata[[i]],2:ncol(
                        codelinkdata[[i]]), paste0(names(codelinkdata[[i]])[-1],i)))
  
  cdlnkdata <- Reduce(function(x, y) merge(x, y, by = c("SYMBOL"),
                                           all=TRUE), cdlnkdata)
  
  return(cdlnkdata)
}



######################################  End  ######################################

## Function for Differential Expression of affymetrix data (Fold Change)
foldchange <- list()
fold_change_affy <- function(fileList, newpath, df1, fold, pval, newnames, conditions,AdditionalInfo,
                        FDRMethod){
  if(length(fileList) < 1){
    stop(" Please select at least one study with two 
         or more experimental conditions to perform differential expression")
  }
  ext <- grepl("*.txt$", fileList)
  for(i in 1:length(fileList))
  { 
    if(ext[i])
    {
      incProgress(0.1, detail = paste("File", i))  
      cond <- conditions
      cat("................................\n")
      
      ## Reading the phenotype files to construct the model matrix 
      lines <- readLines(newnames[i])
      
      pdata <- read.table(text=sub('.*(#.*)', '\\1', lines),
                         check.names=FALSE, stringsAsFactors=FALSE, 
                         comment.char='#',header=TRUE, row.names = NULL)
      print(head(pdata))
      
      raw.data <-  read.affy(fileList[i], path = newpath)
      
      x.rma <- call.exprs(raw.data, "rma")
      cat("testing affy .......\n")
      print(head(x.rma))   
      unique <- unique(pdata[, 2])
      f <- factor(pdata[, 2], levels = unique) 
      
      combination <- grep('(.*)-\\1',unique(as.vector(outer(unique, 
                                                            unique, FUN = paste,
                                                            sep = '-'))), value = TRUE, invert = TRUE)
      
      cat("possible combination:\n")
      conditions <- gsub(" ","", conditions[[1]])
      if(!all(conditions[[1]]%in%combination)){
       stop("Choosen contrasts is not a possible combination, please choose valid contrasts and try again")
      }   
      design <- model.matrix(~0 + f)
      colnames(design) <- unique  
      fit <- lmFit(x.rma, design)
      contrast.matrix <- makeContrasts(contrasts = conditions[[1]], levels = design)
      fit2 <- contrasts.fit(fit, contrast.matrix)
      fit2 <- eBayes(fit2)
      for(i in 1:ncol(contrast.matrix))
      {
        top <- topTable(fit2, coef = i, n = Inf, p.value = pval,
                        lfc = fold, adjust.method = FDRMethod)
        if(nrow(top) > 0)
        {
          setnames(top, c("AlogFC", "AveExpr", "t-statistic", "PVALUE", "AdjustedP-val",
                          "B-statistic"))
          if(all(AdditionalInfo == FALSE)){
            foldPval <- top[, c(1, 2, 4)]
            foldPval <- add_rownames(foldPval)
            names(foldPval)[1] <- "ID"
            #setnames(foldPval,c("ID", "AlogFC", "AveExpr", "PVALUE"))
          }
          else{
            foldPval <- top[, c(1, 2, 4, (which(AdditionalInfo[1] == TRUE)+2),
                                which(AdditionalInfo[2] == TRUE)+4,
                                which(AdditionalInfo[3]==TRUE)+5)]
            cat("checking the extracted columns\n")
            foldPval <- add_rownames(foldPval)
            names(foldPval)[1] <- "ID"
          }
           
          #setnames(foldPval,c("ID", "AlogFC", "AveExpr", "PVALUE"))
          cat("testing differential expressin .........\n")
          print(head(foldPval))
          foldchange[[length(foldchange)+1]] <- foldPval
        }
      }
    }
  }
  if (length(foldchange) < 1){
    stop("No gene found significant at the threshold p-value= ", pval ,
         "and foldchange = ",fc, " please change your pvalue 
         and foldchange and try again")
  }
  return(foldchange)
}

## Function to annotate for fold change, linked to above code
foldchange1 <- list()
affy_combined_fc <- function(df, df1, affy_diff_range){
  if(!is.null(affy_diff_range)){
    vec <- as.numeric(unlist(strsplit(as.character(affy_diff_range), split = "," )))
    print(length(df))
    print(vec)
  }else{
    vec <- length(df)
  }
  ## Through an error if length of range is greater than the
  ## length of annotation files.
  if(length(vec) > length(df1)) 
    stop("Range argument not a multiple of length of annotation")
  
  if(length(vec) < length(df1))
    stop("length of range is not equal to length of annotation ")
  
  if(any(vec== 0))
    stop("Range element is zero")
  
  if(sum(vec) != length(df))
    stop("length of range is not equal to length of input files")
  
  cat("testing more than one gpls ....\n")
  
  affy_fc_annotated <- do.call(c, 
                           Map(function(x,y) lapply(x, function(dat) 
                             merge(y, dat, by = c("ID"), all.y = TRUE)), 
                             split(df, rep(seq_along(vec), vec)),
                             df1))
  
  for(i in 1:length(affy_fc_annotated))
  {
    incProgress(0.1, detail = paste("File", i))
    cat("Testing affymetrix diff expr ......\n")
    print(head(df[i]))
    print(length(affy_fc_annotated))
    print(head(affy_fc_annotated[1]))
    afterAnot <- na.omit(affy_fc_annotated[[i]][, -c(1,2)])
    dt <- data.table(afterAnot)
    setkey(dt, PVALUE) ; indx <- dt[,.I[1L], by = SYMBOL]$V1 
    dataWoDuplic <- dt[indx]
    foldchange1[[length(foldchange1)+1]] <- dataWoDuplic
  }
  foldchange <- lapply(seq_along(foldchange1), function(i) setnames(foldchange1[[i]],
                                                               2:ncol(foldchange1[[i]]), 
                                                               paste0(names(foldchange1[[i]])[-1],i)))
  foldchange <- Reduce(function(x, y) merge(x, y, by = c("SYMBOL"), all=TRUE), foldchange)
}
###########################  End ######################################



## Function for differential expression of codlink data 
coddiff <- list()
codelinkdiff <- function(fileList, newpath, cntrsts, pval, fc, annotation, newnames,
                         AdditionalInfo, FDRMethod){
  ext <- grepl("*.txt$", fileList)
  
  # Removing extra space to avoid error in contrasts.............. 
  cntrsts <- gsub(" ","", cntrsts[[1]])
  if(length(fileList) < 1){
    stop("Please select at least one study with two
         or more experimental conditions to perform differential expression")
    }
  if(!any(ext)){
    stop("No target files found, please check target files and try again")
    }
  if(all(ext)){
    stop("No data files found, please check input data files and try again")
    }
  for(i in 1:length(fileList))
    { 
    if(ext[i])
      {
      cat("Differential expression function called \n")
      incProgress(0.1, detail = paste("File", i))
      ## Reading phenotype data to construct model matrix
      lines <- readLines(newnames[i])
      pdata <- read.table(text=sub('.*(#.*)', '\\1', lines),
                          check.names=FALSE, stringsAsFactors = FALSE, 
                          comment.char='#', header=TRUE)
      #print(pdata)
      unique <- unique(pdata[, 2])
      f <- factor(pdata[, 2], levels = unique)
      cat("Treatment conditions:", unique, "\n")
      combination <- grep('(.*)-\\1',unique(as.vector(outer(unique, unique,
                                                            FUN = paste, sep = '-'
                                                            ))),value = TRUE, invert = TRUE)
      if(!all(cntrsts[[1]]%in%combination)){
        stop("Choosen contrast is not a possible combination,
             please choose the possible combination and try again")
        }
      pdata$FileName <- file.path(newpath, pdata[,1])  ## This part is trick(Files should follow the same pattern)
      cat("setting path done..\n")
      codset <- readCodelinkSet(filename = pdata$FileName) 
      codset = codCorrect(codset, method = "half", offset = 0)
      codset = codNormalize(codset, method = "loess",
                            weights = getWeight(codset), loess.method = "fast")
      design <- model.matrix(~0 + f)
      cat("Experimental design:\n",design)
      colnames(design) <- unique ## c("C","TT","D","DT")
      fit <- lmFit(codset, design,weights = getWeight(codset))
      contrast.matrix <- makeContrasts(contrasts = cntrsts[[1]], levels = design)
      fit2 <- contrasts.fit(fit, contrast.matrix)
      fit2 <- eBayes(fit2)
      for(i in 1:ncol(contrast.matrix))
        {
        top <- topTable(fit2, coef = i, n = Inf, p.value = pval,
                        lfc = fc, adjust.method = FDRMethod)
        cat("I am checking for top\n")
        print(top)
        if(nrow(top) > 0)
          {
          top <- top[,c("probeName","meanSNR", "logFC", "AveExpr", "t", "P.Value" ,"adj.P.Val", "B")]
          if(all(AdditionalInfo == FALSE)){
            foldPval <- top[, c("probeName", "logFC", "AveExpr", "P.Value")]
            setnames(foldPval, c("PROBEID", "ClogFC", "AveExpr", "PVALUE"))
            }
          else{
            top <- top[,c("probeName","meanSNR", "logFC", "AveExpr", "t", "P.Value" ,"adj.P.Val", "B")]
            setnames(top, c("PROBEID", "meanSNR", "ClogFC", "AveExpr", "t-statistic", "PVALUE",
                            "AdjustedP-val", "B-statistic"))
            foldPval <- top[, c(1, 3, 4, 6, (which(AdditionalInfo[1] == TRUE)+1),
                                (which(AdditionalInfo[2] == TRUE)+4),
                                (which(AdditionalInfo[3]==TRUE)+6),
                                (which(AdditionalInfo[4] == TRUE)+7))]
            }
          ## Extracting information from ids to have a match with accnum 
          ids <- foldPval$PROBEID
          test <- grepl("_", ids[1])
          newids <- gsub("\\..*|_PROBE.*", "", ids)
          test1 <- grepl("GE", newids[1])
          if(test)
            {
            foldPval <- cbind(newids, foldPval[, -1])
            colnames(foldPval)[1] <- "PROBEID"
            }
          if(annotation == "hwgcod")
            {
            if(test){
              keys <- AnnotationDbi::select(hwgcod.db, newids, c("SYMBOL"), "ACCNUM")
              }
            else if(test1){
              keys <- AnnotationDbi::select(hwgcod.db, newids, c("SYMBOL"), "PROBEID")
            }
            }
          else if(annotation == "h20kcod")
            {
            if(test){
              keys <- AnnotationDbi::select(h20kcod.db, newids, c("SYMBOL"), "ACCNUM")
              }
            else if(test1){
              keys <- AnnotationDbi::select(h20kcod.db, newids, c("SYMBOL"), "PROBEID")
            }
            }
          else if(annotation == "h10kcod")
            {
            if(test){
              keys <- AnnotationDbi::select(h10kcod.db, newids, c("SYMBOL"), "ACCNUM")
              }
            else if(test1){
              keys <- AnnotationDbi::select(h10kcod.db, newids, c("SYMBOL"), "PROBEID")
            }
            }
          setnames(keys, c("PROBEID", "SYMBOL"))
          cat("keys form database:")
          #print(head(keys))
          fold <- merge(keys, foldPval, by = c("PROBEID"), all.x = TRUE)
          fold1 <- fold[, -1]
          fold1 <- na.omit(fold1)
          cat("fold is:\n")
          #print(head(fold))
          #print(head(fold1))
          dt <- data.table(fold1)
          setkey(dt, PVALUE) ; indx <- dt[,.I[1L], by = SYMBOL]$V1
          dataWoDuplic <- dt[indx]
          coddiff[[length(coddiff)+1]] <- dataWoDuplic
        }
      }   
    }
    }
  if (length(coddiff) == 0){
    stop("No gene found significant at the threshold p-value= ", pval," 
         and foldchange = ", fc,  " please change 
         your pvalue and foldchange and try again")
    }
  if (length(coddiff) == 0){
    stop("Files selected belong to only one study,
         please provide atleast two studies to
         perform meta-analysis of 
         differential expression")
    }
  return(coddiff)
}

### Linked to codelink diff expression (Above code).
codlinkdiff1 <- function(df){
  codfinal <- lapply(seq_along(df), function(i) setnames(df[[i]],
                                                         2:ncol(df[[i]]),
                                                         paste0(names(df[[i]])[-1],i)))
  codfinal <- Reduce(function(x, y) merge(x, y, by =c("SYMBOL"), all = TRUE), codfinal)
  cat("testing the flow\n")
  print(head(codfinal))
  codfinal
}
############################## End ###################################


## Function for illumina processed data to perform meta-analysis
illuminadata <- list()
illumina <- function(df, files, df1, illumina_range){
  if(length(df)<2){
    stop("Less than two files provided,please select 2 or more 
         files to perform meta-analysis")
  }
  
  if(!is.null(illumina_range)){
    vec <- as.numeric(unlist(strsplit(as.character(illumina_range), split = "," )))
    print(vec)
  }else{
    vec <- length(df)
  }
  ## Through an error if length of range is greater than the
  ## length of annotation files.
  if(length(vec) > length(df1)) 
    stop("Range argument not a multiple of length of annotation")
  
  if(length(vec) < length(df1))
    stop("length of range is not equal to length of annotation ")
  
  if(any(vec== 0))
    stop("Range element is zero")
  
  if(sum(vec) != length(df))
    stop("length of range is not equal to length of input files")
  
  cat("testing more than one gpls ....\n")
  
  illumina_combined <- do.call(c, 
                           Map(function(x,y) lapply(x, function(dat) 
                             merge(y, dat, by = c("ID"), all.y = TRUE)), 
                             split(df, rep(seq_along(vec), vec)),
                             df1))
  
  
  for(i in 1:length(illumina_combined)) {
    cat("Illumina in annotation section:\n")
    incProgress(0.1, detail = paste("File", i))
    dataAnot <- na.omit(illumina_combined[[i]][, -c(1, 2)])
    dt <- data.table(dataAnot)
    setkey(dt, PVALUE) ; indx <- dt[,.I[1L], by = SYMBOL]$V1 
    dataWoDuplic <- dt[indx]
    illuminadata[[length(illuminadata) + 1]] <- dataWoDuplic
    }
  ilmnfinal <- lapply(seq_along(illuminadata), 
                      function(i) setnames(illuminadata[[i]],
                                           2:ncol(illuminadata[[i]]),
                                           paste0(names(illuminadata[[i]])[-1],i)))
  ilmnfinal <- Reduce(function(x, y) merge(x, y, by = c("SYMBOL"), all=TRUE), ilmnfinal)
  return(ilmnfinal)
}
########################### End ################################################


## Function to process Illumina raw data
ilm <- list()
illuminaRaw <- function(illumina){
  for(i in 1:length(illumina))
    {
    incProgress(0.1, detail = paste("File", i))
    x.lumi <- lumiR(illumina[i])
    lumiExpr <- lumiExpresso(x.lumi, bg.correct = TRUE, 
                             normalize = TRUE, verbose = TRUE)
    exprs <- exprs(lumiExpr)
    pvalue <- detection(lumiExpr)
    #print(head(exprs))
    #print(head(pvalue))
    if(length(exprs) < 3)
      {
      lumidata <- cbind(exprs, pvalue)
      data.processed <- data.frame(lumidata)
      data.processed <- add_rownames(data.processed, "VALUE")
      colnames(data.processed) <- c("ID", "INTENSITY", "PVALUE")
      ilm[[length(ilm) + 1]] <- data.processed
      }else{
        exprs <- data.frame(exprs)
        #print(head(exprs))
        pvalue <- data.frame(pvalue)
        exprs <- add_rownames(exprs, "ID")
        pvalue <- add_rownames(pvalue, "ID")
        names(exprs)[2:ncol(exprs)] <- paste("INTENSITY", 1:ncol(exprs), sep = "")
        names(pvalue)[2:ncol(pvalue)] <- paste("PVALUE", 1:ncol(pvalue), sep = "")
        length <- ncol(exprs)-1
        for(i in 1:length)
          {
          exprs1 <- exprs[, c(1, i + 1)]
          pvalue1 <- pvalue[, c(1, i + 1)]
          exp.pval <- merge(exprs1, pvalue1, by = "ID")
          setnames(exp.pval, c("ID", "INTENSITY", "PVALUE"))
          ilm[[length(ilm)+1]] <- exp.pval
        }
      }
    }
  cat("Illumina raw data processing done ... ")
  return(ilm)
  }
################################# End ######################################


## Function for differential expression of illumina data
ilmndiff <- list()
illuminaDe <- function(ilmNames, ilmTargets, x, fc, pval, newnames,
                       AdditionalInfo, FDRMethod){
  # Removing extra space to avoid error in contrasts.............. 
  x <- gsub(" ","", x[[1]])
  if(length(ilmNames) && length(ilmTargets) < 1){
    stop("Please select at least one study with two
         or more experimental conditions to 
         perform differential expression")
    }
  if(length(ilmTargets) == 0){
    stop("No target files found, please check target files and try again")
    }
  if(length(ilmNames) == 0){
    stop("No data files found, please check input data files and try again")
    }
  for(i in 1:length(ilmNames))
    {
    incProgress(0.1, detail = paste("File", i))  
    lines <- readLines(ilmTargets[[i]])
    pdata <- read.table(text = sub('.*(#.*)', '\\1', lines), check.names = FALSE,
                        stringsAsFactors = FALSE, comment.char = '#', header = TRUE)
    unique <- unique(pdata[, 2])
    f <- factor(pdata[, 2], levels = unique)
    target <- readTargets(ilmTargets[[i]])
    #cat("Targets:\n")
    ##print(target)
    combination <- grep('(.*)-\\1', unique(as.vector(outer(unique, unique,
                                                           FUN = paste, sep = '-'
                                                           ))), value = TRUE, invert = TRUE)
    if(!all(x[[1]]%in%combination)){
      stop("Choosen contrast is not a possible combination in target file", i ,
           ",please choose the possible combination and try again")
      }
    x.lumi <- lumiR(ilmNames[[i]])
    lumiExpr <- lumiExpresso(x.lumi, bg.correct = TRUE, normalize = TRUE, verbose = TRUE)
    dataMatrix <- exprs(lumiExpr)
    presentCount <- detectionCall(lumiExpr)
    selDataMatrix <- dataMatrix[presentCount > 0,]
    probeList <- rownames(selDataMatrix)
    design <- model.matrix(~0 + f)
    colnames(design) <- unique ## c("C","TT","D","DT")
    fit <- lmFit(selDataMatrix, design)
    cat("fitting the model done...\n")
    contrast.matrix <- makeContrasts(contrasts = x[[1]], levels = design)
    cat("contrast matrix construction done:...\n")
    fit2 <- contrasts.fit(fit, contrast.matrix)
    fit2 <- eBayes(fit2)
    cat("Design contrast done...\n")
    for(i in 1:ncol(contrast.matrix))
      {
      top <- topTable(fit2, coef = i, n = Inf,
                      p.value = pval, lfc = fc, adjust.method = FDRMethod)
      print(top)
      #top <- topTable(fit2, coef=i,number=geneNum, adjust.method="BH")
      if(nrow(top) > 0)
        {
        setnames(top, c("IlogFC", "AveExpr", "t-statistic", "PVALUE", "AdjustedP-val",
                        "B-statistic"))
        print(top)
        if(all(AdditionalInfo==FALSE)){
          foldPval <- top[, c(1, 2, 4)]
          foldPval <- add_rownames(foldPval)
          names(foldPval)[1] <-"ID"
          #setnames(foldPval,c("ID", "IlogFC", "AvgExpr", "PVALUE"))
          }
        else{
          foldPval <- top[, c(1, 2, 4, (which(AdditionalInfo[1] == TRUE)+2),
                              which(AdditionalInfo[2] == TRUE)+4,
                              which(AdditionalInfo[3]==TRUE)+5)]
          cat("checking the extracted columns\n")
          foldPval <- add_rownames(foldPval)
          names(foldPval)[1] <- "ID"
          }
        ilmndiff[[length(ilmndiff) + 1]] <- foldPval
        cat("Illumina differential expression data done...\n")
      }
    } 
    }
  if (length(ilmndiff) == 0){
    stop("No gene found significant at the threshold p-value= ",
         pval , " and foldchange = ", fc, "
         please change your pvalue and 
         foldchange and try again")
    }
  if (length(ilmndiff) < 1){
    stop("please provide atleast one studies with 
         different experimental conditions
         to perform meta-analysis of 
         differential expression")
    }
  cat("ilmdiff is:\n")
  print(ilmndiff)
  return(ilmndiff)
  }



## Function to annotate for fold change illumina data, linked to above code
foldchange2 <- list()
Illumina_combined_fc <- function(df, df1, illumina_diff_range){
  
  if(!is.null(illumina_diff_range)){
    vec <- as.numeric(unlist(strsplit(as.character(illumina_range), split = "," )))
    print(vec)
  }else{
    vec <- length(df)
  }
  ## Through an error if length of range is greater than the
  ## length of annotation files.
  if(length(vec) > length(df1)) 
    stop("Range argument not a multiple of length of annotation")
  
  if(length(vec) < length(df1))
    stop("length of range is not equal to length of annotation ")
  
  if(any(vec== 0))
    stop("Range element is zero")
  
  if(sum(vec) != length(df))
    stop("length of range is not equal to length of input files")
  
  cat("testing more than one gpls ....\n")
  
  illumina_diff_combined <- do.call(c, 
                               Map(function(x,y) lapply(x, function(dat) 
                                 merge(y, dat, by = c("ID"), all.y = TRUE)), 
                                 split(df, rep(seq_along(vec), vec)),
                                 df1))
  
  for(i in 1:length(illumina_diff_combined))
    {
    incProgress(0.1, detail = paste("File", i))
    afterAnot <- na.omit(illumina_diff_combined[[i]][, -c(1, 2)])
    dt <- data.table(afterAnot)
    setkey(dt, PVALUE) ; indx <- dt[,.I[1L], by = SYMBOL]$V1 
    dataWoDuplic <- dt[indx]
    foldchange1[[length(foldchange1) + 1]] <- dataWoDuplic
    }
  foldchange <- lapply(seq_along(foldchange1),
                       function(i) setnames(foldchange1[[i]],
                                            2:ncol(foldchange1[[i]]),
                                            paste0(names(foldchange1[[i]])[-1],i)))
  foldchange<- Reduce(function(x, y) merge(x, y, by = c("SYMBOL"), all = TRUE), foldchange)
  return(foldchange)       
}
########################################## End ########################################


## Illumina function to split the list of files as Names and Targets
names <- list()
targets <- list()
ilmTest <- function(fileList){
  for(i in 1:length(fileList))
    {
    lines <- readLines(fileList[i], n = 20)
    file <- read.table(text = sub('.*(#.*)', '\\1', lines),
                       check.names = FALSE, stringsAsFactors = FALSE,
                       comment.char = '#', fill = TRUE, header = TRUE)
    #print(head(file))
    cols1 <- colnames(file)
    fileTst <- grepl("bead", cols1, ignore.case = TRUE)
    fileTst1 <- grepl("AVG_Signal", cols1, ignore.case = TRUE)
    cat("Testing for target file :\n")
    if(any(fileTst, fileTst1)){
      names[[length(names) + 1]] <- fileList[i]
      }else{
        targets[[length(targets) + 1]] <- fileList[i]
      }
    }
  newList <- list("names" = names, "targets" = targets)
  return(newList)
  }
############################ End #####################################



## Function for codelink data to reorder the dataframes(Automatic 
## Detection of data platform)
fun1 <- function(df){
  #print(df)
  if(ncol(df) >= 2){
    if(ncol(df) == 4)
      {
      df <- df[, -1]
      cat("after removing columns\n")
      }
    indx <- sapply(df, is.character)
    df[indx] <- lapply(df[indx], type.convert)
    means <- colMeans(df)
    df <- df[order(means, decreasing = TRUE)]
    cat("data frame ordering done ...\n")
    if(ncol(df) == 3){
      setnames(df, c("ID", "SIGNALINTENSITY", "SNR"))
      }
    if(ncol(df) == 2){
      setnames(df, c("ID", "SNR"))
      }
    cat("names setting done ... \n")
    df
  }
  }
############################### END ########################################


## Function for Affymetrix data to reorder(Automatic Detection of data platform)
fun2 <- function(data){
  if(ncol(data) > 4){
    stop("Data is not of correct formate, please check the files and try again")
    }
  dim(data)
  indx <- sapply(data, is.character)
  data[indx] <- lapply(data[indx], type.convert)
  d1 <- data[1, , drop = FALSE]
  nums <- d1[, nn <- sapply(d1, is.numeric), drop = FALSE]
  if(length(nums) == 2){
    p <- names(nums[, nums < 1, drop = FALSE])
    val <- setdiff(names(nums), p)
    }
  ch <- d1[, !nn, drop = FALSE]
  id <- names(ch[, grepl('_at$', as.character(unlist(ch))), drop = FALSE])
  if(ncol(data) == 4){
    abs <- setdiff(names(ch), id)
    d <- data[, c(id, val, abs, p)] 
    setnames(d,c("ID", "INTENSITY", "CALL", "PVALUE"))
    }
  else if(ncol(data) == 3)
    {
    if(length(nums) ==2 )
      {
      p <- names(nums[, nums < 1, drop = FALSE])
      d <- data[, c(id, val, p)]
      setnames(d,c("ID", "INTENSITY", "PVALUE"))
      }else{
        abs <- setdiff(names(ch), id)
        d <- data[, c(id, names(nums), abs)]
        setnames(d,c("ID", "INTENSITY", "CALL"))
      }
    }
  else if(ncol(data) == 2){
    p <- names(nums)
    d <- data[, c(id, p)]
    setnames(d, c("ID", "PVALUE"))
    }
  return(d)
  }
##################################### END #######################################


## Function to order illumina data
ilmnOrder <- function(data){
  dim(data)
  indx <- sapply(data, is.character)
  data[indx] <- lapply(data[indx], type.convert)
  d1 <- data[1, , drop = FALSE]
  nums <- d1[, nn <- sapply(d1, is.numeric), drop = FALSE]
  if(length(nums) == 2){
    p <- names(nums[, nums <= 1, drop = FALSE])
    if(length(p) == 0){
      stop("Data is not of illumina formate, please check the data files and try again")
      }
    val <- setdiff(names(nums), p)
    ch <- d1[, !nn, drop = FALSE]
    id <- names(ch)
    d <- data[, c(id, val, p)]
    setnames(d, c("ID", "INTENSITY", "PVALUE"))
    }else{
      ch <- d1[, !nn, drop = FALSE]
      id <- names(ch)
      if(nums > 1){
        stop("Data is not of correct formate, please check the data files and try again")
        }
      p <- names(nums)
      d <- data[, c(id, p)]
      setnames(d, c("ID", "PVALUE"))
      }
  return(d)
}
#################################### End ##############################


## Function to read the processed data to store into a list
ready_data <- list()
proc_data <- function(path, leng){
  for (i in 1 : leng)
    {  
    incProgress(0.1, detail = paste("File",i)) 
    process <- fread(path[[i]], data.table = FALSE, stringsAsFactors = FALSE)
    ready_data[[length(ready_data) + 1]] = process
    }
  print(length(ready_data))
  return(ready_data)
  }
############################# End ######################################

 
shinyServer(function(input, output,session) {
  filenames <- list.files(path = "data", pattern="\\.txt$")
  names(filenames) <- gsub(pattern = "\\.txt$", "", filenames)
  
  
  ## Dropdown box for chosing and loading annotation file
  output$AnnotationFile = renderUI({
    if(input$checkbox == T && input$radio != 2 | input$checkbox != T){
      wellPanel(
        h5("Upload Annotation File"),
        selectInput('dataset', "Choose platform annotation file",
                    c("Please select a file" = '', filenames), multiple = TRUE))
      }
    })
  
  
  output$annotation_range = renderUI({
    if(length(input$dataset) > 1){
      wellPanel(
        textInput("Range", label = "Range")
      )
      }
  })
  
  
  observe({
    print(input$Range)
    print(str(input$Range))
    rang <- as.numeric(unlist(strsplit(as.character(input$Range), split = "," )))
    print(rang)
    print(str(rang))
    print(length(input$dataset))
    print(input$dataset)
  })
  
  
  
  ## Reset input files
  observeEvent(input$Reset_Input,{
    reset("files")
    print(input$num)
    cat("checkbox.............\n")
    print(input$checkbox1)
    })
  
  
  
  ## Reset submit button and annotation files loaded
  observe({
    input$Reset_Input
    #reset("Submit")
    updateActionButton(session, "Submit", label = "Submit")
    updateSelectInput(session, "dataset", selected = FALSE)
    })
  
  
  
  ## Reset other parameters ***************************************************
  observe({
    input$Reset_Input
    updateNumericInput(session, "num", value = "")
    updateNumericInput(session, "num1", value = "")
    updateTextInput(session,"text", value = "Normal-Prostate")
    updateTextInput(session,"text1", value = "spermatogenesis2-spermatogenesis5")
    updateSelectizeInput(session, "select", selected = "Choose DB")
    updateNumericInput(session, "num3", value = "")
    updateNumericInput(session, "num4", value = "")
    updateNumericInput(session, "num5", value = "")
    updateNumericInput(session, "num6", value = "")
    updateRadioButtons(session, "radio", selected = 1)
    updateTabsetPanel(session, "MamgedTabs", selected = paste0("panel"))
    updateCheckboxInput(session, "checkbox1", value = "")
    updateCheckboxInput(session, "checkbox2", value = "")
    updateCheckboxInput(session, "checkbox3", value = "")
    updateRadioButtons(session, "radio1", selected = 'BH')
    })
  
  
  
  ## Display message on submit without loading input files ********************
  observeEvent(input$Submit,{
    if(is.null(input$files)){
      session$sendCustomMessage(type = 'testmessage',
                                message = 'No files to submit')
      }
    })
  
  
  ## Display message on reset ************************************************* 
  observeEvent(input$Reset_Input, {
    if(is.null(input$files) && (input$num == "" || is.null(input$num) ||is.na(input$num))
       && (input$num1 == "" || is.null(input$num1) ||is.na(input$num1))
       && (input$num3 == "" || is.null(input$num3) ||is.na(input$num3))
       && (input$num4 == "" || is.null(input$num4) ||is.na(input$num4))
       && (input$num5 == "" || is.null(input$num5) ||is.na(input$num5))
       && (input$num6 == "" || is.null(input$num6) ||is.na(input$num6))
       && (input$checkbox1 == "" || is.null(input$checkbox1) || !input$checkbox1)
       && (input$checkbox2 == "" || is.null(input$checkbox2) || !input$checkbox2)
       && (input$checkbox3 == "" || is.null(input$checkbox3) || !input$checkbox3)){
      session$sendCustomMessage(type = 'testmessage',
                                message = 'Nothing to reset')
      }else{
      session$sendCustomMessage(type = 'testmessage',
                              message = 'All parameters are reset')
    }
  })
  
  
  ## Display Differential Expression Parameters *******************************
  output$hide=renderUI({ 
    if(input$checkbox == T && input$radio == 1){
      wellPanel(
        h5("Affymetrix Differential Expression"),
        list(numericInput("num", label = h5("Fold Change"),
                          value = NULL, min = 0.5, max = 5, step = 0.5),
             numericInput("num1", label = h5("P-Value"),
                          value = NULL, min = 0.01, max = 0.05, step = 0.01),
             textInput("text", label = h5("Make Contrasts"),
                       value = "Normal-Prostate")))
      }else if(input$checkbox == T && input$radio == 2){
        wellPanel(
          h5("Codelink Differential Expression"),
          list(
            textInput("text1", label = h5("Make Contrasts"), 
                      value = "spermatogenesis2-spermatogenesis5"),
            selectizeInput("select", label = h5("Choose Annotation"),
                           choices = list("h10kcod.db" = "h10kcod",
                                          "h20kcod.db" = "h20kcod",
                                          "hwgcod.db" = "hwgcod"),
                           options = list(placeholder = 'Choose DB',
                                        onInitialize = I('function()
                                                         { this.setValue(""); }'))),
            numericInput("num3", label = h5("Fold Change"),
                         value = NULL, min = 0.5, max = 5, step = 0.5),
            numericInput("num4", label = h5("P-Value"),
                         value = NULL, min = 0.01, max = 0.05, step = 0.01)))
        }else if (input$checkbox == T && input$radio == 3){
          wellPanel(
            h5("Illumina Differential Expression"),
            list(
              textInput("text3", label = h5("Make Contrasts"), value = "Enter conditions.."),
              numericInput("num5", label = h5("Fold Change"),
                           value = NULL, min = 0.5, max = 5, step = 0.5),
              numericInput("num6", label = h5("P-Value"),
                           value = NULL, min = 0.01, max = 0.05, step = 0.01)))
          }
    })
  
  
  
  ## Display options for additional information for affymetrix data ***********
  output$FDRAffymetrix <- renderUI({
    if(input$checkbox == T && (input$radio == 1 | input$radio == 3)){
      wellPanel(
        h5("Additional information you want in supplementary file"),
        list(
          checkboxInput("checkbox1", label = "t-statistic", value = FALSE),
          checkboxInput("checkbox2", label = "adjusted p-value", value = FALSE),
          checkboxInput("checkbox3", label = "B-statistic", value = FALSE),
          actionButton("selectAll", label = "Select All"),
          actionButton("deselectAll", label = "Deselect All"),
          br(),
          br(),
          #checkboxInput("checkbox4", label = "All", value = FALSE),
          radioButtons("radio1", label = h5("Choose False Discovery Control"),
                       choices = list("Benjamini & Hochberg (BH)" = "BH",
                                      "Benjamini & Yekutieli (BY)" = "BY",
                                      "holm" = "holm",
                                      "hochberg" = "hochberg",
                                      "hommel" = "hommel",
                                      "bonferroni" = "bonferroni",
                                      "fdr" = "fdr",
                                      "none" = "none"),
                       selected = NULL)))
      }
    })
  
  
  
  ## Display options for additional information for codelink data *************
  output$FDRCodlink <- renderUI({
    if(input$checkbox == T && input$radio == 2){
      wellPanel(
        h5("Additional information you want in supplementary file"),
        list(
          checkboxInput("checkbox4", label = "meanSNR", value= FALSE),
          checkboxInput("checkbox5", label = "t-statistic", value = FALSE),
          checkboxInput("checkbox6", label = "adjusted p-value", value = FALSE),
          checkboxInput("checkbox7", label = "B-statistic", value = FALSE),
          actionButton("selectAll1", label = "Select All"),
          actionButton("deselectAll1", label = "Deselect All"),
          br(),
          br(),
          #checkboxInput("checkbox5", label = "All", value = FALSE),
          radioButtons("radio1", label = h5("Choose False Discovery Control"),
                       choices = list("Benjamini & Hochberg (BH)" = "BH",
                                      "Benjamini & Yekutieli (BY)" = "BY",
                                      "holm" = "holm",
                                      "hochberg" = "hochberg",
                                      "hommel" = "hommel",
                                      "bonferroni" = "bonferroni", 
                                      "fdr" = "fdr",
                                      "none" = "none"),
                       selected = NULL)))
      }
    })
  
  
  
  
  #* For Supplementary information for affymetrix and illumina data platforms
  #* This observer will update checkboxes 1 - 4 to TRUE whenever selectAll is clicked 
  observeEvent(
    eventExpr = input$selectAll,
    handlerExpr = 
      {
        lapply(paste0("checkbox", 1:4),
               function(x)
                 {
                 updateCheckboxInput(session = session, 
                                     inputId = x, 
                                     value = TRUE)
                 }
        )
        }
    )
  
  
  
  #* This observer will update checkboxes 1 - 4 to FALSE whenever deselectAll is clicked
  observeEvent(
    eventExpr = input$deselectAll,
    handlerExpr = 
      {
        lapply(paste0("checkbox", 1:4),
               function(x)
                 {
                 updateCheckboxInput(session = session,
                                     inputId = x, 
                                     value = FALSE)
                 }
        )
        }
    )
  
  
  #* For supplementary information choosen by user for codelink data platform.
  #* This observer will update checkboxes 1 - 5 to TRUE whenever selectAll is clicked 
  observeEvent(
    eventExpr = input$selectAll1,
    handlerExpr = 
    {
      lapply(paste0("checkbox", 4:7),
             function(x)
             {
               updateCheckboxInput(session = session, 
                                   inputId = x, 
                                   value = TRUE)
             }
      )
    }
  )
  
  
  
  #* This observer will update checkboxes 1 - 5 to FALSE whenever deselectAll is clicked
  observeEvent(
    eventExpr = input$deselectAll1,
    handlerExpr = 
    {
      lapply(paste0("checkbox", 4:7),
             function(x)
             {
               updateCheckboxInput(session = session, 
                                   inputId = x, 
                                   value = FALSE)
             }
      )
    }
  )
  
  
  
  
  ## Ths reactive function accepts input files and stores in a variable.
  infiles <- eventReactive(input$Submit,{ 
    if (is.null(input$files)){
      return(NULL)
      }
    input$files
    })
  
  
  
  ## This renderDataTable displays list of files uploaded
  output$file <- renderDataTable({
    if (is.null(input$files)){
      return(NULL)
      }
    infiles()
    })
  
  
  
  ## This reactive function reads the file and process the uploaded data files
  list_files <- eventReactive(input$Submit,{
    fileList <-  input$files[['name']]  
    path1 <- input$files[['datapath']][1]
    lines <- readLines(path1, n = 10)
    IlluminaCheck <- read.table(text=sub('.*(#.*)', '\\1', lines),
                        check.names = FALSE, stringsAsFactors = FALSE, 
                        comment.char = '#', header = TRUE,
                        fill = TRUE, row.names = NULL)
    
    #* Check for illumina data 
    na.omit(IlluminaCheck)
    names <- names(IlluminaCheck)
    IlluminaCols <- ncol(IlluminaCheck)
    n <- grepl("bead", names, ignore.case = TRUE)
    IsAvgSignal <- grepl("AVG_Signal", names, ignore.case = TRUE)
    IsDetectionPval <- grep('^(detection|Pval)', names, ignore.case = TRUE)
    CheckLength <- all(length(names) > 4, length(IsDetectionPval) >= 2)
    cat("Is this a bead array:", any(n, IsAvgSignal, CheckLength), "\n")
    an <- any(n, IsAvgSignal, CheckLength)
    
    ##* Checking if files are CEL 
    CheckforCEL <- grepl("*.cel$", fileList,ignore.case = TRUE)
  
    ##* Checking if files are text
    CheckforTXT <- grepl("*.TXT$", fileList)
  
    ##* Checking if files are both CEL and text
    CheckforCEL_TXT <- all(sapply(c('\\.CEL$', '\\.txt$'),
                        function(x) any(grepl(x, fileList, ignore.case = TRUE))))
    
    ##* check for if the files belong to DE for codelink
    TextforCodelink <- c('TXT', 'txt')
    CheckforCodelinkDE <- all(TextforCodelink %in% sub('.*\\.', '', fileList) )
     
  
    ##* checking if the files belong to DE of illumina bead array
    CheckforIlluminaDE <- grepl("*.txt$", fileList)
    chkbox <- input$checkbox
  
    ##* To check if the checkbox is needed for differential expression
    if(any(CheckforCEL_TXT, CheckforCodelinkDE)){
      if(!(chkbox)){
        stop("You seem to be interested in differential expression,
             please make use of differential expression
             checkbox and try again")
      }
      }
  
    
  ## Test to check if check box is needed for illumina data. If the user loads differential
  ## expression data but forgot to check the differential expression box, the program will 
  ## try to identify the differential data here by dividing the list into names and targets
  if(all(CheckforIlluminaDE)){
    ilmpath = input$files[['datapath']]
    nameLst <- ilmTest(ilmpath)
    
    ##* creating seperate list for data files and phenotype files
    ilmNames <- nameLst$names
    ilmTargets <- nameLst$targets
    if(length(ilmNames) > 0 && length(ilmNames) < length(ilmpath)){
      if(!chkbox){ 
        stop("You seem to be interested in differential expression,
             please make use of differential expression
             checkbox and try again")
      }
    }  
    }
    cat("Test for affymetrix DE:", CheckforCEL_TXT, "\n")
    cat("Test for codelink DE:", CheckforCodelinkDE, "\n")
    cat("Test for illumina DE:", all(CheckforIlluminaDE, chkbox), "\n")
    withProgress(message = 'Processing Data . . .',{
      if(all(CheckforCEL) && !chkbox)
        {
        if(input$radio == 2){
          stop("Data files uploaded are not of codelink format,
               please check the data files and try again")
          }
        if(input$radio == 3){ 
          stop("Data files uploaded are not of illumina format,
               please check the data files and try again")
          }
        genAffyData = ReadAffy(filenames = input$files[['datapath']])
        cat("Affymetrix Raw data selected\n")
        print(genAffyData)
        df <- processed_data(genAffyData)
        #print(head(df))
        df
        }
      else if(all(CheckforTXT) && !chkbox)
        {
        if(input$radio == 1){
          stop("Data uploaded is not of Affymetrix format,
               please check the data files and try again")
          }
        if(input$radio == 3){ 
          stop("Data files uploaded is not of Illumina format,
               please check the data files and try again")
          }
        cat("Codelink raw data selected\n")
        codelinkpath = input$files[['datapath']]
        code <- codelink_data(codelinkpath)
        code <- lapply(code, na.omit)
        }
      else if(CheckforCEL_TXT && chkbox)
        {
        if(input$radio == 2){
          stop("Data uploaded is not of Codelink format,
               please check the data files and try again")
          }
        if(input$radio == 3){
          stop("Data files uploaded is not of illumina format,
               please check the data files and try again")
          }
        cat("Checking for making contrasts of affymetrix data:")
        fold <- input$num
        if(is.na(fold)){
          stop("Foldchange is missing, please choose and try again")
          }
        pval <- input$num1
        if(is.na(pval)){
          stop("p-value is missing, please choose and try again")
          }
        conditions <- input$text
        if(is.na(conditions)){
          stop("Comparision information is missing, please provide and try again")
          }
        conditions <- strsplit(conditions, split = ",")
        cat("Disease conditions choosen to compare:\n")
        cat("Differential expression choosen\n")
        
        ## Renaming the files, as it gets changed on uploading the files
        from <- input$files[['datapath']]
        newnames <- file.path(dirname(from), basename(input$files[['name']]))
        newone <- newnames[1]
        
        ##Removing the file name to set the path
        newpath <- sub("/[^/]*$", "", newone) 
        n <- file.rename(from, newnames)
        
        ## Checking for additional information and FDR method choosen
        AdditionalInfo <- c(input$checkbox1,input$checkbox2,input$checkbox3)
        cat("checking for additional information")
        FDRMethod <- input$radio1
        
        ## Calling affymetrix differential expression function
        affy_fc <- fold_change_affy(fileList, newpath, df1, fold, pval, newnames,
                           conditions, AdditionalInfo, FDRMethod)
        }
      else if(CheckforCodelinkDE && chkbox)
        {
        if(input$radio == 1){
          stop("Data uploaded is not of Affymetrix format,
               please check the data files and try again")
          }
        if(input$radio == 3){
          stop("Data uploaded is not of illumina format,
               please check the data files and try again")
          }
        cat("Checking for making contrasts of codelink data:")
        cntrst <- input$text1
        if(is.na(cntrst)){
          stop("Comparison information is missing, please provide and try again")
          }
        fc <- input$num3
        if(is.na(fc)){
          stop("Foldchange is missing, please provide and try again")
          }
        pval <- input$num4
        if(is.na(pval)){
          stop("p-value is missing, please choose and try again")
          }
        cntrsts <- strsplit(cntrst, split = ",")
        cat("checking for annotation DB:")
        annotation <- input$select
        if(is.na(annotation)){
          stop("Annotation data base is missing, please choose and try again")
          }
        from <- input$files[['datapath']]
        cat("From datapath:", from, "\n")
        newnames <- file.path(dirname(from), basename(input$files[['name']]))
        cat("New names:",newnames,"\n")
        newone <- newnames[1]
        
        ##Removing the file name to set the path
        newpath <- sub("/[^/]*$", "", newone) 
        cat("New path:",newpath,"\n")
        n <- file.rename(from, newnames)
        
        ## Checking for additional information and FDR method choosen
        AdditionalInfo <- c(input$checkbox4,input$checkbox5,
                            input$checkbox6,input$checkbox7)
        FDRMethod <- input$radio1
        codelink <- codelinkdiff(fileList, newpath, cntrsts, pval, fc, annotation,
                             newnames, AdditionalInfo, FDRMethod)
        }
      else if(all(CheckforIlluminaDE) && chkbox)
        {
        cat("you are here in illumina")
        if(length(ilmNames) > 0 && length(ilmNames) < length(ilmpath)){
          if(input$radio == 2){
            stop("Data uploaded is not of codelink format,
                 please check the data files and try again")
            }
          if(input$radio == 1){
            stop("Data files uploaded is not of Affymetrix format,
                 please check the data files and try again")
          }
          
          ## To read the phenotype file to make sure the model matrix 
          ##construction is done correctly
          from <- input$files[['datapath']]
          newnames <- file.path(dirname(from), basename(input$files[['name']]))
          x <- input$text3
          cat("contrasts choosen:", x, "\n")
          if(is.na(x)){
            stop("Make contrasts information is missing, please provide and try again")
            }
          cat("Spliting names of list:\n")
          x <- strsplit(x, split = ",")
          fc <- input$num5
          if(is.na(fc)){
            stop("Foldchange information is missing, please choose and try again")
            }
          pval <- input$num6
          if(is.na(pval)){
            stop("pvalue information is missing, please provide and try again")
          }
          
          ## Both the phenotype files and raw data files are of .txt format, so we divide the
          ## list of files into phenotype and raw data files.
          ilmpath = input$files[['datapath']]
          nameLst <- ilmTest(ilmpath)
          
          ## creating seperate list for data files and phenotype files
          ilmNames <- nameLst$names
          ilmTargets <- nameLst$targets
          cat("Name list done:\n")
          cat("Target list done:\n")
          
          ## Checking for user input of additional data and FDR control
          AdditionalInfo <- c(input$checkbox1,input$checkbox2,input$checkbox3,input$checkbox4)
          cat("checking for additional information")
          FDRMethod <- input$radio1
          illumina <- illuminaDe(ilmNames, ilmTargets, x, fc, pval, newnames,
                                 AdditionalInfo, FDRMethod )
          cat("Illumina df's for DE:\n")
          }
        else if(length(ilmNames) == length(ilmpath)){
          stop("No target files found, please check the files and try again")
          }
        else if(length(ilmNames) == 0){ 
          stop("No file with bead information found,
               please check your data files and try again")
        }
        }
      else if(any(n, IsAvgSignal, CheckLength)) 
        {
        if(chkbox){         
          stop("No file compartible for differential expression,
               please check the data file and try again")
          }
        if(input$radio == 1){
          stop("Data uploaded is not of Affymetrix format,
               please check the data files and try again")
          }
        if(input$radio == 2){
          stop("Data uploaded is not of Codelink format,
               please check the data files and try again")
        }
        ## Illumina raw data processing
        cat("Illumina raw data processing:\n")
        illumina <- input$files[['datapath']]
        illuminaDat <- illuminaRaw(illumina)
        }
      else if(all(CheckforCEL)|all(CheckforTXT) && chkbox){
        stop("No target files found, please check target files and try again")
        }
      else if(all(CheckforIlluminaDE) && chkbox){
        stop("No data files found, please check data files and try again")
        }else{
          cat("you are here\n")
          if(chkbox){
            stop("Are you sure of perfoming differential expression
                 on the loaded data files,as the files are not
                 compatible for differential expression, 
                 please check your data and try again")
            }
          cat("Processed data  choosen\n")
          path <- input$files[['datapath']]
          leng <- length(path)
          
          ## Reading all the processed files and storing in a list
          procs <- proc_data(path, leng)
          procs <- lapply(procs, na.omit)
          vec <- procs[[1]][1, ]
          aff <- sapply(vec, function(x) grepl("_at$", x))
          ilm <- sapply(vec, function(x) grepl("ILMN_", x))
          if('TRUE'%in%aff)
            {
            if(input$radio == 2){ 
              stop("Data uploaded is not of codelink format,
                   please check the data files and try again")
              }
            if(input$radio == 3){
              stop("Data files uploaded is not of Illumina format,
                   please check the data files and try again")
              }
            cat("Affymetrix data choosen: \n")
            procss <- lapply(procs, fun2)
            }
          else if(!any(aff, ilm))##input$radio==2)
            {
            if(input$radio == 1){
              stop("Data uploaded is not of Affymetrix format,
                   please check the data files and try again")
              }
            if(input$radio == 3){ 
              stop("Data uploaded is not of Illumina format,
                   please check the data files and try again")
              } 
            cat("Codelink data choosen\n")
            procss <- lapply(procs, fun1)
            if(length(procss) == 0){
              stop("Data files seem not of correct formate,
                   please check the data files and try again")
              }
            if(length(procss) < length(procs)){
              warning("Some file are skipped")
              }
            procss <- Filter(function(x) !is.null(x), procss)
            procss <- lapply(procss, na.omit)
            } 
          else if('TRUE'%in%ilm)##input$radio==3)
            {
            if(input$radio == 2){
              stop("Data uploaded is not of codelink format,
                   please check the data files and try again")
              }
            if(input$radio == 1){ 
              stop("Data uploaded is not of Affymetrix format,
                   please check the data files and try again")
              }
            cat("Illumina data choosen\n")
            procs <- lapply(procs, ilmnOrder)
            procs <- lapply(procs, na.omit)
          }
          }
      })
    })
  
  
  ## renderDataTable function to display content of one file ******************
  output$sourced <- renderDataTable({ 
    input$Submit
    if (is.null(input$files)){
      return(NULL)
      }
    first_file <- list_files()[[1]]
    })
  
  
  ## Reactive to select a file from the data-base (Annotation File)
  annotation_datafiles_list <- eventReactive(input$Submit,{
    cat("Testing gpl file.........\n")
    print(input$dataset)
    annotation_file_list <- list()
    annotation <- input$select
    cat("print annotation files selected...\n", input$dataset, "\n")
    lapply(input$dataset, function(i){
      if (grepl("[/\\\\]", i)) {
        stop("Invalid dataset")
      }
    })
    if(input$radio == 2 && input$checkbox == T){
      if(annotation == "h10kcod" |
         annotation == "h20kcod" |
         annotation == "hwgcod"){
        cat("No annotation file to display,bioconductor annotation choosen\n")
        }
      }else{
        annotation_data <- lapply(input$dataset, function(i){
           annotation_file_list[[length(annotation_file_list) + 1]] <- read.csv(file.path("data", i),
                                                                                sep = "\t",comment.char = "#",
                                                                                check.names=FALSE, stringsAsFactors=FALSE)
          })
        cat("testing annotation file data......\n")
        col_names <- c("ID", "NAME", "SYMBOL")
        annotation_data_with_colnames <- lapply(annotation_data, setNames, col_names)
        }
    })
  
  
  
  ## This reactive clean the annotation file *****************************
  cleaned_annotation_files <- eventReactive(input$Submit,{
    junk_free_annotation_files <- list()
    if (is.null(annotation_datafiles_list())){
      return(NULL)
    }
    annotation_files_store <- lapply(annotation_datafiles_list(), function(i){
      df1 <- i
      df1[] <- lapply(i, as.character)
      df1 <- df1[df1[, 3] != '', ]
      df1 <- df1[df1[,3] != '---',]
      df1 <- df1[!grepl('[/]', df1$SYMBOL), ]
      df1 <- df1[!grepl('NM_', df1$SYMBOL), ]
      df1 <- df1[, c(1:3)]
      junk_free_annotation_files[[length(junk_free_annotation_files) + 1]] <- na.omit(df1)
    })
    cat("Testing junk free GPL file done ...\n")
    print(head(annotation_files_store[[1]]))
    annotation_files_store
    }) 
  
  
  
  # This renderDataTable is to display the annotation file *******************
  output$annotation <- renderDataTable({
    withProgress(message = 'Loading Data . . .', {
      if (is.null(input$files)){
        return(NULL) 
      }
      cat("testing annotation display.....\n")
      cleaned_annotation_files()[[1]]
      })
    })
  
  
  
     
  ## This reactive function generated the data by calling various functions
  data <- eventReactive(input$Submit,{
    if (is.null(input$files)){
      return(NULL)
      }
    #if(is.null(input$dataset)){
      #return(NULL)
      #}
    df1 <-  cleaned_annotation_files()
    withProgress(message = 'Annotating . . .',value = 0.01,{
      rout <- input$radio
      df <- list_files()
      print(length(df))
      files <- length(df)
      cat("You are in the main program... \n")
      cat("Sample frame generation done...\n")
      forFold <- grepl('AlogFC', colnames(df[[1]]))[2]
      forFoldIlmn <- grepl('IlogFC', colnames(df[[1]]))[2]
      if(forFold)
        {
        foldchange <- affy_combined_fc(df, df1, input$Range)
        if(nrow(foldchange) == 0){
          stop("Data files uploaded and annotation does not matching,
               please check the files and try again")
          }
        foldchange
        }
      else if(forFoldIlmn)
        {
        cat(" Processing Illumina differential data: ")
        foldchange <- Illumina_combined_fc(df, df1)
        if(nrow(foldchange) == 0){
          stop("Data files uploaded and annotation does not matching,
               please check the files and try again")
          }
        foldchange
        }
      else if(rout == 1 && any(length(df[[1]]) == 4,
                               length(df[[1]]) == 3,
                               length(df[[1]] == 2)
                               )
              )
        {
        cat("Affymetrix  data selected \n")
        affyFin <- affymetrix(df, files, df1, input$Range)
        affyFin <- data.frame(affyFin, stringsAsFactors = FALSE) 
        if(nrow(affyFin) == 0){
          stop("Data files uploaded and annotation does not matching,
               please check the files and try again")
          }
        affyFin 
        }
      else if(any(length(df[[1]]) == 3,length(df[[1]] == 2)) && rout == 2 )
        {
        cat("Codelink data selected \n")
        symbols <-  grepl("SYMBOL", colnames(df[[1]]))
        if(any(symbols)){
          cdlnkFin <- codlinkdiff1(df)
          }else{
            cdlnkFin <- codelink(df, files, df1, input$Range)
            } 
        cdlnkFin <- data.frame(cdlnkFin, stringsAsFactors = FALSE)
        if(nrow(cdlnkFin) == 0){ 
          stop("Data files uploaded and annotation does not matching,
               please check the files and try again")
          }
        cdlnkFin 
        }
      else if(any(length(df[[1]]) == 3,length(df[[1]] == 2)) && rout == 3)
        {
        cat("Illumina data selected \n")
        ilmdata <- illumina(df, files, df1, input$Range)
        ilmdata <- data.frame(ilmdata, stringsAsFactors = FALSE)
        if(nrow(ilmdata) == 0){
          stop("Data files uploaded and annotation does not matching,
               please check the files and try again") 
          }
        ilmdata 
        }
      })
    })
  
  
  
   ## This reactive function to summarize the data *************************** 
  fulldata <- eventReactive(input$Submit,{ 
    withProgress(message = 'Generating Final Summary . . .',value = 0.01,{
      cat ("Data Summary done ...\n")
      affy <- data()
      })
    })
  
  
  
  ## This reactive summarize the results and assign score, and then cummulative score.
  summary <- eventReactive(input$Submit,{
    final <- fulldata()
    if (is.null(final)){
      return(NULL)
      } 
    if(ncol(final)<3){ 
      stop("Less than two data files are selected,
           please provide two or more file and try again")
      }
    logical <- grepl('PVALUE', colnames(final))
    IlluminaCols <- ncol(final)-1
    sm <- sum(logical, na.rm = TRUE)
    logical3 <- grepl('INTENSITY', colnames(final))
    logical4 <- grepl('CALL', colnames(final))
    logical1 <- grepl('SNR', colnames(final))[3]
    logical2 <- grepl('AlogFC|ClogFC|IlogFC', colnames(final))[2]
    chkbox <- input$checkbox
    if(any(logical) && any(logical3) | sm == IlluminaCols)
      {
      cat("Final data generation reached:\n")
      final <- data.frame(final)
      finaldata <- final[grep('^(SYMBOL|PVALUE)', names(final))]
      finaldata[is.na(finaldata)] <- 0
      print(head(finaldata))
      
      ## converting the column factors
      indx <- sapply(finaldata[, -1], is.factor)[1]
      if(indx){
        finaldata[, 2:ncol(finaldata)] = lapply(finaldata[, -1],
                                                function(x)
                                                  as.numeric(levels(x))[x])
        }
      finaldata$CummuScore <- rowSums(ifelse
                                      (finaldata[,-1] == 0, 0,
                                        ifelse
                                        (finaldata[, -1] <= 0.05, 2,
                                          ifelse
                                          (finaldata[,-1] >= 0.065, -2, 0)
                                          )
                                        )
                                      )
      finaldata$CummuScore <- as.numeric(finaldata$CummuScore)
      finaldata = finaldata[ order(finaldata$CummuScore, decreasing = TRUE),]
      finaldata$FinalCall <- ifelse(finaldata$CummuScore > 0, 'P',
                                    ifelse(finaldata$CummuScore < 0, 'A', 0))
      cat("Affymetrix final data score assignment done ...\n")
      finaldata
      }
    else if(!any(logical) && any(logical4))
      {
      cat("Final data has CALL and INTENSITY present\n")
      final <- data.frame(final)
      finaldata <- final[grep('^(SYMBOL|CALL)', names(final))]
      finaldata[] = lapply(finaldata, as.character)
      finaldata[is.na(finaldata)] <- 0
      finaldata$CummuScore <- rowSums(ifelse(finaldata[, -1] == 'P', 2, 
                                             ifelse(finaldata[, -1] == 'A', -2, 0)))
      finaldata$CummuScore <- as.numeric(finaldata$CummuScore)
      #print(head(finaldata,100))
      finaldata = finaldata[ order(finaldata$CummuScore, decreasing = TRUE),]
      finaldata$FinalCall <- ifelse(finaldata$CummuScore > 0, 'P',
                                    ifelse(finaldata$CummuScore < 0, 'A', 0))
      cat("Affymetrix final data score assignment done ...\n")
      finaldata
      }
    else if(logical1)
      {
      finaldata <- final[grep('^(SYMBOL|SNR)', names(final))]
      finaldata[is.na(finaldata)] <- 0
      indx <- sapply(finaldata[, -1], is.factor)[1]
      if(indx){
        finaldata[, 2:ncol(finaldata)] = lapply(finaldata[, -1],
                                                function(x)
                                                  as.numeric(levels(x))[x])
        }
      finaldata$CummuScore <- rowSums(ifelse(finaldata[, -1] == 0, 0, 
                                             ifelse(finaldata[, -1] >= 1, 2, -2)))
      finaldata$CummuScore <- as.numeric(finaldata$CummuScore)
      finaldata = finaldata[ order(finaldata$CummuScore,decreasing = TRUE),]
      finaldata$FinalCall <- ifelse(finaldata$CummuScore > 0, 'P',
                                    ifelse(finaldata$CummuScore < 0, 'A', 0))
      cat("Codelink final data score assignment done ... \n")
      finaldata
      }
    else if(logical2 && chkbox)
      {
      cat("Performing Differential Expression\n")
      chk <- input$checkbox
      fold <- data.frame(final)
      
      ##* Checking if there are many comparisions (comparision for more than one conditions
      ##* that generates many fold changes and other attributes)
      howmanyFC <- grep("logFC", names(fold))
      if(length(howmanyFC) == 1){
        finaldata <- fold[, c(1,2,4)] 
        finaldata[is.na(finaldata)] <- 0 
        cat("one comparision condition choosen")
        finaldata$CummuScore <- ifelse(finaldata[, 2] == 0, 0, 
                                     ifelse(finaldata[, 2] > 0, 2, -2))
        }else{
          finaldata <- fold[grep('^(SYMBOL|AlogFC|ClogFC|IlogFC|PVALUE)', names(fold))]
          #finaldata <- fold[,names(fold)[col.names]]
          #finaldata <- fold[, c(1, 2, 4, seq(5, ncol(fold), 3))]
          finaldata[is.na(finaldata)] <- 0
          finaldata$CummuScore <- rowSums(ifelse(finaldata[, grep("logFC",names(finaldata))] == 0, 0,
                                                 ifelse(finaldata[, grep("logFC",names(finaldata))] > 0, 2, -2)
                                                 )
                                          )
        }
      finaldata$CummuScore <- as.numeric(finaldata$CummuScore)
      finaldata = finaldata[order(finaldata$CummuScore, decreasing = TRUE), ] 
      finaldata$FinalCall <- ifelse(finaldata$CummuScore > 0, 'P',
                                    ifelse(finaldata$CummuScore < 0, 'A', 0))
      finaldata
      }
    else if(input$radio == 3)
      {
      
      # Extracting pvalue columns only
      finaldata <- final[, c(1, seq(3, ncol(final), 2))]
      finaldata[is.na(finaldata)] <- 0
      
      ## Assigning values and summing up.
      finaldata$CummuScore <- rowSums(ifelse(finaldata[, -1] == 0, 0, 
                                           ifelse(finaldata[, -1] <= 0.05, 2, 
                                                  ifelse(finaldata[,- 1] >= 0.065, -2, 0)
                                                  )
                                           )
                                      )
      ## Sorting Cummulative score in decreasing order
      finaldata$CummuScore <- as.numeric(finaldata$CummuScore)
      finaldata = finaldata[order(finaldata$CummuScore,decreasing = TRUE),] 
      finaldata$FinalCall <- ifelse(finaldata$CummuScore > 0, 'P',
                                  ifelse(finaldata$CummuScore < 0, 'A', 0)
      )
      }
    finaldata
    })
    
  
  
  ## This displays the final data as an autput
  observeEvent(input$Submit,{
    output$final <- renderDataTable({
      if (is.null(data())){
        return(NULL) 
      }
      print(summary())
      summary()
      })
    })
  
  
  
  ## This displays data in table with full information including p-value and signal intensity
  output$full <-renderDataTable({
    if (is.null(fulldata())) {
      return(NULL)
      }
    alldata <- fulldata()
    print(head(alldata, 20))
    logical <- grepl('INTENSITY', colnames(alldata))
    logical1 <- grepl("CALL", colnames(alldata))
    logical2 <- grepl("PVALUE", colnames(alldata))
    if(any(logical1) && !any(logical2)) {
      finaldata <- alldata[grep('^(SYMBOL|INTENSITY|CALL)', names(alldata))]
      }
    else if(any(logical)) {
      cat("U R almost done...")
      finaldata <- alldata[grep('^(SYMBOL|INTENSITY|SIGNALINTENSITY|SNR|PVALUE)',
                                names(alldata))]
      }else{
        alldata
        }
    })
  
  
  
  ## This display the data with fold change
  output$fold <- renderDataTable({
    if (is.null(summary1())){
      return(NULL)
      }
    summary1()  
    })
  
  
  
  ## This download handler is to download Affymetrix Sample Data
  output$AffymetrixData <- downloadHandler( 
    filename <- "Affymetrix Data.zip",
    content <- function(file){
      file.copy("www/Affymetrix Data.zip", file)
      })
  outputOptions(output, "AffymetrixData", suspendWhenHidden=FALSE)
  
  
  
  ## This download handler is to download Codelink Sample Data
  output$CodelinkData <- downloadHandler( 
    filename <- "Codelink Data.zip",
    content <- function(file){
      file.copy("www/Codelink Data.zip", file)
      })
  outputOptions(output, "CodelinkData", suspendWhenHidden=FALSE)
  
  
  
  ## This download handler is to download Illumina Sample Data
  output$IlluminaData <- downloadHandler( 
    filename <- "Illumina Data.zip",
    content <- function(file){
      file.copy("www/Illumina Data.zip", file)
      })
  outputOptions(output, "IlluminaData", suspendWhenHidden=FALSE)
  
  
  
  ## This download handler is to download  summary data
  output$downloadData1 <- downloadHandler(
    filename <- function() { paste('Summary', '.csv', sep = '')},
    content <- function(file) {
    write.csv(summary(), file)
    })
  
  
  ## This download handler is to download Supplementary file
  output$downloadData2 <- downloadHandler(
    filename <- function() { paste('Supplementary', '.csv', sep = '')},
    content <- function(file){
    write.csv(fulldata(), file)
    })
})

##################### End ###################################