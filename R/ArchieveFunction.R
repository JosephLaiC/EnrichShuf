#' Input the bed format data.frame and convert to grange
#'
#' @param data input the data.frame or the path to bed file.
#' @param name assign the fourth column name and keep it exsist in grange
#' @param strand if strand assign as TRUE means input data contain the strand information at column 6, if assign NULL strand will be *, otherwise could be "+" or "-"
#'
#'  
bedfromfile <- function(data, name="name", strand=NULL){

  #library(GenomicFeatures)

  column.names <- c("chr", "start", "end", name)

  if (is.data.frame(data)){
    table <- data
  } else if (file.exists(data)){
    table <- readr::read_tsv(data, col_names=FALSE, show_col_types=FALSE)
  } else {
    stop("Check the input data format or check the file exsist in path")
  }

  if (ncol(table) < 4){
    stop("Check the input bed table or file conatained reegion name information")
  }

  if (isTRUE(strand)){

    if (ncol(table) < 6){
      stop("Check the input bed table or file conatained strand information")
    }

    table     <- table[,1:6]
    table[,5] <- NULL
    colnames(table) <- c(column.names, "strand")
  } else {
    table     <- table[,1:4]
    colnames(table) <- column.names
  }

  grange <- GenomicRanges::makeGRangesFromDataFrame(
    table, keep.extra.columns=TRUE, starts.in.df.are.0based=TRUE)

  if (any(isTRUE(strand), is.null(strand))){

    return(grange)

  } else if (strand%in%c("+", "-")){

    GenomicRanges::strand(grange) <- strand
    return(grange)

  } else {
    stop("Check strand must be TRUE, NULL, + or -.")
  }

}

#' Convert bed format data.frame or grange to txdb object
#'
#' @param data input data
#' @param name assign the fourth column name and keep it exsist in grange
#' @param strand if strand assign as TRUE means input data contain the strand information at column 6, if assign NULL strand will be *, otherwise could be "+" or "-"
#'
#'  
txdbfromBed <- function(data, name="name", strand=NULL){
  #library(dplyr)

  if (class(data)[[1]]=="GRanges"){

    if (name%in%colnames(data.frame(data))){
      grange <- data
    } else {
      stop("Input data format is grange and please ckeck it contains, the name you assigned..")
    }

  } else {
    grange <- bedfromfile(data, name=name, strand=strand)
  }

  grange$type          <- c("transcript")
  grange$transcript_id <- data.frame(grange)[,name]
  grange$gene_id       <- data.frame(grange)[,name]

  txdb <- GenomicFeatures::makeTxDbFromGRanges(grange)
  return(txdb)

}


#' Output the distance of factor X to nearest element A.
#'
#' Extract the correlated information from factor X to txdb object made by element A.
#'
#' @param data Input factor X as bed format file, data.frame or grange.
#' @param txdb Element A as txdb format used for X annotation
#' @param tag tag before the names of column of output data
#' @param name ID of factor X used for annotation
#' @param strand if strand assign as TRUE means input data contain the strand information at column 6, if assign NULL strand will be *, otherwise could be "+" or "-"
#' @param verbose if TRUE, output the detail of processing
#' @param strand.INFO if TRUE, means the name of txdb contained strand information ("/+" and "/-")
#'
#'
#'  
ExtractCorrelate <- function(
    data, txdb=txdb, tag=NULL, name="name", strand=NULL, verbose=TRUE, strand.INFO=FALSE){

  #library(data.table); library(dplyr); library(ChIPseeker); library(stringr)

  if (is.null(tag)){ tag <- "Region" }

  ## Check the input bed format
  EXTRACT.NAME <- c("transcriptId", "annotation", "distanceToTSS")
  FINAL.NAME   <- c(name, tag, paste0(tag, "_annotation"), paste0(tag, "_dist"))

  if (class(data)[[1]]=="GRanges"){

    if (name%in%colnames(data.frame(data))){
      grange <- data
    } else {
      stop("Input data format is grange and please ckeck it contains, the name you assigned..")
    }

  } else {

    grange <- bedfromfile(data, name=name, strand=strand)

  }

  TABLE <- ChIPseeker::annotatePeak(
    grange, TxDb=txdb, verbose=verbose,
    genomicAnnotationPriority = c("Exon", "Intron", "Downstream", "Intergenic", "5UTR", "3UTR", "Promoter")) %>% 
    data.frame()

  if (!all(EXTRACT.NAME%in%colnames(TABLE))){
    stop("Check the ChIPseeker version..")
  }

  TABLE[TABLE$annotation%like%"Exon","distanceToTSS"] <- 0

  TABLE           <- TABLE[,c(name, EXTRACT.NAME)]
  colnames(TABLE) <- FINAL.NAME

  if (isTRUE(strand.INFO)){

    PLUS.IDX  <- stringr::str_split_fixed(TABLE[,2], "/", 2)[,2]%in%"+" %>% which()
    MINUS.IDX <- stringr::str_split_fixed(TABLE[,2], "/", 2)[,2]%in%"-" %>% which()


    if (length(PLUS.IDX)+length(MINUS.IDX)==nrow(TABLE)){

      TABLE[         ,2] <- stringr::str_split_fixed(TABLE[,2], "/", 2)[,1]
      TABLE[MINUS.IDX,4] <- TABLE[MINUS.IDX,4]*-1

      TABLE[TABLE[,4] ==0,3] <- "overlay"
      TABLE[TABLE[,4] < 0,3] <- "upstream"
      TABLE[TABLE[,4] > 0,3] <- "downstream"

    } else {

      stop("Check the information of column2 of region: function ExtractCorrelate")

    }


  } else {

    TABLE[TABLE[,4]==0,3] <- "overlay"
    TABLE[TABLE[,4]!=0,3] <- "no_overlay"

  }

  TABLE[,4] <- abs(TABLE[,4])

  return(TABLE)

}

#' Annotate the  distance of factor X to nearest element A.
#'
#' Extract the correlated information from factor X to txdb object made by element A.
#'
#' @param data Factor X as bed format file, data.frame or grange.
#' @param region Element A as bed format file, data.frame or grange. used for X annotation
#' @param tag tag before the names of column of output data
#' @param data.name ID of factor X used for annotation
#' @param region.name ID of element A used for annotation
#' @param verbose if TRUE, output the detail of processing
#'
#'  
RegionAnnoFromData <- function(
    data, region=region, tag=FALSE, data.name="name", region.name="name", verbose=TRUE){
  #library(data.table); library(dplyr); library(ChIPseeker)

  bed.plus   <- bedfromfile(region, name=region.name, strand="+") %>% data.frame()
  bed.minus  <- bedfromfile(region, name=region.name, strand="-") %>% data.frame()

  bed.plus$width  <- NULL
  bed.plus$name   <- paste(bed.plus$name,  "+",  sep="/")
  bed.minus$width <- NULL
  bed.minus$name  <- paste(bed.minus$name, "-",  sep="/")

  txdb <- GenomicRanges::makeGRangesFromDataFrame(
    rbind(bed.plus, bed.minus), keep.extra.columns=TRUE,
    starts.in.df.are.0based=FALSE, strand.field="strand") %>% txdbfromBed()

  TABLE <- bedfromfile(data, name=data.name, strand=NULL) %>%
    ExtractCorrelate(txdb=txdb, tag=tag, name=data.name, verbose=verbose, strand.INFO=TRUE)

  return(TABLE)

}



#' Radomly select peaks from a GRanges object.
#'
#' @param grange GRanges object.
#' @param seed Seed number.
#' @param n Number of peaks to select.
#' @param origin.n Number of peaks in grange.
#' 
#'  
randomPeak <- function(grange, seed=1, n=NULL, origin.n=NULL) {

  if (is.null(n)) {
    stop("Please assign a number to n")
  }

  if (is.null(origin.n)) {
    origin.n <- length(grange)
  }

  set.seed(seed)
  return(grange[sample(origin.n, n, replace=FALSE)])

}

#' Calculate specific peak's enrichment factor by compare with random peaks.
#' 
#' @param subject.name A character vector of peak name.
#' @param query A GRanges object of query peaks.
#' @param factor A GRanges object of factor peaks.
#' @param random.num Number of random times to get expect result.
#' @param pval P-value cutoff.
#' @param parallel If TRUE, will use parallel to calculate.
#' @param peak.name The column name of peak name in query.
#' 
#'  
PeakFactorStat <- function(
  subject.name, query=query, factor=factor, random.num=10000, pval=0.05, parallel=FALSE, peak.name="name") {

  if (!is.character(subject.name)) {
    stop("Please assign a character to subject.name")
  }

  if (!class(query)=="GRanges") {
    stop("Please assign a GRanges object in query")
  }

  if (!class(factor)=="GRanges") {
    stop("Please assign a GRanges object in factor")
  }

  if (!length(subject.name)==length(unique(subject.name))) {
    stop("Please assign a unique subject.name")
  }

  idx <- which(data.frame(query)[,peak.name]%in%subject.name)

  if (!length(idx)==length(subject.name)) {
    stop("Please assign a valid subject.name")
  } else {
    subject.n <- length(subject.name)
  }

  subject <- query[idx]

  observe <- GenomicRanges::findOverlaps(subject, factor) %>% 
    { data.frame(.)[,1] } %>% unique() %>% length()

  if (isTRUE(parallel)) {
    expect <- BiocParallel::bplapply(1:random.num, function(x) 
      randomPeak(query, seed=x, n=subject.n, origin.n=length(query)) %>% 
        { GenomicRanges::findOverlaps(., factor) } %>% 
        { data.frame(.)[,1] } %>% unique() %>% length()
    ) %>% unlist()
  } else {
    expect <- lapply(1:random.num, function(x) 
      randomPeak(query, seed=x, n=subject.n, origin.n=length(query)) %>% 
        { GenomicRanges::findOverlaps(., factor) } %>% 
        { data.frame(.)[,1] } %>% unique() %>% length()
    ) %>% unlist()
  }
  
  expect.mean <- mean(expect)
  expect.sd   <- sd(expect)

  log2FC  <- log2(observe/expect.mean)
  upper.p <- pnorm(observe, mean=expect.mean, sd=expect.sd, lower.tail=FALSE, log.p=TRUE)
  lower.p <- pnorm(observe, mean=expect.mean, sd=expect.sd, lower.tail=TRUE , log.p=TRUE)

  if (exp(lower.p) < pval) {
    type <- "lower"
  } else if (exp(upper.p) < pval) {
    type <- "upper"
  } else {
    type <- "none"
  }

  result <- list(
    type    = type,
    observe = observe,
    expect  = expect,
    log2FC  = log2FC,
    upper.p = upper.p,
    lower.p = lower.p)

  return(result)

}


# foreach(i=1:random.num, .combine="c") %dopar% {
#       randomPeak(subject, seed=i, n=observe) %>% 
#         GenomicRanges::findOverlaps(., factor) %>% 
#         { data.frame(.)[,1] } %>% unique() %>% length()


#' Calculate specific peak's (intersect with subject) enrichment factor by compare with random peaks.
#' 
#' @param subject A GRanges object of subject peaks.
#' @param query A GRanges object of query peaks.
#' @param factor A GRanges object of factor peaks.
#' @param random.num Number of random times to get expect result.
#' @param pval P-value cutoff.
#' @param parallel If TRUE, will use parallel to calculate.
#' @param peak.name The column name of peak name in query.
#' 
#'  
PeaktoPeakIntersect <- function(
  subject, query=query, factor=factor, random.num=10000, pval=0.05, parallel=FALSE, peak.name="name") {

  if (is.character(subject)) {
    if (file.exists(subject)) {
      message("Read subject from file")
      subject <- bedfromfile(subject)
    } else {
      stop("Please assign a valid file path")
    }
  }

  if (is.character(query)) {
    if (file.exists(query)) {
      message("Read query from file")
      query <- bedfromfile(query)
    } else {
      stop("Please assign a valid file path")
    }
  }

  if (is.character(factor)) {
    if (file.exists(factor)) {
      message("Read factor from file")
      factor <- bedfromfile(factor)
    } else {
      stop("Please assign a valid file path")
    }
  }

  if (!class(subject)=="GRanges") {
    stop("Please assign a GRanges object in subject")
  }

  if (!class(query)=="GRanges") {
    stop("Please assign a GRanges object in query")
  }

  if (!class(factor)=="GRanges") {
    stop("Please assign a GRanges object in factor")
  }

  ## Intersect the query to subject and get the number of intersected peaks
  subject.name <- GenomicRanges::findOverlaps(query, subject) %>% 
    { data.frame(.)[,1] } %>%
    unique() %>% { data.frame(query)[,peak.name][.] }

  result <- PeakFactorStat(
    subject.name=subject.name, query=query, factor=factor, 
    random.num=random.num, pval=pval, parallel=parallel, peak.name=peak.name)
  
  return(result)
}

#' Input shuffle peaks correlated (annotation) information as a list object.
#'
#' @param dir Directory of shuffle peak sets
#' @param parallel If set TRUE, will compile the shuffle peak sets with multiple threads
#' @param shuffle.n Number of peaks input into list. Reflect the flag (1~n) of shuffle_$flag.txt.gz
#' @param shuffle.prefix Name of prefix of peak sets.
#' @param file.ext Extension name of peak sets.
#' @param intersect If assign TRUE, output result will contained the number intersect peaks ($tag_dist==0) at first row
#' @param condition list of two number saperate with "-", will be used to count the number of factor A within the distance of numbers.
#'
#'
#'  
compileCorrelation <- function(
    dir, parallel=FALSE, shuffle.n=10000, shuffle.prefix="shuffle_", file.ext=".txt.gz",
    intersect=TRUE, condition=c("0-3000", "3000-10000", "10000-20000", "20000-30000", "40000-50000")){

  #library(dplyr)

  file.list <- paste0(shuffle.prefix, rep(1:shuffle.n), file.ext)

  if (!all(file.list%in%list.files(dir))){
    stop("Check the intput files in expected (shuffle peaksets) directory")
  }

  if (isTRUE(parallel)){

    expect.list <- BiocParallel::bplapply(
      file.path(dir, file.list), CountCorrelation, intersect=intersect,
      condition=condition)


  } else {

    expect.list <- NULL
    for (i in 1:length(file.list)){
      expect.list <- CountCorrelation(
        file.path(dir, file.list[i]), condition=condition)
    }

  }

  return(expect.list)

}


#' Input shuffle peaks correlated (annotation) information as a list object.
#'
#' @param dir Directory of shuffle peak sets
#' @param parallel If set TRUE, will compile the shuffle peak sets with multiple threads
#' @param shuffle.n Number of peaks input into list. Reflect the flag (1~n) of shuffle_$flag.txt.gz
#' @param shuffle.prefix Name of prefix of peak sets.
#' @param file.ext Extension name of peak sets.
#' @param intersect If assign TRUE, output result will contained the number intersect peaks ($tag_dist==0) at first row
#' @param condition list of two number saperate with "-", will be used to count the number of factor A within the distance of numbers.
#'
#'  
compileCorrelationByBin <- function(
    dir, parallel=FALSE, shuffle.n=10000, shuffle.prefix="shuffle_", file.ext=".txt.gz",
    intersect=TRUE, bin=1000, min=0, max=1000000, count.type="continue"){

  condition <- BinsDefine(bin=bin, min=min, max=max, type=count.type)
  result    <- compileCorrelation(
    dir, parallel=parallel, shuffle.n=shuffle.n, shuffle.prefix=shuffle.prefix,
    file.ext=file.ext, intersect=intersect, condition=condition)
  return(result)
}
#' Input shuffle peaks correlated (annotation) information as a list object.
#'
#' @param dir Directory of shuffle peak sets
#' @param parallel If set TRUE, will compile the shuffle peak sets with multiple threads
#' @param shuffle.n Number of peaks input into list. Reflect the flag (1~n) of shuffle_$flag.txt.gz
#' @param shuffle.prefix Name of prefix of peak sets.
#' @param file.ext Extension name of peak sets.
#' @param intersect If assign TRUE, output result will contained the number intersect peaks ($tag_dist==0) at first row
#' @param condition list of two number saperate with "-", will be used to count the number of factor A within the distance of numbers.
#'
#'
#'  
compileCorrelation <- function(
    dir, parallel=FALSE, shuffle.n=10000, shuffle.prefix="shuffle_", file.ext=".txt.gz",
    intersect=TRUE, condition=c("0-3000", "3000-10000", "10000-20000", "20000-30000", "40000-50000")){

  #library(dplyr)

  file.list <- paste0(shuffle.prefix, rep(1:shuffle.n), file.ext)

  if (!all(file.list%in%list.files(dir))){
    stop("Check the intput files in expected (shuffle peaksets) directory")
  }

  if (isTRUE(parallel)){

    expect.list <- BiocParallel::bplapply(
      file.path(dir, file.list), CountCorrelation, intersect=intersect,
      condition=condition)


  } else {

    expect.list <- NULL
    for (i in 1:length(file.list)){
      expect.list <- CountCorrelation(
        file.path(dir, file.list[i]), condition=condition)
    }

  }

  return(expect.list)

}


#' Input shuffle peaks correlated (annotation) information as a list object.
#'
#' @param dir Directory of shuffle peak sets
#' @param parallel If set TRUE, will compile the shuffle peak sets with multiple threads
#' @param shuffle.n Number of peaks input into list. Reflect the flag (1~n) of shuffle_$flag.txt.gz
#' @param shuffle.prefix Name of prefix of peak sets.
#' @param file.ext Extension name of peak sets.
#' @param intersect If assign TRUE, output result will contained the number intersect peaks ($tag_dist==0) at first row
#' @param condition list of two number saperate with "-", will be used to count the number of factor A within the distance of numbers.
#'
#'  
compileCorrelationByBin <- function(
    dir, parallel=FALSE, shuffle.n=10000, shuffle.prefix="shuffle_", file.ext=".txt.gz",
    intersect=TRUE, bin=1000, min=0, max=1000000, count.type="continue"){

  condition <- BinsDefine(bin=bin, min=min, max=max, type=count.type)
  result    <- compileCorrelation(
    dir, parallel=parallel, shuffle.n=shuffle.n, shuffle.prefix=shuffle.prefix,
    file.ext=file.ext, intersect=intersect, condition=condition)
  return(result)
}
