#' Input the bed format data.frame and convert to grange
#'
#' @param data input the data.frame or the path to bed file.
#' @param name assign the fourth column name and keep it exsist in grange
#' @param strand if strand assign as TRUE means input data contain the strand information at column 6, if assign NULL strand will be *, otherwise could be "+" or "-"
#'
#' @export
bedfromfile <- function(data, name="name", strand=NULL){

  library(GenomicFeatures)

  column.names <- c("chr", "start", "end", name)

  if (is.data.frame(data)){
    table <- data
  } else if (file.exists(data)){
    table <- read.table(data, header=FALSE)
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

  grange <- makeGRangesFromDataFrame(
    table, keep.extra.columns=TRUE, starts.in.df.are.0based=TRUE)

  if (any(isTRUE(strand), is.null(strand))){

    return(grange)

  } else if (strand%in%c("+", "-")){

    strand(grange) <- strand
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
#' @export
txdbfromBed <- function(data, name="name", strand=NULL){
  library(dplyr)

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

  txdb <- makeTxDbFromGRanges(grange)
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
#' @param
ExtractCorrelate <- function(data, txdb=txdb, tag=NULL, name="name", strand=NULL, verbose=TRUE){

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

  TABLE <- annotatePeak(grange, tssRegion=c(0,0), TxDb=txdb, verbose=verbose) %>% data.frame()

  if (!all(EXTRACT.NAME%in%colnames(TABLE))){
    stop("Check the ChIPseeker version..")
  }

  TABLE[TABLE$annotation%like%"Exon","distanceToTSS"] <- 0

  TABLE           <- TABLE[,c(name, EXTRACT.NAME)]
  colnames(TABLE) <- FINAL.NAME
  return(TABLE)

}

RegionAnnoFromData <- function(
    data, region=region, tag=FALSE, data.name="name", region.name="name", verbose=TRUE){
  library(data.table); library(dplyr); library(ChIPseeker)

  ## Import the tadb
  TXDB.plus  <- txdbfromBed(region, name=region.name, strand="+")
  TXDB.minus <- txdbfromBed(region, name=region.name, strand="-")

  grange <- bedfromfile(data, name=data.name, strand=NULL)

  TABLE.plus  <- ExtractCorrelate(grange, txdb=TXDB.plus,  tag="plus",  name=data.name, verbose=verbose)
  TABLE.minus <- ExtractCorrelate(grange, txdb=TXDB.minus, tag="minus", name=data.name, verbose=verbose)
  TABLE.all   <- merge(TABLE.plus, TABLE.minus, by=data.name)

  ### Table process ###
  OVERLAP     <- subset(TABLE.all, TABLE.all[,"plus_dist"]==0 | TABLE.all[,"minus_dist"]==0) %>%
    { .[,1:2] } %>%
    mutate(region=.[,2], type="overlay",    dist=0)
  OVERLAP[,2] <- NULL

  PLUS        <- subset(TABLE.all, abs(TABLE.all[,"plus_dist"]) <  abs(TABLE.all[,"minus_dist"])) %>%
    { .[,1:4] } %>%
    mutate(region=.[,2], type="no_overlay", dist=.[,4])
  PLUS[,2:4]  <- NULL

  MINUS       <- subset(TABLE.all, abs(TABLE.all[,"plus_dist"]) >  abs(TABLE.all[,"minus_dist"])) %>%
    { .[,c(1,5:7)] } %>%
    mutate(region=.[,2], type="no_overlay", dist=-.[,4])
  MINUS[,2:4] <- NULL

  EQUAL       <- subset(TABLE.all, abs(TABLE.all[,"plus_dist"]) == abs(TABLE.all[,"minus_dist"]) & abs(TABLE.all[,"plus_dist"]) > 0 & abs(TABLE.all[,"minus_dist"]) > 0) %>%
    { .[,1:4] } %>%
    mutate(region=.[,2], type="no_overlay", dist=.[,4])
  EQUAL[,2:4] <- NULL

  RESULT      <- rbind(OVERLAP, PLUS, MINUS, EQUAL)
  ### Table process ###

  if (is.character(tag)){
    names <- paste(tag, colnames(RESULT)[2:4], sep = "_")
  } else {
    names <- colnames(RESULT)[2:4]
  }

  colnames(RESULT)[2:4] <- names

  return(RESULT)

}
