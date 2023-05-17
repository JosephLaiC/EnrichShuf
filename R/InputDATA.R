#' Input the bed format data.frame and convert to grange
#'
#' @param data input the data.frame or the path to bed file.
#' @param name assign the fourth column name and keep it exsist in grange
#' @param strand if strand assign as TRUE means input data contain the strand information at column 6, if assign NULL strand will be *, otherwise could be "+" or "-"
#'
#' @export
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
#' @export
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


#' Input factors and associate to elements. Output data.frame conatined the information of factors to nearest elements with the distance.
#' 
#' @param factor input the factor data.frame or the path to bed file.
#' @param element input the element data.frame or the path to bed file.
#' @param strand.INFO if assign as TRUE means input data contain the strand information at column 6, if assign NULL strand will be *, otherwise could be "+" or "-"
FactorElementCorrelate <- function(
  factor, element=element, strand.INFO=TRUE) {
  
  ## Check factor format and input
  if (is.character(factor)) {
    
    if (!file.exists(factor)) {
      stop("Check the factor file exsist in path")
    }
    
    factor <- valr::bed_read(factor, n_fields=4)
    colnames(factor)[4] <- "factor_name" 
    
  } else if (is.data.frame(factor)) {
    
    factor <- factor[,1:4]
    colnames(factor) <- c("chrom", "start", "end", "factor_name")
    
  } else if (class(factor)[[1]]=="GRanges") {
    
    factor <- data.frame(factor)
    
    if (!all(c("seqnames", "start", "end", "name")%in%colnames(factor))) {
      stop("Check the factor grange contains the name information")
    }
    
    factor <- factor[,c("seqnames", "start", "end", "name")]
    colnames(factor) <- c("chrom", "start", "end", "factor_name")
    
  } else {
    stop("Check the factor format")
  }
  
  ## Check element format and input
  if (is.character(element)) {
    
    if (!file.exists(element)) {
      stop("Check the element file exsist in path")
    }
    
    element <- valr::bed_read(element, n_fields=4)
    colnames(element)[4] <- "element_name" 
    
  } else if (is.data.frame(element)) {
    
    element <- element[,1:4]
    colnames(element) <- c("chrom", "start", "end", "element_name")
    
  } else if (class(element)[[1]]=="GRanges") {
    
    element <- data.frame(element)
    
    if (!all(c("seqnames", "start", "end", "name")%in%colnames(element))) {
      stop("Check the element grange contains the name information")
    }
    
    element <- element[,c("seqnames", "start", "end", "name")]
    colnames(element) <- c("chrom", "start", "end", "element_name")
    
  } else {
    stop("Check the element format")
  }
  
  table <- valr::bed_closest(factor, element)[
    ,c("factor_name.x", "element_name.y", ".dist", ".overlap")] %>% data.frame()
  colnames(table) <- c("factor_name", "element_name", "distance", "overlap")
  
  intersect <- table[table$distance==0,] %>% 
    { .[order(.$overlap, decreasing=TRUE),] }
  
  multi.overlap <- table(intersect$factor_name) %>% 
    { .[.>1] } %>% names()
  
  if (length(multi.overlap) > 0) {
    
    multi.info <- lapply(multi.overlap, function(x) 
      intersect[intersect$factor_name==x,] %>% { .[1,] }) %>%
      Reduce(rbind, .)
    
    intersect.info <- rbind(
      intersect[!intersect$factor_name%in%multi.overlap,], multi.info) 
    
  } else {
    
    intersect.info <- intersect
    
  }
  
  associate <- table[!table$factor_name%in%intersect.info$factor_name,]
  multi.associate <- table(associate$factor_name) %>% 
    { .[.>1] } %>% names()
  
  if (length(multi.associate) > 0) {
    
    multi.info <- lapply(multi.associate, function(x) 
      associate[associate$factor_name==x,] %>% { .[1,] }) %>%
      Reduce(rbind, .)
    
    associate.info <- rbind(
      associate[!associate$factor_name%in%multi.associate,], multi.info) 
    
  } else {
    
    associate.info <- associate
    
  }
  
  result <- rbind(intersect.info, associate.info)
  return(result)
  
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
#' @export
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
#' @export
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
