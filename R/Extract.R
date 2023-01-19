


ExtractStartEnd <- function(grange, strand=NULL){

  if (length(grange)>1){
    stop("Query peak must only contained 1 peak: ExtractStartEnd")
  }

  if (is.null(strand)){
    strand <- "+"
  }

  if (isTRUE(strand)){
    strand <- GenomicRanges::strand(grange)@values
  }

  if (strand=="+"){

    result <- data.frame(
      chr   = c(GenomicRanges::seqnames(grange), GenomicRanges::seqnames(grange)),
      start = c(GenomicRanges::start(grange),    GenomicRanges::end(grange)),
      end   = c(GenomicRanges::start(grange),    GenomicRanges::end(grange)))

  } else if (strand=="-"){

    result <- data.frame(
      chr   = c(GenomicRanges::seqnames(grange), GenomicRanges::seqnames(grange)),
      start = c(GenomicRanges::end(grange),      GenomicRanges::start(grange)),
      end   = c(GenomicRanges::end(grange),    GenomicRanges::start(grange)))

  } else {
    stop("Strand onlu can specify NULL, TRUE, +, -, strand is ", strand)
  }

  return(GenomicRanges::makeGRangesFromDataFrame(result))

}

DistCalculate <- function(query, subject=subject, name="name", DistIn=1000000){

  if (nrow(data.frame(query))>1){
    stop("Query peak must only contained 1 peak")
  }

  overlap.idx    <- data.frame(GenomicRanges::findOverlaps(query, subject))[,2]
  overlap.result <- rep(0, length(overlap.idx))

  if (length(overlap.idx) > 0){
    subject <- subject[-overlap.idx]
  }

  if (GenomicRanges::width(query)==1){
    query <- ExtractStartEnd(query+1, strand=TRUE)
  } else {
    query <- ExtractStartEnd(query, strand=TRUE)
  }

  table      <- data.frame(GenomicRanges::distanceToNearest(subject, query))
  upstream   <- table[table[,"subjectHits"]==1,]
  downstream <- table[table[,"subjectHits"]==2,]

  up.res <- upstream[,"distance"]
  names(up.res) <- data.frame(subject)[upstream[,"queryHits"],"name"]

  down.res <- downstream[,"distance"]
  names(down.res) <- data.frame(subject)[downstream[,"queryHits"],"name"]

  result <- c(overlap.result, -up.res, down.res)

  if (is.null(DistIn)){
    result <- result
  } else {
    result <- result[abs(result) < DistIn]
  }

  return(result)

}

CompilePeak <- function(
  query, subject=subject, name="name", DistIn=1000000, parallel=FALSE, save=NULL, strand=FALSE){

  if (class(query)=="GRanges"){

    next

  } else if (is.data.frame(query)){

    query <- bedfromfile(query)

  } else if (is.character(query)){

    if (file.exists(query)){
      query <- bedfromfile(query)
    } else {
      stop("File doesn't exist")
    }

  } else {
    stop("Check input query format")
  }

  if (isTRUE(parallel)){

    ## Transform grange to list
    query.list <- split(query, 1:length(query))
    result     <- BiocParallel::bplapply(query.list, DistCalculate, )

  } else {
    
  }

}
