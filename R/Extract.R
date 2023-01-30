


ExtractStartEnd <- function(grange){

  if (length(grange)>1){
    stop("Query peak must only contained 1 peak: ExtractStartEnd")
  }

  strand <- GenomicRanges::strand(grange)@values

  if (strand=="*"){
    strand <- "*"
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
      end   = c(GenomicRanges::end(grange),      GenomicRanges::start(grange)))

  } else {
    stop("Strand onlu can specify NULL, TRUE, +, -, strand is ", strand)
  }

  return(GenomicRanges::makeGRangesFromDataFrame(result))

}


DistCalculate <- function(query, subject=subject, name="name", DistIn=1000000, strand=FALSE){

  if (nrow(data.frame(query))>1){
    stop("Query peak must only contained 1 peak")
  }

  if (isTRUE(strand)){

    overlap.idx    <- data.frame(GenomicRanges::findOverlaps(query, subject))[,2]
    overlap.result <- rep(0, length(overlap.idx))

    if (length(overlap.idx) > 0){
      subject <- subject[-overlap.idx]
    }

    if (GenomicRanges::width(query)==1){
      query <- ExtractStartEnd(query+1, strand=strand)
    } else {
      query <- ExtractStartEnd(query, strand=strand)
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

  } else {

    result <- GenomicRanges::distance(subject, query)
    names(result) <- data.frame(subject)[,name]

  }

  return(result)

}

DistPassSub <- function(
  list, query=query, subject=subject, name="name", DistIn=1000000, strand=FALSE){

  if (!all(names(list)%in%c("idx", "list"))){
    stop("")
  }

}

CompilePeak <- function(
  query, subject=subject, name="name", DistIn=1000000, parallel=FALSE, save=NULL, strand=FALSE){

  if (class(query)=="GRanges"){

    next

  } else if (is.data.frame(query)){

    query <- bedfromfile(query, name=name, strand=strand)

  } else if (is.character(query)){

    if (file.exists(query)){
      query <- bedfromfile(query, name=name, strand=strand)
    } else {
      stop("Query file doesn't exist")
    }

  } else {
    stop("Check input query format")
  }

  if (class(subject)=="GRanges"){

    next

  } else if (is.data.frame(subject)){

    subject <- bedfromfile(subject, name=name)

  } else if (is.character(subject)){

    if (file.exists(query)){
      query <- bedfromfile(query, name=name, strand=strand)
    } else {
      stop("Subject file doesn't exist")
    }

  } else {
    stop("Check input subject format")
  }

  if (!any(colnames(data.frame(subject))))

  OVERLAP.RESULT <- data.frame(GenomicRanges::findOverlaps(query+DistIn, subject))

  SIGN <- NULL; RESULT <- NULL
  for (i in 1:nrow(OVERLAP.RESULT)){
    if (is.null(SIGN)){
      Q.flag    <- OVERLAP.RESULT[i,1]
      S.numbers <- OVERLAP.RESULT[i,2]
      SIGN      <- TRUE
    } else {

      if (Q.flag==OVERLAP.RESULT[i,1]){
        S.numbers <- c(S.numbers, OVERLAP.RESULT[i,2])
      } else {
        LIST <- list(idx=Q.flag, list=S.numbers)
        RESULT[[paste("flag", Q.flag, sep="-")]] <- LIST
        Q.flag <- OVERLAP.RESULT[i,1]
        S.numbers <- OVERLAP.RESULT[i,2]
      }

    }
  }


  if (isTRUE(parallel)){

    ## Transform grange to list
    query.list <- split(query, 1:length(query))
    result     <- BiocParallel::bplapply(query.list, DistCalculate, )

  } else {

  }

}
