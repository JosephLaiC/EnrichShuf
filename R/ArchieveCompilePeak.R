#' Extract start and end position from grange
#'
#' @param grange input grange
#'
#' @export
ExtractStartEnd <- function(grange, strand=NULL){

  if (length(grange)>1){
    stop("Query peak must only contained 1 peak: ExtractStartEnd")
  }

  strand <- GenomicRanges::strand(grange)@values

  if (is.null(strand)){
    strand <- "+"
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
    stop("Strand onlu can specify *, +, -, strand is ", strand)
  }

  return(GenomicRanges::makeGRangesFromDataFrame(result))

}

#' Extract distance between query and subject within a distance threshold
#'
#' @param query input query peak
#' @param subject input subject peak
#' @param name assign the fourth column name and keep it exsist in grange
#' @param DistIn assign the distance threshold
#' @param strand if strand assign as TRUE means input data contain the strand information at column 6, if assign NULL strand will be *, otherwise could be "+" or "-"
#'
#' @export
DistCalculate <- function(
    query, subject=subject, name="name", DistIn=1000000, OVERLAY.ONLY=FALSE, strand=NULL, verbose=FALSE){

  if (!nrow(data.frame(query))==1) {
    stop("Query peak must only contained 1 peak")
  } else if (!is.numeric(DistIn)) {
    stop("Distance threshold must be numeric")
  }

  ## Search the global overlap region
  idx <- data.frame(GenomicRanges::findOverlaps(query+DistIn, subject))[,2]

  if (length(idx) > 0){
    subject <- subject[idx]
  } else {
    return(NULL)
  }

  overlap.idx    <- data.frame(GenomicRanges::findOverlaps(query, subject))[,2]
  overlap.result <- rep(0, length(overlap.idx)) ## 0 is the overlap result, number is equal to the length of overlap index

  if (any(isTRUE(OVERLAY.ONLY), DistIn==0)){

    if (isTRUE(verbose)){
      message("OVERLAY.ONLY is TRUE or DistIn is 0, only return the overlap region")
    }

    if (length(overlap.idx)==0) {
      return(NULL)
    } else {
      return(overlap.result)
    }


  }

  ## remove the overlap region from subject
  if (length(overlap.idx) > 0) {
    subject <- subject[-overlap.idx]
  }


  if (isTRUE(strand)){

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

  } else {

    result <- GenomicRanges::distance(subject, query)
    names(result) <- data.frame(subject)[,name]

  }

  if (is.null(DistIn)){
    result <- result
  } else {
    result <- result[abs(result) < DistIn]
  }

  return(result)

}


# DistPassSub <- function(
#   list, query=query, subject=subject, name="name", DistIn=1000000, strand=FALSE){

#   if (!all(names(list)%in%c("idx", "list"))){
#     stop("")
#   }

# }

#' Compile the distance information of factors to elements within distance
#'
#' @param query input query peak
#' @param subject input subject peak
#' @param name assign the fourth column name and keep it exsist in grange
#' @param DistIn assign the distance threshold
#' @param parallel if assign TRUE, the function will run in parallel
#' @param save if assign TRUE, the function will save the result
#' @param strand if strand assign as TRUE means input data contain the strand information at column 6, if assign NULL strand will be *, otherwise could be "+" or "-"
#'
#' @export
CompilePeak <- function(
  query, subject=subject, name="name", DistIn=1000000, parallel=FALSE, save=NULL, strand=NULL){

  # Import query
  if (is.data.frame(query)) {

    query <- bedfromfile(query, name=name, strand=strand)

  } else if (is.character(query)) {

    if (file.exists(query)) {
      query <- bedfromfile(query, name=name, strand=strand)
    } else {
      stop("Query file doesn't exist")
    }

  }

  if (!class(query)=="GRanges"){
    stop("Check input query format")
  }

  # Import subject
  if (is.data.frame(subject)) {

    subject <- bedfromfile(subject, name=name)

  } else if (is.character(subject)) {

    if (file.exists(subject)) {
      subject <- bedfromfile(subject, name=name)
    } else {
      stop("Subject file doesn't exist")
    }

  }

  if (!class(subject)=="GRanges") {
    stop("Subject must be GRanges")
  }

  # SIGN <- NULL; RESULT <- NULL
  # for (i in 1:nrow(OVERLAP.RESULT)){
  #   if (is.null(SIGN)){
  #     Q.flag    <- OVERLAP.RESULT[i,1]
  #     S.numbers <- OVERLAP.RESULT[i,2]
  #     SIGN      <- TRUE
  #   } else {

  #     if (Q.flag==OVERLAP.RESULT[i,1]){
  #       S.numbers <- c(S.numbers, OVERLAP.RESULT[i,2])
  #     } else {
  #       LIST <- list(idx=Q.flag, list=S.numbers)
  #       RESULT[[paste("flag", Q.flag, sep="-")]] <- LIST
  #       Q.flag <- OVERLAP.RESULT[i,1]
  #       S.numbers <- OVERLAP.RESULT[i,2]
  #     }

  #   }
  # }

  if (isTRUE(parallel)){

    ## Transform grange to list
    query.list <- split(query, 1:length(query))
    library(BiocParallel)
    result     <- BiocParallel::bplapply(
      query.list, DistCalculate, subject=subject, name=name, DistIn=DistIn, strand=strand)

  } else {

    result <- lapply(1:length(query), function(i){
      DistCalculate(query[i], subject=subject, name=name, DistIn=DistIn, strand=strand)
    })

  }

  names(result) <- data.frame(query)[,name]
  return(result)

}

#' Compare the observe compile information with expect compile information by normal distribution
#' 
#' @param observe A numeric vector of observe data
#' @param expect.data A list of numeric vector of expect data
#' @param parallel The number of cores to use
#' @param parallel.type  Could be specify one of: \cr
#' \cr
#' "mclapply" - Use mclapply to run in parallel\cr
#' \cr
#' "bplapply" - Use BiocParallel to run in parallel
#'
#' @export
normalDistPeakCompile <- function(
    observe, expect.data=expect.data, 
    parallel=1, parallel.type="mclapply"){
  
  if (!is.numeric(observe)) {
    stop("Check the observe data")
  }
  
  lapply(expect.data, function(x){
    
    if (!is.numeric(x)) {
      stop("Check the expect data")
    }
    
    if (!identical(names(x), names(observe))) {
      stop("Check the names of expect and observe data")
    }
    
  })
  
  if (parallel > 1) {

    idx.num <- seq(1, length(observe), ceiling(length(observe)/parallel))
    idx <- data.frame(
      start = idx.num,
      end = c(idx.num[-1]-1, length(observe)))
    
    if (parallel.type=="mclapply") {
      
      gc(verbose = FALSE)
      STAT.INFO <- parallel::mclapply(1:nrow(idx), function(x){
        
        start <- idx[x,1]; end <- idx[x,2]
        lapply(start:end, function(y){
          
          lapply(expect.data, function(z){z[[y]]}) %>% unlist() %>%
            { list(mean=mean(.), sd=sd(.)) }
          
        })
        
      }, mc.cores = parallel) %>% Reduce(c, .)
      
    }
    
  } else {
    
    STAT.INFO <- lapply(1:length(observe), function(x){
      
      lapply(expect.data, function(y){y[[x]]}) %>% unlist() %>%
        { list(mean=mean(.), sd=sd(.)) }
      
    })
    
  }
  
  
  
  if (parallel > 1) {
    
    idx.num <- seq(1, length(observe), ceiling(length(observe)/parallel))
    idx <- data.frame(
      start = idx.num,
      end = c(idx.num[-1]-1, length(observe)))
    
    if (parallel.type=="mclapply") {
      
      gc(verbose = FALSE)
      result <- parallel::mclapply(1:nrow(idx), function(x){
        
        start <- idx[x,1]; end <- idx[x,2]
        lapply(start:end, function(y){
          data.frame(
            name    = names(observe)[y],
            observe = observe[y],
            expect  = STAT.INFO[[y]]$mean,
            log2FC  = log2(observe[y]/STAT.INFO[[y]]$mean),
            upper.p = pnorm(
              observe[y], mean=STAT.INFO[[y]]$mean, sd=STAT.INFO[[y]]$sd, lower.tail=FALSE),
            lower.p = pnorm(
              observe[y], mean=STAT.INFO[[y]]$mean, sd=STAT.INFO[[y]]$sd, lower.tail=TRUE))
        }) %>% Reduce(rbind, .)
        
      }, mc.cores = parallel) %>% Reduce(rbind, .)
      
    } else if (parallel.type=="bplapply") {
      
      gc(verbose = FALSE)
      result <- parallel::bplapply(1:nrow(idx), function(x){
        
        start <- idx[x,1]; end <- idx[x,2]
        lapply(start:end, function(y){
          data.frame(
            name    = names(observe)[y],
            observe = observe[y],
            expect  = expect.data[[y]]$mean,
            log2FC  = log2(observe[y]/expect.data[[y]]$mean),
            upper.p = pnorm(
              observe[y], mean=expect.data[[y]]$mean, sd=expect.data[[y]]$sd, lower.tail=FALSE),
            lower.p = pnorm(
              observe[y], mean=expect.data[[y]]$mean, sd=expect.data[[y]]$sd, lower.tail=TRUE))
        }) %>% Reduce(rbind, .)
        
      }) %>% Reduce(rbind, .)
      
    } else {
      
      stop("Check the parallel.type")
      
    }
    
    
  } else {
    
    result <- lapply(1:length(observe), function(x){
      data.frame(
        name    = names(observe)[x],
        observe = observe[x],
        expect  = expect.data[[x]]$mean,
        log2FC  = log2(observe[x]/expect.data[[x]]$mean),
        upper.p = pnorm(
          observe[x], mean=expect.data[[x]]$mean, sd=expect.data[[x]]$sd, lower.tail=FALSE),
        lower.p = pnorm(
          observe[x], mean=expect.data[[x]]$mean, sd=expect.data[[x]]$sd, lower.tail=TRUE))
    }) %>% Reduce(rbind, .)
    
  }
  
  return(result)
  
}
