#' Transform the list of per shuffle with condition to per condition with shuffle
#'
#' @param list Input the characters of conditions (type column in correlation table)
#' @param data Correlation table list of expect shuffle peak sets.
#' @param type The column name of conditions of range of distance.
#' @param number The column name of number of enriched regions.
#' @param names If TRUE, will name the output list with type/condition.
#'
#' @export
transformList <- function(list, data=data, type="type", number="number", names=TRUE){

  result <- NULL
  for (i in 1:length(list)){
    LIST <- NULL
    idx <-  which(data[[1]]$type%in%list[i])
    for (j in 1:length(data)){

      if (!all(colnames(data[[j]])%in%c(type, number))){
        stop("Check the input observe info")
      } else if (data[[j]][idx,type]==list[i]){
        LIST <- c(LIST, data[[j]][idx,number])
      } else {
        stop("Check data format in expect list...")
      }

    }; result[[i]] <- LIST
  }

  if (isTRUE(names)){
    names(result) <- list
  }

  return(result)

}

#' Create normal distribution table with input numbers
#'
#' @param data Input numbers
#'
#' @export
NormDistribute <- function(data){

  if (!is.numeric(data)){
    stop("Check the input data: func NormDistribute")
  }

  MEAN  <- mean(data)
  SD    <- sd(data)
  x     <- seq(-4, 4, length = 1000) * SD + MEAN
  y     <- dnorm(x, MEAN, SD)
  table <- data.frame(x=x, y=y)

  return(table)
}

#' Apply the statistic model to observe and expect peaks
#'
#' Use expected peak set result to build the empirical model and compare observe number to model to calculate the significance.
#'
#' @param observe Information of observed peaks
#' @param expect Information of expect peak sets list.
#' @param stat.type Statistic method. \cr
#' \cr
#' norm: Apply p.norm to calculate the significance \cr
#' \cr
#' ecdf : Apply ecdf function to calculate probability of upeer or lower. \cr
#' \cr
#' both : Both apply norm and ecdf.
#' @param tail Could assign as upper, lower or both.
#' @param plotINFO.save If assign as TRUE, will save the essential information to do some plot.
#' @param log.p If TRUE, will log2 the p.value (Only could apply on p.norm).
#' @param parallel If TRUE, will apply bapply to run the process.
#'
#' @export
compileStat <- function(
    observe=observe, expect=expect, stat.type="both", tail="both", plotINFO.save=FALSE, log.p=TRUE, parallel=FALSE){

  ## check the input data
  if (!all(colnames(observe)%in%c("type", "number"))){
    stop("Check the input observe info")
  } else {
    condition <- observe[,"type"]
  }

  if (isTRUE(parallel)){

    PRE.LIST       <- BiocParallel::bplapply(condition, transformList, data=expect, names=FALSE)
    EXPECT.LIST    <- NULL
    for (i in 1:length(PRE.LIST)){
      EXPECT.LIST[[i]] <- PRE.LIST[[i]][[1]]
    }

  } else {
    EXPECT.LIST <- transformList(condition, data=expect, names=FALSE)
  }

  if (length(condition)==length(EXPECT.LIST)){
    CONDITION.NUM <- 1:length(condition)
    OBSERVE.NUM   <- observe[,"number"]
  } else {
    stop("number of condition error, stop...")
  }

  result <- data.frame(type=condition)

  log2FC <- NULL
  for (i in CONDITION.NUM){
    log2FC <- c(log2FC, log2(OBSERVE.NUM[i]/mean(EXPECT.LIST[[i]])))
  }

  result$log2FC <- log2FC

  if (stat.type%in%c("norm", "both")){

    EXPECT.MEAN <- NULL
    EXPECT.SD   <- NULL
    for (i in CONDITION.NUM){
      EXPECT.MEAN <- c(EXPECT.MEAN, mean(EXPECT.LIST[[i]]))
      EXPECT.SD   <- c(EXPECT.SD,   sd(EXPECT.LIST[[i]]))
    }

    if (tail%in%c("both", "lower")){
      norm.p <- NULL
      for (i in CONDITION.NUM){
        norm.p <- c(
          norm.p, pnorm(OBSERVE.NUM[i], EXPECT.MEAN[i], EXPECT.SD[i], lower.tail=TRUE,  log.p=log.p))
      }; result$Norm_lower_P <- norm.p
    }

    if (tail%in%c("both", "upper")){
      norm.p <- NULL
      for (i in CONDITION.NUM){
        norm.p <- c(
          norm.p, pnorm(OBSERVE.NUM[i], EXPECT.MEAN[i], EXPECT.SD[i], lower.tail=FALSE, log.p=log.p))
      }; result$Norm_upper_P <- norm.p
    }

  }

  if (stat.type%in%c("ecdf", "both")){

    ECDF.LIST <- NULL
    for (i in CONDITION.NUM){
      ECDF.LIST[[i]] <- ecdf(EXPECT.LIST[[i]])
    }

    if (tail%in%c("both", "lower")){
      ecdf.p <- NULL
      for (i in CONDITION.NUM){
        ecdf.p <- c(ecdf.p, ECDF.LIST[[i]](OBSERVE.NUM[i]))
      }; result$ECDF_lower_P <- ecdf.p
    }

    if (tail%in%c("both", "upper")){
      ecdf.p <- NULL
      for (i in CONDITION.NUM){
        ecdf.p <- c(ecdf.p, 1-ECDF.LIST[[i]](OBSERVE.NUM[i]))
      }; result$ECDF_upper_P <- ecdf.p
    }

  }

  result$observe.num <- OBSERVE.NUM
  result$expect.mean <- EXPECT.MEAN

  if (isTRUE(plotINFO.save)){
    result.list <- list(
      statistic = result,
      expectNum = EXPECT.LIST)
    return(result.list)
  } else {
    return(result)
  }

}

#' Define the observe and expect annotation info to do the statistic.
#'
#' @param dir Directory contained observe and expect info files.
#' @param observe PATH to observe info.
#' @param expect.dir Directory contained expect info files.
#' @param parallel If TRUE, will apply bapply to run the process.
#' @param shuffle.n Times of shuffle.
#' @param shuffle.prefix Prefix name of shuffle info.
#' @param observe.prefix Prefix name of observe info.
#' @param file.ext Extension name of all files.
#' @param intersect If TRUE, result will contained intersect number.
#' @param condition Range of distance to nearest factor. Two number saperate by "-".
#' @param stat.type  Statistic method. \cr
#' \cr
#' norm: Apply p.norm to calculate the significance \cr
#' \cr
#' ecdf : Apply ecdf function to calculate probability of upeer or lower. \cr
#' \cr
#' both : Both apply norm and ecdf.
#' @param tail Could assign as upper, lower or both.
#' @param plotINFO.save If assign as TRUE, will save the essential information to do some plot.
#' @param log.p If TRUE, will log2 the p.value (Only could apply on p.norm).
#'
#' @export
ObsExpCompare <- function(
    dir=NULL, observe=NULL, expect.dir=NULL, parallel=FALSE, shuffle.n=10000, shuffle.prefix="shuffle_", observe.prefix="observe", file.ext=".txt.gz",
    intersect=TRUE, condition=c("0-3000", "3000-10000", "10000-20000", "20000-30000", "40000-50000"),
    stat.type="both", tail="both", plotINFO.save=FALSE, log.p=TRUE){

  #library(dplyr)


  if (!is.null(dir)){

    observe    <- file.path(dir, paste0(observe.prefix, file.ext))
    expect.dir <- dir
    message(
          "dir is not null, suppose observe=", observe,
          ", shuffle n =", shuffle.n, ", under ", dir)
      } else {

    if (any(is.null(observe), is.null(expect.dir))){
      stop("Please assign observe file and expect.dir or specify dir output by EnrichShuf.sh")
    }

  }

  observe.res <- CountCorrelation(observe, intersect=intersect, condition=condition)
  expect.res  <- compileCorrelation(
     expect.dir, parallel=parallel, shuffle.n=shuffle.n, shuffle.prefix=shuffle.prefix,
     file.ext=file.ext, intersect=intersect, condition=condition)

   result <- compileStat(
     observe=observe.res, expect=expect.res, stat.type=stat.type, tail=tail,
     plotINFO.save=plotINFO.save, log.p=log.p, parallel=parallel)

   return(result)

}

#' Define the observe and expect annotation info to do the statistic by each bin.
#'
#' @param dir Directory contained observe and expect info files.
#' @param observe PATH to observe info.
#' @param expect.dir Directory contained expect info files.
#' @param parallel If TRUE, will apply bapply to run the process.
#' @param shuffle.n Times of shuffle.
#' @param shuffle.prefix Prefix name of shuffle info.
#' @param observe.prefix Prefix name of observe info.
#' @param file.ext Extension name of all files.
#' @param intersect If TRUE, result will contained intersect number.
#' @param bin Interval size of each window.
#' @param min Minimum number of conditions
#' @param max Maximum number of conditions
#' @param type Could be specify: \cr
#' \cr
#' "continue" - All conditions will continues from min+interval*(max/bin(n)-1) to min+interval(max/bin(n))\cr
#' \cr
#' "within"   - All conditions will start from 0 to each intervals {0+interval(max/bin(n))}
#' @param stat.type  Statistic method. \cr
#' \cr
#' norm: Apply p.norm to calculate the significance \cr
#' \cr
#' ecdf : Apply ecdf function to calculate probability of upeer or lower. \cr
#' \cr
#' both : Both apply norm and ecdf.
#' @param tail Could assign as upper, lower or both.
#' @param plotINFO.save If assign as TRUE, will save the essential information to do some plot.
#' @param log.p If TRUE, will log2 the p.value (Only could apply on p.norm).
#'
#' @export
ObsExpCompareByBin <- function(
    dir=NULL, observe=NULL, expect.dir=NULL, parallel=FALSE, shuffle.n=10000, shuffle.prefix="shuffle_", observe.prefix="observe", file.ext=".txt.gz",
    intersect=TRUE, bin=1000, min=0, max=1000000, count.type="continue",
    stat.type="both", tail="both", plotINFO.save=FALSE, log.p=TRUE){

  condition <- BinsDefine(bin=bin, min=min, max=max, type=count.type)
  result    <- ObsExpCompare(
    dir=dir, observe=observe, expect.dir=expect.dir, parallel=parallel, shuffle.n=shuffle.n, shuffle.prefix=shuffle.prefix,
    observe.prefix=observe.prefix, file.ext=file.ext, intersect=intersect, condition=condition,
    stat.type=stat.type, tail=tail, plotINFO.save=plotINFO.save, log.p=log.p)

  return(result)
}

#' Compile the statistic result with intersect numbers from ObsExpCompare.
#' @param data Result from ObsExpCompare.
#' 
#' @export
ObsExpIntersectMerge <- function(data) {
  
  if (!all(names(data)%in%c("statistic", "expectNum"))) {
    stop("Please assign the result from ObsExpCompare")
  }
  
  type.list <- data$statistic[,"type"]
  
  if (!any(type.list%in%"intersect")) {
    stop("Didn't find intersect number in statistic result")
  }
  
  Num <- length(type.list)
  ## Expect number process
  expectNum <- lapply(1:Num, function(x) {
    if (x==1) {
      data$expectNum[[x]]
    } else {
      data$expectNum[[x]] + data$expectNum[[1]]
    }
  })
  
  ## Extract observe numbers
  observeNum <- lapply(1:Num, function(x) {
    if (x==1) {
      data$statistic[,"observe.num"][x]
    } else {
      data$statistic[,"observe.num"][x] + data$statistic[,"observe.num"][1]
    }
  }) %>% unlist()

  EXPECT.MEAN <- lapply(expectNum, function(x) mean(x)) %>% unlist()
  EXPECT.SD   <- lapply(expectNum, function(x) sd(x))   %>% unlist()
  
  log2FC  <- lapply(1:length(observeNum), function(x) 
    log2(observeNum[x]/EXPECT.MEAN[[x]])) %>% unlist()
  upper.p <- lapply(1:length(observeNum), function(x) 
    pnorm(observeNum[x], mean=EXPECT.MEAN[x], sd=EXPECT.SD[x], lower.tail=FALSE, log.p=TRUE)) %>% unlist()
  lower.p <- lapply(1:length(observeNum), function(x) 
    pnorm(observeNum[x], mean=EXPECT.MEAN[x], sd=EXPECT.SD[x], lower.tail=TRUE, log.p=TRUE))  %>% unlist()
  
  statistic <- data.frame(
    type         = type.list,
    observe.num  = observeNum,
    expect.mean  = EXPECT.MEAN,
    expect.sd    = EXPECT.SD,
    log2FC       = log2FC,
    Norm_upper_P = upper.p,
    Norm_lower_P = lower.p)

  result <- list(
    statistic = statistic,
    expectNum = expectNum)

  return(result)
  
}


#' Randomly select factors from a list.
#' 
#' @param list A character vector.
#' @param seed Seed number.
#' @param n Number of factors to select.
#' 
#' @export
randomFactor <- function(list, seed=1, n=NULL) {

  if (!is.character(list)) {
    stop("Please assign a character vector to list")
  }

  if (is.null(n)) {
    stop("Please assign a number to n")
  }

  set.seed(seed)
  return(sample(list, n, replace=TRUE)) ## true means could be sample will be put back to list

}

#' Associate factors to elements in a total list and calculate the significance
#' 
#' @param factor A character vector of factors.
#' @param total A character vector of total elements.
#' @param element A character vector of elements to be associated.
#' @param random.num Number of random times to get expect result.
#' @param log.p If TRUE, will log2 the p.value.
#' @param parallel If TRUE, will use parallel to do the calculation.
#' 
#' @export
TargetFactorSTAT <- function(
  factor, total=NULL, element=NULL, random.num=10000, log.p=FALSE, parallel=FALSE) {

  if (is.null(total)) {
    stop("Please assign a character to total")
  }

  if (is.null(element)) {
    stop("Please assign a character to element")
  }

  if (!is.character(factor)) {
    stop("Please assign a character vector to factor")
  }

  if (!is.character(total)) {
    stop("Please assign a character vector to total")
  }

  if (!is.character(element)) {
    stop("Please assign a character vector to element")
  }

  if (!all(element%in%total)) {
    message("Some elements are not in total... will remove them")
    element <- element[element%in%total]
  }


  total.factor     <- total[total%in%factor]
  total.factor.num <- length(total.factor)

  ## Associate factor with element
  observe.num <- sum(total.factor%in%element)

  ## Randomly select elements from total
  if (isTRUE(parallel)) {

    expect.num <- BiocParallel::bplapply(1:random.num, function(x) 
      sum(randomFactor(total, seed=x, n=total.factor.num)%in%element)) %>% unlist()

  } else {

    expect.num <- lapply(1:random.num, function(x) 
      sum(randomFactor(total, seed=x, n=total.factor.num)%in%element)) %>% unlist()

  }

  ## Calculate the statistic
  result <- data.frame(
    observe      = observe.num,
    expect       = mean(expect.num),
    log2FC       = log2(observe.num/mean(expect.num)),
    upper_pval   = pnorm(
      observe.num, mean=mean(expect.num), sd=sd(expect.num), lower.tail=FALSE, log.p=log.p),
    lower_pval   = pnorm(
      observe.num, mean=mean(expect.num), sd=sd(expect.num), lower.tail=TRUE , log.p=log.p))

  return(result)
}




#' Radomly select peaks from a GRanges object.
#'
#' @param grange GRanges object.
#' @param seed Seed number.
#' @param n Number of peaks to select.
#' @param origin.n Number of peaks in grange.
#' 
#' @export
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
#' @export
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
#' @export
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

