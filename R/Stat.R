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
    dir=NULL, observe=NULL, expect.dir=NULL, parallel=FALSE, shuffle.n=10000, shuffle.prefix="shuffle", observe.prefix="observe", file.ext=".txt.gz",
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
    dir=NULL, observe=NULL, expect.dir=NULL, parallel=FALSE, shuffle.n=10000, shuffle.prefix="shuffle", observe.prefix="observe", file.ext=".txt.gz",
    intersect=TRUE, bin=1000, min=0, max=1000000, count.type="continue",
    stat.type="both", tail="both", plotINFO.save=FALSE, log.p=TRUE){

  condition <- BinsDefine(bin=bin, min=min, max=max, type=count.type)
  result    <- ObsExpCompare(
    dir=dir, observe=observe, expect.dir=expect.dir, parallel=parallel, shuffle.n=shuffle.n, shuffle.prefix=shuffle.prefix,
    observe.prefix=observe.prefix, file.ext=file.ext, intersect=intersect, condition=condition,
    stat.type=stat.type, tail=tail, plotINFO.save=plotINFO.save, log.p=log.p)

  return(result)
}

