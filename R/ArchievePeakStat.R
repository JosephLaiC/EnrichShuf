

#' Combine the observe and expect info to do the statistic.
#' 
#' @param observe Observe info.
#' @param expect Expect info.
#' 
#' @noRd
ObsExpCompile <- function(expect, observe=observe) {

  peak.set <- observe
  observe  <- lapply(observe, function(x) length(x)) %>% unlist()
  expect   <- lapply(expect , function(x) length(x)) %>% unlist()

  result <- list(
    PeakSet = peak.set,
    observe = observe,
    expect  = expect,
    log2FC  = log2(observe+1)-log2(expect+1))

}

#' Combine the observe and multiple expect info to do the statistic.
#' 
#' @param observe Observe info.
#' @param expect Expect info.
#' @param parallel If TRUE, will apply bapply to run the process.
#' @param expect.name list of name for expect results. If NULL, will assign as "Shuffle_1", "Shuffle_2"...
#' 
#' @noRd
ObsExpCompileList <- function(
  expect, observe=observe, parallel=FALSE, expect.name=NULL) {
  
  peak.set <- observe
  observe  <- lapply(observe, function(x) length(x)) %>% unlist()

  if (isTRUE(parallel)) {

    expect <- BiocParallel::bplapply(expect, function(x)
      lapply(x, function(y) length(y)) %>% unlist())
    log2FC <- BiocParallel::bplapply(expect, function(x) log2(observe+1)-log2(x+1))

  } else {
      
    expect <- lapply(expect, function(x)
      lapply(x, function(y) length(y)) %>% unlist())
    log2FC <- lapply(expect, function(x) log2(observe+1)-log2(x+1))

  }

  if (is.null(expect.name)) {
    names(expect) <- paste0("Shuffle_", 1:length(expect))
    names(log2FC) <- paste0("Shuffle_", 1:length(log2FC))
  } else {
    names(expect) <- expect.name
    names(log2FC) <- expect.name
  }

  result <- list(
    PeakSet = peak.set,
    expect  = expect,
    observe = observe,
    log2FC  = log2FC)

}

#' Extract the distance information from a vector.
#' 
#' @param data Input data.
#' @param dist Distance to calculate the number of peaks.
#' @param intersect If TRUE, will extract the distance smaller than dist. If FALSE, will extract the distance larger than 0 and smaller than dist.
#' 
#' @noRd
ExtractDistance <- function(data, dist=dist, intersect=TRUE) {

  if (isTRUE(intersect)) {
    data <- data[data<=dist]
  } else {
    data <- data[data>0 & data<=dist]
  }

  return(data)
}


#' Combine the observe and expect info from RDS files to do the statistic.
#' 
#' @param shuffleRDS PATH to shuffle info. Could assign as a list of RDS files.
#' @param observeRDS PATH to observe info.
#' @param dist Distance to calculate the number of peaks.
#' @param parallel If TRUE, will apply bapply to run the process.
#' @param expect.name list of name for expect results. If NULL, will assign as "Shuffle_1", "Shuffle_2"...
#' @param intersect If TRUE, will extract the distance smaller than dist. If FALSE, will extract the distance larger than 0 and smaller than dist.
#' 
#' @noRd
ObsExpRDS <- function(
  shuffleRDS, observeRDS=observeRDS, dist=1000000, parallel=FALSE, expect.name=NULL, intersect=TRUE) {

  library(dplyr)

  if (!is.character(observeRDS)) {
    stop("The observeRDS should be character.")
  }

  if (!is.character(shuffleRDS)) {
    stop("The shuffleRDS should be character.")
  }
    
  if (file.exists(observeRDS)) {
    observe <- readRDS(observeRDS) %>% 
      { lapply(., function(x) ExtractDistance(x, dist=dist, intersect=intersect)) }
  } else {
    stop("The observeRDS file is not exists.")
  }

  if (length(shuffleRDS)==1) {
      
    if (file.exists(shuffleRDS)) {

      expect <- readRDS(shuffleRDS) %>%
        { lapply(., function(x) ExtractDistance(x, dist=dist, intersect=intersect)) }

    } else {
      stop("The shuffleRDS file is not exists.")
    }

    result <- ObsExpCompile(expect, observe=observe)

  } else {

    lapply(shuffleRDS, function(x) {
          
      if (!file.exists(x)) {
        stop("The shuffleRDS file is not exists.")
      } 
          
    })

    if (isTRUE(parallel)) {

      expect <- BiocParallel::bplapply(shuffleRDS, function(x) 
        readRDS(x) %>% { lapply(., function(x) ExtractDistance(x, dist=dist, intersect=intersect)) } )
      

    } else {

      expect <- lapply(shuffleRDS, function(x) 
        readRDS(x) %>% { lapply(., function(x) ExtractDistance(x, dist=dist, intersect=intersect)) } )

    }

    result <- ObsExpCompileList(
      expect, observe=observe, parallel=parallel, expect.name=expect.name)

  }

  return(result)
}

#' Perform the statistic for the observe and expect info by exact testing based on edgeR.
#' 
#' @param data The data should be a list with observe, expect and log2FC, export by ObsExpCompile or ObsExpRDS.
#' @param type The type of statistic. Could be "binomial", "exact_2tail", "exact_sample_compare", "exact_deviance", "exact_binomial".
#' @param parallel If TRUE, will apply bapply to run the process.
#' 
#' @noRd
ObsExpCompileSTAT <- function(data, type=NULL, parallel=FALSE) {

  if (!all(c("observe", "expect", "log2FC")%in%names(data))) {
    stop("The data should be a list with observe, expect and log2FC, export by ObsExpCompile or ObsExpRDS.")
  }

  if (is.character(type)) {

    type.list <- c(
      "binomial", "exact_2tail", "exact_sample_compare", "exact_deviance", "exact_binomial")

    if (!all(type%in%type.list)) {
      stop("The type should be one of ", paste(type.list, collapse=", "))
    }

  }

  pval <- NULL
  if (is.list(data[["expect"]])) {
    
    if (isTRUE(parallel)) {

      if ("binomial"%in%type) {
        pval[["binomial"]] <- BiocParallel::bplapply(data[["expect"]], function(x) 
          edgeR::binomTest(data[["observe"]], x))
      } 
      
      if ("exact_2tail"%in%type) {
        pval[["exact_2tail"]] <- BiocParallel::bplapply(data[["expect"]], function(x) 
          edgeR::exactTestDoubleTail(data[["observe"]], x))
      } 
      
      if ("exact_sample_compare"%in%type) {
        pval[["exact_sample_compare"]] <- BiocParallel::bplapply(data[["expect"]], function(x) 
          edgeR::exactTestBySmallP(data[["observe"]], x))
      } 
      
      if ("exact_deviance"%in%type) {
        pval[["exact_deviance"]] <- BiocParallel::bplapply(data[["expect"]], function(x) 
          edgeR::exactTestByDeviance(data[["observe"]], x))
      } 
      
      if ("exact_binomial"%in%type) {
        pval[["exact_deviance"]] <- BiocParallel::bplapply(data[["expect"]], function(x) 
          edgeR::exactTestBetaApprox(data[["observe"]], x))
      }

    } else {
       
      if ("binomial"%in%type) {
        pval[["binomial"]] <- lapply(data[["expect"]], function(x) 
          edgeR::binomTest(data[["observe"]], x))
      } 
      
      if ("exact_2tail"%in%type) {
        pval[["exact_2tail"]] <- lapply(data[["expect"]], function(x) 
          edgeR::exactTestDoubleTail(data[["observe"]], x))
      } 
      
      if ("exact_sample_compare"%in%type) {
        pval[["exact_sample_compare"]] <- lapply(data[["expect"]], function(x) 
          edgeR::exactTestBySmallP(data[["observe"]], x))
      } 
      
      if ("exact_deviance"%in%type) {
        pval[["exact_deviance"]] <- lapply(data[["expect"]], function(x) 
          edgeR::exactTestByDeviance(data[["observe"]], x))
      } 
      
      if ("exact_binomial"%in%type) {
        pval[["exact_deviance"]] <- lapply(data[["expect"]], function(x) 
          edgeR::exactTestBetaApprox(data[["observe"]], x))
      }

    }

  } else {

    if ("binomial"%in%type) {
      pval[["binomial"]] <- edgeR::binomTest(data[["observe"]], data[["expect"]])
    } 
      
    if ("exact_2tail"%in%type) {
      pval[["exact_2tail"]] <- edgeR::exactTestDoubleTail(data[["observe"]], data[["expect"]])
    } 
      
    if ("exact_sample_compare"%in%type) {
      pval[["exact_sample_compare"]] <- edgeR::exactTestBySmallP(data[["observe"]], data[["expect"]])
    } 
      
    if ("exact_deviance"%in%type) {
      pval[["exact_deviance"]] <- edgeR::exactTestByDeviance(data[["observe"]], data[["expect"]])
    } 
      
    if ("exact_binomial"%in%type) {
      pval[["exact_deviance"]] <- edgeR::exactTestBetaApprox(data[["observe"]], data[["expect"]])
    }

  }

  data[["ExactSTAT"]] <- pval

  return(data)
}

#' Get up or down info from the info of observe and expect to perform binomial test.
#' 
#' @param data The data should be a list with observe, expect and log2FC, export by ObsExpCompile or ObsExpRDS.
#' 
#' @noRd
ObsExpBinomTable <- function(data) {

  if (!all(c("observe", "expect", "log2FC")%in%names(data))) {
    stop("The data should be a list with observe, expect and log2FC, export by ObsExpCompile or ObsExpRDS.")
  }

  if (!is.list(data[["expect"]])) {
    stop("The expect should be a list, to generate the table, please use ObsExpBinomTableList.")
  }

  up.res   <- lapply(data[["expect"]], function(x) data[["observe"]] > x) %>% data.frame()
  down.res <- lapply(data[["expect"]], function(x) data[["observe"]] < x) %>% data.frame()

  colnames(up.res)   <- names(data[["expect"]])
  colnames(down.res) <- names(data[["expect"]])
  
  data$BinomTable <- list(up=up.res, down=down.res)
  return(data)
}


#' Perform binomial test for the up and down targets.
#' 
#' @param data The data should be a list with observe, expect and log2FC, export by ObsExpCompile or ObsExpRDS and contsins the BinomTable.
#' @param n The number of bins to perform the binomial test.
#' @param parallel If TRUE, use BiocParallel to perform the binomial test.
#' 
#' @noRd
ObsExpBinomial <- function(data, n=10, parallel=FALSE) {

  if (!"BinomTable"%in%names(data)) {
    message("The BinomTable is not found, running ObsExpBinomTable...")
    data <- ObsExpBinomTable(data)
  }

  if (isTRUE(parallel)) {
    up.res   <- BiocParallel::bplapply(n, function(x) 
      lapply(rowSums(data[["BinomTable"]][["up"]][,1:x]), function(y)
      binom.test(y, x, alternative="greater")$p.value) %>% unlist())
    down.res <- BiocParallel::bplapply(n, function(x)
      lapply(rowSums(data[["BinomTable"]][["down"]][,1:x]), function(y)
      binom.test(y, x, alternative="greater")$p.value) %>% unlist())
  } else {
    up.res   <- lapply(n, function(x) 
      lapply(rowSums(data[["BinomTable"]][["up"]][,1:x]), function(y)
      binom.test(y, x, alternative="greater")$p.value) %>% unlist())
    down.res <- lapply(n, function(x)
      lapply(rowSums(data[["BinomTable"]][["down"]][,1:x]), function(y)
      binom.test(y, x, alternative="greater")$p.value) %>% unlist())
  }

  list.res <- lapply(n, function(x) colnames(data$BinomTable[[1]])[1:x])

  names(up.res)   <- paste0("Binom_", n)
  names(down.res) <- paste0("Binom_", n)
  names(list.res) <- paste0("Binom_", n)

  data$Binom_Pval <- list(
    up       = up.res,
    down     = down.res,
    list     = list.res)
  return(data)
}

#' Compile the binomial test results.
#' 
#' @param data The data should be a list with observe, expect and log2FC, export by ObsExpCompile or ObsExpRDS and contsins the BinomTable.
#' @param p.adjust The p.adjust method, default is NULL.
#' 
#' @noRd
ObsExpBinomCompile <- function(data, p.adjust=NULL) {

  if (!"Binom_Pval"%in%names(data)) {
    stop("The Binom_Pval is not found, please run ObsExpBinomial...")
  }

  observe      <- data$observe
  expect.table <- data.frame(data$expect)
  col.list     <- data$Binom_Pval$list

  result.list <- NULL
  for (i in names(col.list)) {
    list <- col.list[[i]]

    result.list[[i]] <- data.frame(
      log2FC        = log2(observe+1) - log2(rowMeans(expect.table[,list])+1),
      pVal_up       = data$Binom_Pval$up[[i]],
      pVal_down     = data$Binom_Pval$down[[i]])
  }

  if (!is.null(p.adjust)) {
    for (i in names(result.list)) {
      result.list[[i]]$FDR_up       <- p.adjust(result.list[[i]]$pVal_up, method=p.adjust)
      result.list[[i]]$FDR_down     <- p.adjust(result.list[[i]]$pVal_down, method=p.adjust)
    }
  }

  data$Binom_compile <- result.list
  return(data)
}



#' Transform the list of per shuffle with condition to per condition with shuffle
#'
#' @param list Input the characters of conditions (type column in correlation table)
#' @param data Correlation table list of expect shuffle peak sets.
#' @param type The column name of conditions of range of distance.
#' @param number The column name of number of enriched regions.
#' @param names If TRUE, will name the output list with type/condition.
#'
#' @noRd
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
#' @noRd
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
#' @noRd
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
#' @noRd
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
#' @noRd
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
#' @noRd
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
