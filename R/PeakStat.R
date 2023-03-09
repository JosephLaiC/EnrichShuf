

#' Combine the observe and expect info to do the statistic.
#' 
#' @param observe Observe info.
#' @param expect Expect info.
#' 
#' @export
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
#' @export
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
#' @export
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
#' 
#' @export
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
#' @export
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
#' @export
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
#' 
#' @export
ObsExpBinomial <- function(data, n=10) {

  if (!"BinomTable"%in%names(data)) {
    message("The BinomTable is not found, running ObsExpBinomTable...")
    data <- ObsExpBinomTable(data)
  }

  up.res   <- lapply(n, function(x) 
    lapply(rowSums(data[["BinomTable"]][["up"]][,1:x]), function(y)
      binom.test(y, x, alternative="greater")$p.value) %>% unlist())
  down.res <- lapply(n, function(x)
    lapply(rowSums(data[["BinomTable"]][["down"]][,1:x]), function(y)
      binom.test(y, x, alternative="less")$p.value) %>% unlist())
  two.res  <- lapply(n, function(x)
    lapply(rowSums(data[["BinomTable"]][["up"]][,1:x]), function(y)
      binom.test(y, x, alternative="two.sided")$p.value) %>% unlist())
  list.res <- lapply(n, function(x) colnames(data$BinomTable[[1]])[1:x])

  names(up.res)   <- paste0("Binom_", n)
  names(down.res) <- paste0("Binom_", n)
  names(two.res)  <- paste0("Binom_", n)
  names(list.res) <- paste0("Binom_", n)

  data$Binom_Pval <- list(
    up       = up.res,
    down     = down.res,
    two.tail = two.res,
    list     = list.res)
  return(data)
}

#' Compile the binomial test results.
#' 
#' @param data The data should be a list with observe, expect and log2FC, export by ObsExpCompile or ObsExpRDS and contsins the BinomTable.
#' @param p.adjust The p.adjust method, default is NULL.
#' 
#' @export
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
      pVal_down     = data$Binom_Pval$down[[i]],
      pVal_two_tail = data$Binom_Pval$two.tail[[i]])
  }

  if (!is.null(p.adjust)) {
    for (i in names(result.list)) {
      result.list[[i]]$FDR_up       <- p.adjust(result.list[[i]]$pVal_up, method=p.adjust)
      result.list[[i]]$FDR_down     <- p.adjust(result.list[[i]]$pVal_down, method=p.adjust)
      result.list[[i]]$FDR_two_tail <- p.adjust(result.list[[i]]$pVal_two_tail, method=p.adjust)
    }
  }

  data$Binom_compile <- result.list
  return(data)
}