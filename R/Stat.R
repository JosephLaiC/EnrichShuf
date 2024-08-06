
#' Calculate the p-value by comparing the observed value to an empirical model constructed from 
#' expected values with specific range of distances (condition).
#' 
#' @param data The input object contained the information of observed and expected number.
#' @param name A character used to extract a specific object for calculating significance by comparing observed 
#' values with empirical models built from expected values.
#' @param log.p If set to TRUE, the logarithm of the p-value will be returned.
#' 
#' @export
ObsExpSTATbyName <- function(data, name=name, log.p=FALSE) {

  observe <- data$observe
  expect  <- data$expect

  expect.number <- lapply(1:length(expect), function(x) {
    expect[[x]][name]
  }) %>% unlist()

  mean <- mean(expect.number)
  sd   <- sd(expect.number)

  result <- data.frame(
    condition = name,
    observe   = observe[name],
    expect    = mean,
    z_score   = (observe[name]-mean)/sd,
    log2FC    = log2(observe[name]/mean),
    upper.p   = pnorm(observe[name], mean=mean, sd=sd, lower.tail=FALSE, log.p=log.p),
    lower.p   = pnorm(observe[name], mean=mean, sd=sd, lower.tail=TRUE,  log.p=log.p))

  return(result)

}

#' Calculate the p-value by comparing the observed value to an empirical model constructed 
#' from expected values with each range of distances (condition).
#' 
#' @param data The input object contained the information of observed and expected number.
#' @param log.p If TRUE, the log of p-value will be returned.
#' @param parallel If a number greater than 1 is assigned, the function will run in parallel.
#' 
#' @export
ObsExpSTAT <-  function(
  data, log.p=FALSE, parallel=1) {
  
  if (!all(names(data)%in%c("observe", "expect"))) {
    stop("Please assign the result from ObsExpObj")
  }

  if (!is.numeric(parallel)) {
    stop("Please assign the number of cores to run the process")
  }

  observe <- data$observe
  expect  <- data$expect

  ## Check the names of observe and expect are the same
  lapply(1:length(expect), function(x){

    if (!identical(names(observe), names(expect[[x]]))) {
    
      message("Names of observed and expected objects do not match.")

      if (!all(names(observe) %in% names(expect))) {
        stop("Please ensure all observed names are in the expected list.")
      }

    }

  })

  ## Get charaters of names
  name_chr <- names(observe)

  if (parallel==1) {

    result <- lapply(name_chr, function(x){

      ObsExpSTATbyName(data, name=x, log.p=log.p)

    }) %>% Reduce(rbind, .)

  } else if (parallel > 1) {

    # // regist parallel
    gc(verbose = FALSE)
    doParallel::registerDoParallel(parallel)
    # // get index for parallel (list)
    if (length(name_chr) < parallel) {
      split_n <- split(1:length(name_chr), 1:length(name_chr))
    } else {
      split_n <- split(1:length(name_chr), cut(1:length(name_chr), parallel))
    }

    # // apply parallel
    result <- foreach(n = split_n, .combine=rbind) %dopar% {
      lapply(n, function(x){
        ObsExpSTATbyName(data, name=name_chr[x], log.p=log.p)
      }) %>% Reduce(rbind, .)
    }
    doParallel::stopImplicitCluster()

  } else {
    stop("Please assign the number over 1 of cores to run the process")

  }

  return(result)

}

#' Randomly select a character from a character list.
#' 
#' @param list The input list contained characters.
#' @param seed A number used to initialize the random character generator in R, ensuring consistent results.
#' @param n Number of characters to select.
#' 
#' @export
randomFactor <- function(list, seed=1, n=NULL) {

  # if (!is.character(list)) {
  #   stop("Please assign a character vector to list")
  # }

  if (is.null(n)) {
    stop("Please assign a number to n")
  }

  set.seed(seed)
  return(sample(list, n, replace=TRUE)) ## true means could be sample will be put back to list

}

#' Determine significance by associating factors with elements in the complete list.
#' 
#' @param factor Characters determined as factors.
#' @param total The complete list is utilized to extract randomized factors for calculating significance.
#' @param element Characters determined as element.
#' @param random.num Times of randomization.
#' @param log.p If set to TRUE, the logarithm of the p-value will be returned.
#' @param parallel If a number greater than 1 is assigned, the function will run in parallel.
#' 
#' @export
TargetFactorSTAT <- function(
  factor, total=NULL, element=NULL, random.num=10000, 
  log.p=FALSE, parallel=1) {

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


  total_dat <- data.table::data.table(
    total   = total,
    feature = total %in% factor,
    target  = total %in% element)
  
  ## Associate factor with element
  observe.num     <- total_dat[feature == TRUE, sum(target)]
  feature_len_num <- total_dat[, sum(feature)]
  total_idx       <- total_dat[, .I]  

  ## Randomly select elements from total
  if (parallel==1) {
    
    expect.num <- sapply(
      1:random.num, 
      function(x) {
        set.seed(x)
        sampled_indices <- sample(total_idx, feature_len_num, replace = TRUE)
        total_dat[sampled_indices, sum(target)]
      }
    )
    
  } else {

    gc(verbose = FALSE)
    doParallel::registerDoParallel(parallel)
    if (random.num < parallel) {
      split_n <- split(1:random.num, 1:random.num)
    } else {
      split_n <- split(1:random.num, cut(1:random.num, parallel))
    }

    expect.num <- foreach(
      n = split_n, 
      .combine = c
    ) %dopar% {

      sapply(
        n, 
        function(x) {
          set.seed(x)
          sampled_indices <- sample(total_idx, feature_len_num, replace = TRUE)
          total_dat[sampled_indices, sum(target)]
        }
      )

    }
    doParallel::stopImplicitCluster()

  }

  mean <- mean(expect.num)
  sd   <- sd(expect.num)

  ## Calculate the statistic
  result <- data.frame(
    observe      = observe.num,
    expect       = mean,
    log2FC       = log2(observe.num/mean),
    z_score      = (observe.num-mean)/sd,
    upper_pval   = pnorm(
      observe.num, mean=mean, sd=sd, lower.tail=FALSE, log.p=log.p),
    lower_pval   = pnorm(
      observe.num, mean=mean, sd=sd, lower.tail=TRUE , log.p=log.p))

  return(result)
}


#' Associate two interval-based factors on the genome within elements, 
#' compute their correlation, and determine their statistical significance.
#' 
#' @param factorA A path to a bed file or an intervals information, determine as factorA.
#' @param factorA.min Minimum distance of factorA to element. 
#' If factorA falls within this range of distances, it will be included.
#' @param factorA.max Maximum distance of factorA to element. If factorA falls within this range of distances, it will be included.
#' @param factorB A path to a bed file or an intervals information, determine as factorB.
#' @param factorB.min Minimum distance of factorB to element. 
#' If factorB falls within this range of distances, it will be included.
#' @param factorB.max Maximum distance of factorB to element. 
#' If factorB falls within this range of distances, it will be included.
#' @param element A path to a bed file or an intervals information, determine as element.
#' @param random.num Times of randomization.
#' @param log.p If set to TRUE, the logarithm of the p-value will be returned.
#' @param parallel If a number greater than 1 is assigned, the function will run in parallel.
#' 
#' @export
twoFactorElementSTAT <- function(
  factorA=NULL, factorA.min=0, factorA.max=0, factorB=NULL, factorB.min=0, factorB.max=0, 
  element=NULL, random.num=10000, log.p=FALSE, parallel=1) {

  if (is.null(factorA)) {
    stop("Please assign a character to factorA")
  }

  if (is.null(factorB)) {
    stop("Please assign a character to factorB")
  }

  if (is.null(element)) {
    stop("Please assign a character to element")
  }

  ##  Check the number of min and max
  if (!all(is.numeric(factorA.min), is.numeric(factorA.max), is.numeric(factorB.min), is.numeric(factorB.max))) {
    stop("Please assign a number to factorA.min, factorA.max, factorB.min, factorB.max")
  }

  if (all(factorA.min==0, factorA.max==0)) {
    
    factorA.intersect <- TRUE
 
  } else {

    if (factorA.min > factorA.max) {
      stop("factorA.min should be smaller than factorA.max")
    } 

    factorA.intersect <- FALSE

  }

  if (all(factorB.min==0, factorB.max==0)) {

    factorB.intersect <- TRUE
  
  } else {

    if (factorB.min > factorB.max) {
      stop("factorB.min should be smaller than factorB.max")
    }

    factorB.intersect <- FALSE

  }

  if (is.character(element)) {

    if (file.exists(element)) {
      element <- valr::read_bed(element)
    } else {
      stop("Please assign a valid file path")
    }

  } else {

    if (!is.data.frame(element)) {
      stop("Please check the input element")
    }

  }

  element.list <- unique(data.frame(element)[,4])

  element.factorA <- FactorElementCorrelate(
    factor  = element, 
    element = factorA, 
    tag     = "A")

  element.factorB <- FactorElementCorrelate(
    factor  = element, 
    element = factorB, 
    tag     = "A")

  if (isTRUE(factorA.intersect)) {

    factorA.list <- unique(element.factorA[element.factorA[,4]==0,1])

  } else {

    factorA.list <- unique(element.factorA[
      element.factorA[,4]>factorA.min & element.factorA[,4]<=factorA.max,1])

  }

  if (isTRUE(factorB.intersect)) {

    factorB.list <- unique(element.factorB[element.factorB[,4]==0,1])

  } else {

    factorB.list <- unique(element.factorB[
      element.factorB[,4]>factorB.min & element.factorB[,4]<=factorB.max,1])

  }


  result <- TargetFactorSTAT(
    factor        = factorA.list, 
    total         = element.list, 
    element       = factorB.list, 
    random.num    = random.num, 
    log.p         = log.p, 
    parallel      = parallel)

  return(result)

}







