
#' Calculate the p-value of observed number of peaks in the expect distribution
#' 
#' @param data The input data contained the information of expect distribution or observe number.
#' @param name The name of condition to be calculated.
#' @param log.p If TRUE, the log of p-value will be returned.
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
    log2FC    = log2(observe[name]/mean),
    upper.p   = pnorm(observe[name], mean=mean, sd=sd, lower.tail=FALSE, log.p=log.p),
    lower.p   = pnorm(observe[name], mean=mean, sd=sd, lower.tail=TRUE,  log.p=log.p))

  return(result)

}

#' Calculate the p-value of observed number of peaks in the expect distribution
#' 
#' @param data The input data contained the information of expect distribution or observe number.
#' @param log.p If TRUE, the log of p-value will be returned.
#' @param parallel If assign number > 1, the function will run in parallel
#' @param parallel.type  Could be specify one of: \cr
#' \cr
#' "mclapply" - Use mclapply to run in parallel\cr
#' \cr
#' "bplapply" - Use BiocParallel to run in parallel 
#' 
#' @export
ObsExpSTAT <-  function(
  data, log.p=FALSE, parallel=1, parallel.type="mclapply") {

  if (!all(names(data)%in%c("observe", "expect"))) {
    stop("Please assign the result from ObsExpObj")
  }

  if (!is.numeric(parallel)) {
    stop("Please assign the number of cores to run the process")
  }

  observe <- data$observe
  expect  <- data$expect

  if (parallel==1) {
    
    result <- lapply(names(observe), function(x) {
      numbers <-  lapply(1:length(expect), function(y)
        expect[[y]][x]) %>% unlist()
      mean <-  mean(numbers)
      sd   <-  sd(numbers)
      data.frame(
        condition = x,
        observe   = observe[x],
        expect    = mean,
        log2FC    = log2(observe[x]/mean),
        upper.p   = pnorm(observe[x], mean=mean, sd=sd, lower.tail=FALSE, log.p=log.p),
        lower.p   = pnorm(observe[x], mean=mean, sd=sd, lower.tail=TRUE , log.p=log.p))
    }) %>% Reduce(rbind, .)
  
  } else if (parallel > 1) {

    if (parallel.type=="mclapply") {

      result <- parallel::mclapply(names(observe), function(x) {
        numbers <-  lapply(1:length(expect), function(y)
          expect[[y]][x]) %>% unlist()
        mean <-  mean(numbers)
        sd   <-  sd(numbers)
        data.frame(
          condition = x,
          observe   = observe[x],
          expect    = mean,
          log2FC    = log2(observe[x]/mean),
          upper.p   = pnorm(observe[x], mean=mean, sd=sd, lower.tail=FALSE, log.p=log.p),
          lower.p   = pnorm(observe[x], mean=mean, sd=sd, lower.tail=TRUE , log.p=log.p))
      }, mc.cores=parallel) %>% Reduce(rbind, .)

    } else if (parallel.type=="bplapply") {

      result <- BiocParallel::bplapply(names(observe), function(x) {
        numbers <-  lapply(1:length(expect), function(y)
          expect[[y]][x]) %>% unlist()
        mean <-  mean(numbers)
        sd   <-  sd(numbers)
        data.frame(
          condition = x,
          observe   = observe[x],
          expect    = mean,
          log2FC    = log2(observe[x]/mean),
          upper.p   = pnorm(observe[x], mean=mean, sd=sd, lower.tail=FALSE, log.p=log.p),
          lower.p   = pnorm(observe[x], mean=mean, sd=sd, lower.tail=TRUE , log.p=log.p))
      }) %>% Reduce(rbind, .)

    } else {
      stop("Please assign parallel.type as mclapply or bplapply")
    }

  } else {
    stop("Please assign the number over 1 of cores to run the process")
  }

  # if (isTRUE(parallel)) {

  #   result <- BiocParallel::bplapply(names(observe), function(x) {
  #     numbers <-  lapply(1:length(expect), function(y)
  #       expect[[y]][x]) %>% unlist()
  #     mean <-  mean(numbers)
  #     sd   <-  sd(numbers)
  #     data.frame(
  #       condition = x,
  #       observe   = observe[x],
  #       expect    = mean,
  #       log2FC    = log2(observe[x]/mean),
  #       upper.p   = pnorm(observe[x], mean=mean, sd=sd, lower.tail=FALSE, log.p=log.p),
  #       lower.p   = pnorm(observe[x], mean=mean, sd=sd, lower.tail=TRUE , log.p=log.p))
  #   }) %>% Reduce(rbind, .)

  # } else {
  
  #   result <- lapply(names(observe), function(x) {
  #     numbers <-  lapply(1:length(expect), function(y)
  #       expect[[y]][x]) %>% unlist()
  #     mean <-  mean(numbers)
  #     sd   <-  sd(numbers)
  #     data.frame(
  #       condition = x,
  #       observe   = observe[x],
  #       expect    = mean,
  #       log2FC    = log2(observe[x]/mean),
  #       upper.p   = pnorm(observe[x], mean=mean, sd=sd, lower.tail=FALSE, log.p=log.p),
  #       lower.p   = pnorm(observe[x], mean=mean, sd=sd, lower.tail=TRUE , log.p=log.p))
  #   }) %>% Reduce(rbind, .)

  # }

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
#' @param parallel If TRUE, will use parallel to calculate.
#' @param parallel.type  Could be specify one of: \cr
#' \cr
#' "mclapply" - Use mclapply to run in parallel\cr
#' \cr
#' "bplapply" - Use BiocParallel to run in parallel
#' 
#' @export
TargetFactorSTAT <- function(
  factor, total=NULL, element=NULL, random.num=10000, 
  log.p=FALSE, parallel=1, parallel.type="mclapply") {

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
  if (parallel==1) {

    expect.num <- lapply(1:random.num, function(x) 
      sum(randomFactor(total, seed=x, n=total.factor.num)%in%element)) %>% unlist()

  } else if (parallel > 1) {

    if (parallel.type=="mclapply") {

      expect.num <- parallel::mclapply(1:random.num, function(x) 
        sum(randomFactor(
          total, seed=x, n=total.factor.num)%in%element), mc.cores=parallel) %>% unlist()

    } else if (parallel.type=="bplapply") {

      BiocParallel::register(BiocParallel::MulticoreParam(workers = parallel))
      expect.num <- BiocParallel::bplapply(1:random.num, function(x) 
        sum(randomFactor(total, seed=x, n=total.factor.num)%in%element)) %>% unlist()
      BiocParallel::register(BiocParallel::SerialParam())

    } else {
      stop("Please assign parallel.type as mclapply or bplapply")
    }

  } else {
    stop("Please assign the number over 1 of cores to run the process")
  }

  if (isTRUE(parallel)) {

    gc(verbose=FALSE)
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


# FactorElementSTAT <- function(
#   factor=NULL, factor.min=0, factor.max=0, 
#   element=NULL, name.list=NULL, random.num=10000, log.p=FALSE, parallel=FALSE) {

#   if (is.null(factor)) {
#     stop("Please assign a character to factor")
#   }


#   if (is.null(element)) {
#     stop("Please assign a character to element")
#   }

#   ##  Check the number of min and max
#   if (!all(is.numeric(factor.min), is.numeric(factor.max))) {
#     stop("Please assign a number to factorA.min, factorA.max, factorB.min, factorB.max")
#   }

#   if (all(factor.min==0, factor.max==0)) {
    
#     factor.intersect <- TRUE
 
#   } else {

#     if (factor.min > factor.max) {
#       stop("factorA.min should be smaller than factorA.max")
#     } 

#     factor.intersect <- FALSE

#   }

#   if (is.character(element)) {

#     if (file.exists(element)) {
#       element <- valr::read_bed(element)
#     } else {
#       stop("Please assign a valid file path")
#     }

#   } else {

#     if (!is.data.frame(element)) {
#       stop("Please check the input element")
#     }

#   }

#   element.list <- unique(data.frame(element)[,4])

#   element.factor <- FactorElementCorrelate(
#     factor  = element, 
#     element = factor, 
#     tag     = "A")

#   if (isTRUE(factor.intersect)) {

#     factor.list <- unique(element.factor[element.factor[,4]==0,1])

#   } else {

#     factor.list <- unique(element.factor[
#       element.factor[,4]>factor.min & element.factor[,4]<=factor.max,1])

#   }

#   result <- TargetFactorSTAT(
#     factor     = factor.list, 
#     total      = element.list, 
#     element    = name.list, 
#     random.num = random.num, 
#     log.p      = log.p, 
#     parallel  = parallel)

#   return(result)

# }




#' Associate two factors in elements and calculate the significance
#' 
#' @param factorA A bed information of factorA.
#' @param factorA.min Minimum distance of element to factorA.
#' @param factorA.max Maximum distance of element to factorA
#' @param factorB A bed information of factorB.
#' @param factorB.min Minimum distance of element to factorB.
#' @param factorB.max Maximum distance of element to factorB.
#' @param element A bed information of elements.
#' @param random.num Number of random times to get expect result.
#' @param log.p If TRUE, will log2 the p.value.
#' @param parallel If TRUE, will use parallel to calculate.
#' 
#' @export
twoFactorElementSTAT <- function(
  factorA=NULL, factorA.min=0, factorA.max=0, factorB=NULL, factorB.min=0, factorB.max=0, 
  element=NULL, random.num=10000, log.p=FALSE, parallel=1, parallel.type="mclapply") {

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
    parallel      = parallel,
    parallel.type = parallel.type)

  return(result)

}







