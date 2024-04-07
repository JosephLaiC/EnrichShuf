
#' Associate the factors to each element with specified distance, and compile the information to an object.
#' 
#' @param element Input the element data.frame or the path to bed file.
#' @param factor Input the factor data.frame or the path to the bed file.
#' @param dist Distance information for associating factors to each element.
#' @param strand If set to TRUE, it means that the input element contains strand information in column 6, 
#' and the analysis will take the strand information into consideration.
#' @param parallel If a number greater than 1 is assigned, the function will run in parallel.
#' @param parallel.type Could be specify one of: \cr
#' \cr
#' "mclapply" - Apply malapply to perform the run in parallel.\cr
#' \cr
#' "bplapply" - Apply bplapply to perfrom the run in parallel.
#' 
#' @export
FactorElementCorObj <- function(
  element, factor=factor, dist=1000000, strand=FALSE,
  parallel=1, parallel.type="mclapply") {
  
  if (any(!is.numeric(parallel))) {
    stop("Check the parallel input is numeric")
  }

  if (!parallel > 0) {
    stop("Check the parallel input is larger than 0")
  }

  if (!is.numeric(dist)) {
    stop("Check the dist input is numeric")
  }

  ## Check factor format and input
  if (is.character(factor)) {
    
    if (!file.exists(factor)) {
      stop("Check the factor file exsist in path")
    }
    
    factor <- valr::read_bed(factor, n_fields=4)[,1:4] 
    
  } else if (is.data.frame(factor)) {
    
    factor <- factor[,1:4]
    
  } else if (class(factor)[[1]]=="GRanges") {
    
    stop("Element format is GRanges, please convert it to data.frame with bed format")
    
  } else {
    stop("Check the factor format")
  }
  
  colnames(factor) <- c("chrom", "start", "end", "factor_name")
  
  ## Check element format and input
  if (is.character(element)) {
    
    if (!file.exists(element)) {
      stop("Check the element file exsist in path")
    }
    
    if (isTRUE(strand)) {
      
      element <- valr::read_bed(element, n_fields=6)[,1:6]
      
    } else {
      
      element <- valr::read_bed(element, n_fields=4)[,1:4]
      
    }
    
  } else if (is.data.frame(element)) {
    
    if (isTRUE(strand)) {
      
      element <- element[,1:6]
      
    } else {
      
      element <- element[,1:4]
      
    }
    
  } else if (class(element)[[1]]=="GRanges") {
    
    stop("Element format is GRanges, please convert it to data.frame with bed format")
    
  } else {
    stop("Check the element format")
  }
  
  if (isTRUE(strand)) {
    
    colnames(element) <- c("chrom", "start", "end", "element_name", "score", "strand")
    
  } else {
    
    colnames(element) <- c("chrom", "start", "end", "element_name")
    
  }

  gc(verbose=FALSE)
  if (isTRUE(strand)) {
      
    plus.element  <- element[element$strand=="+",]
    minus.element <- element[element$strand=="-",]
      
    if (parallel > 1) {
        
      ## Create the index for parallel
      plus.idx.num <- seq(1, nrow(plus.element), ceiling(nrow(plus.element)/parallel))
      plus.idx <- data.frame(
        start = plus.idx.num,
        end = c(plus.idx.num[-1]-1, nrow(plus.element)))

      minus.idx.num <- seq(1, nrow(minus.element), ceiling(nrow(minus.element)/parallel))
      minus.idx <- data.frame(
        start = minus.idx.num,
        end = c(minus.idx.num[-1]-1, nrow(minus.element)))

      if (parallel.type=="mclapply") {
          
        plus.result <- parallel::mclapply(1:nrow(plus.idx), function(x) {
          start <- plus.idx[x,1]; end <- plus.idx[x,2]
          lapply(start:end, function(y){

            table <- valr::bed_closest(factor, plus.element[y,]) %>% 
              filter(abs(.dist) <= dist) %>% select(factor_name.x, .dist)
            number <- table$.dist
            names(number) <- table$factor_name.x
            number

          })
        }, mc.cores=parallel) %>% Reduce(c, .)
          

        minus.result <- parallel::mclapply(1:nrow(minus.idx), function(x) {
          start <- minus.idx[x,1]; end <- minus.idx[x,2]
          lapply(start:end, function(y){

            table <- valr::bed_closest(factor, minus.element[y,]) %>% 
              filter(abs(.dist) <= dist) %>% select(factor_name.x, .dist)
            number <- table$.dist * -1
            names(number) <- table$factor_name.x
            number

          })
        }, mc.cores=parallel) %>% Reduce(c, .)
        
      } else if (parallel.type=="bplapply") {
          
        BiocParallel::register(BiocParallel::MulticoreParam(workers = parallel))
        plus.result <- BiocParallel::bplapply(1:nrow(plus.element), function(x) {
          table <- valr::bed_closest(factor, plus.element[x,]) %>% 
            filter(abs(.dist) <= dist) %>% select(factor_name.x, .dist)
          number <- table$.dist
          names(number) <- table$factor_name
          number
        })
          
        minus.result <- BiocParallel::bplapply(1:nrow(minus.element), function(x) {
          table <- valr::bed_closest(factor, minus.element[x,]) %>% 
            filter(abs(.dist) <= dist) %>% select(factor_name.x, .dist)
          number <- table$.dist
          names(number) <- table$factor_name
          number
        })
        BiocParallel::register(BiocParallel::SerialParam())
          
      } 
        
    } else if (parallel==1) {
        
      plus.result <- lapply(1:nrow(plus.element), function(x) {
        table <- valr::bed_closest(factor, plus.element[x,]) %>% 
          filter(abs(.dist) <= dist) %>% select(factor_name.x, .dist)
        number <- table$.dist
        names(number) <- table$factor_name
        number
      })
        
      minus.result <- lapply(1:nrow(minus.element), function(x) {
        table <- valr::bed_closest(factor, minus.element[x,]) %>% 
          filter(abs(.dist) <= dist) %>% select(factor_name.x, .dist)
        number <- table$.dist * -1
        names(number) <- table$factor_name
        number
      })
        
    }
      
    result <- c(plus.result, minus.result)
    
  } else if (!isTRUE(strand)) {
    
    if (parallel > 1) {
      
      idx.num <- seq(1, nrow(element), ceiling(nrow(element)/parallel))
      idx <- data.frame(
        start = idx.num,
        end = c(idx.num[-1]-1, nrow(element)))

      if (parallel.type=="mclapply") {
        
        result <- parallel::mclapply(1:nrow(idx), function(x) {
          start <- idx[x,1]; end <- idx[x,2]
          lapply(start:end, function(y){

            table <- valr::bed_closest(factor, element[y,]) %>% 
              filter(abs(.dist) <= dist) %>% select(factor_name.x, .dist)
            number <- table$.dist
            names(number) <- table$factor_name.x
            number

          })
        }, mc.cores=parallel) %>% Reduce(c, .)
        
        
      } else if (parallel.type=="bplapply") {
        
        BiocParallel::register(BiocParallel::MulticoreParam(workers = parallel))
        result <- BiocParallel::bplapply(1:nrow(element), function(x) {
          table <- valr::bed_closest(factor, element[x,]) %>% 
            filter(abs(.dist) <= dist) %>% select(factor_name.x, .dist) %>% data.frame()
          number <- table[,1]
          names(number) <- table[,2]
          number
        })
        BiocParallel::register(BiocParallel::SerialParam())
        
      }

      
    } else {
      
      result <- lapply(1:nrow(element), function(x) {
        table <- valr::bed_closest(factor, element[x,]) %>% 
          filter(abs(.dist) <= dist) %>% select(factor_name.x, .dist)
        number <- table$.dist
        names(number) <- table$factor_name.x
        number
      })
      
    }
  }
  
  names(result) <- data.frame(element)[,4]
  return(result)
  
}

#' Shuffle factors across the genome, then associate them to each element with specified distance, 
#' and compile the information to an object.
#' 
#' @param element Input the element data.frame or the path to bed file.
#' @param factor Input the factor data.frame or the path to the bed file.
#' @param dist Distance information for associating factors to each element.
#' @param strand If set to TRUE, it means that the input element contains strand information in column 6, 
#' and the analysis will take the strand information into consideration.
#' @param parallel If a number greater than 1 is assigned, the function will run in parallel.
#' @param parallel.type Could be specify one of: \cr
#' \cr
#' "mclapply" - Apply malapply to perform the run in parallel.\cr
#' \cr
#' "bplapply" - Apply bplapply to perfrom the run in parallel.
#' @param genome This parameter specifies the data or file path containing the names and sizes of chromosomes or contigs. 
#' Each name-size pair should be listed on a separate line and delimited by a tab.
#' @param incl The interval information to include the input regions.
#' @param excl The interval information to exclude the input regions.
#' @param seed 	A number used to initialize the random character generator in R, ensuring consistent results.
#' 
#' @export
ShufFactorElementCorObj <- function(
  element, factor=factor, dist=1000000, strand=FALSE,
  parallel=1, parallel.type="mclapply",
  genome=genome, incl=NULL, excl=NULL, seed=1) {

  if (is.character(factor)) {
    
    if (!file.exists(factor)) {
      stop("Check the factor file exsist in path")
    }
    
    factor <- valr::read_bed(factor, n_fields=4)[,1:4] 
    
  } else if (is.data.frame(factor)) {
    
    factor <- factor[,1:4]
    
  } else if (class(factor)[[1]]=="GRanges") {
    
    stop("Element format is GRanges, please convert it to data.frame with bed format")
    
  } else {
    stop("Check the factor format")
  }
  
  colnames(factor) <- c("chrom", "start", "end", "factor_name")

  if (is.character(genome)) {

    if (!file.exists(genome)) {
      stop("Check the genome file exsist in path")
    }

    genome <- valr::read_genome(genome)

  } else if (is.data.frame(genome)) {

    genome <- genome[,1:2]

  } else {
    stop("Check the genome format")
  }

  if (all(!is.null(incl), !is.null(excl))) {
    
    stop("Check the incl and excl, only could specify one")

  } else if (!is.null(excl)) {

    incl    <- data.frame(
      chrom = data.frame(genome)[,1],
      start = 0,
      end   = data.frame(genome)[,2]) %>% valr::bed_subtract(valr::read_bed(excl, n_fields=3))
    shuffle <- valr::bed_shuffle(factor, genome, seed=seed, incl=incl)

  } else if (!is.null(incl)) {

    shuffle <- valr::bed_shuffle(factor, genome, seed=seed, incl=incl)

  } else {

    shuffle <- valr::bed_shuffle(factor, genome, seed=seed)

  }

  result <- FactorElementCorObj(
    element       = element, 
    factor        = shuffle, 
    dist          = dist, 
    strand        = strand, 
    parallel      = parallel, 
    parallel.type = parallel.type)

  return(result)

}

#' Count the number of elements associated with factors within the specified distance.
#'
#' @param numbers The data included factors associated with each element, along with their distances.
#' @param dist Distance to include associating factors to each element.
#' @param intersect If set to TRUE, results will include the factor intersect with elements.
#' @param include Could be specified one of the: \cr
#' \cr
#' "all" - Include all factors associated with elements. \cr
#' \cr
#' "upstream" - Include all factors associated with elements at upstream. \cr
#' \cr
#' "downstream" -  Include all factors associated with elements at downstream.
#'
#' @export
CountNumber <- function(
  numbers, dist=1000000, intersect=FALSE, include="all") {

  if (is.null(numbers)) {
    return(0)
  }

  if (dist==0) {

    result <- lapply(numbers, function(x){
      sum(abs(x) == dist)
    })

  } else {
    
    if (isTRUE(intersect)) {

      if (include=="all") {

        result <- lapply(numbers, function(x){
          sum(abs(x) < dist)
        })

      } else if (include=="upstream") {

        result <- lapply(numbers, function(x){
          sum(x <= 0 & abs(x) < dist)
        })

      } else if (include=="downstream") {

        result <- lapply(numbers, function(x){
          sum(x >= 0 & abs(x) < dist)
        })

      } else {
        stop("Check the include parameter")
      }

    } else {

      if (include=="all") {

        result <- lapply(numbers, function(x){
          sum(abs(x) <= dist & abs(x) > 0)
        })

      } else if (include=="upstream") {

        result <- lapply(numbers, function(x){
          sum(x < 0 & abs(x) <= dist)
        })

      } else if (include=="downstream") {

        result <- lapply(numbers, function(x){
          sum(x > 0 & abs(x) <= dist)
        })

      } else {
        stop("Check the include parameter")
      }

    }

  }

  return(result)

}

#' Transform the data, which includes associated factors for each element along with their distances, into a data list.
#' 
#' @param data The data included factors associated with each element, along with their distances.
#' @param dist Distance to include associating factors to each element.
#' @param intersect If set to TRUE, results will include the factor intersect with elements.
#' @param include Could be specified one of the: \cr
#' \cr
#' "all" - Include all factors associated with elements. \cr
#' \cr
#' "upstream" - Include all factors associated with elements at upstream. \cr
#' \cr
#' "downstream" -  Include all factors associated with elements at downstream.
#' 
#' @export
CompileInfo <- function(data, dist=1000000, intersect=FALSE, include="all") {

  result <- lapply(data, function(x){
    CountNumber(x, dist=dist, intersect=intersect, include=include)
  }) %>% unlist()

  names(result) <- names(data)
  return(result)

}

#' Compare the observed compilation information with the expected compilation information using a binomial distribution.
#' 
#' @param observe The observe compile information.
#' @param expect.data The expect compile information. Typically will be a list.
#' @param parallel If a number greater than 1 is assigned, the function will run in parallel.
#' @param parallel.type Could be specify one of: \cr
#' \cr
#' "mclapply" - Apply malapply to perform the run in parallel.\cr
#' \cr
#' "bplapply" - Apply bplapply to perfrom the run in parallel.
#'
#' @export
binomialPeakCompile <- function(
  observe, expect.data=expect.data, 
  parallel=1, parallel.type="mclapply") {

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

  increase.logic <- lapply(expect.data, function(x){observe > x})
  decrease.logic <- lapply(expect.data, function(x){observe < x})

  if (parallel > 1) {

    if (parallel.type=="mclapply") {

      idx.num <- seq(1, length(observe), ceiling(length(observe)/parallel))
      idx <- data.frame(
        start = idx.num,
        end = c(idx.num[-1]-1, length(observe)))

      gc(verbose = FALSE)
      increase.number <- parallel::mclapply(1:nrow(idx), function(x){

        start <- idx[x,1]; end <- idx[x,2]
        lapply(start:end, function(y){
          lapply(increase.logic, function(z){z[[y]]}) %>% unlist() %>% sum()
        })
      }, mc.cores = parallel) %>% unlist()

      gc(verbose = FALSE)
      decrease.number <- parallel::mclapply(1:nrow(idx), function(x){

        start <- idx[x,1]; end <- idx[x,2]
        lapply(start:end, function(y){
          lapply(decrease.logic, function(z){z[[y]]}) %>% unlist() %>% sum()
        })
      }, mc.cores = parallel) %>% unlist()

    }

  } else {

    increase.number <- lapply(1:length(observe), function(x){
        
        lapply(increase.logic, function(y){y[[x]]}) %>% unlist() %>% sum()
  
    }) %>% unlist()

    decrease.number <- lapply(1:length(observe), function(x){
        
        lapply(decrease.logic, function(y){y[[x]]}) %>% unlist() %>% sum()
  
    }) %>% unlist()

  }

  if (parallel > 1) {

    idx.num <- seq(1, length(observe), ceiling(length(observe)/parallel))
    idx <- data.frame(
      start = idx.num,
      end = c(idx.num[-1]-1, length(observe)))
    
    if (parallel.type=="mclapply") {

      gc(verbose = FALSE)
      shuffle.mean <- parallel::mclapply(1:nrow(idx), function(x){

        start <- idx[x,1]; end <- idx[x,2]
        lapply(start:end, function(y){
          lapply(expect.data, function(z){z[y]}) %>% unlist() %>% mean()
        }) %>% unlist()

      }, mc.cores = parallel) %>% unlist()

      gc(verbose = FALSE)
      result <- parallel::mclapply(1:nrow(idx), function(x){

        start <- idx[x,1]; end <- idx[x,2]
        lapply(start:end, function(y){
          data.frame(
            name    = names(observe)[y],
            observe = observe[y],
            expect  = shuffle.mean[y],
            log2FC  = log2(observe[y]/shuffle.mean[y]),
            upper.p = binom.test(
              increase.number[y], length(expect.data), 0.5, alternative="greater")$p.value,
            lower.p = binom.test(
              decrease.number[y], length(expect.data), 0.5, alternative="greater")$p.value)
        }) %>% Reduce(rbind, .)

      }, mc.cores = parallel) %>% Reduce(rbind, .)

    }

  } else {

    shuffle.mean <- lapply(1:length(observe), function(x){
        
        lapply(expect.data, function(y){y[x]}) %>% unlist() %>% mean()
  
    }) %>% unlist()

    result <- lapply(1:length(observe), function(x){
        
        data.frame(
          name    = names(observe)[x],
          observe = observe[x],
          expect  = shuffle.mean[x],
          log2FC  = log2(observe[x]/shuffle.mean[x]),
          upper.p = binom.test(
            increase.number[x], length(expect.data), 0.5, alternative="greater")$p.value,
          lower.p = binom.test(
            decrease.number[x], length(expect.data), 0.5, alternative="greater")$p.value)
  
    }) %>% Reduce(rbind, .)

  }

  return(result)

}



