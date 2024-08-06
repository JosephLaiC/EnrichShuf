


# FeatRegionChr <- function(
#   feature, region=region, genome=genome, parallel=1
# ) {
#   if (!is.data.frame(feature)) {
#     stop("Check the feature format")
#   }

#   if (!is.data.frame(region)) {
#     stop("Check the region format")
#   }

#   if (!is.data.frame(genome)) {
#     stop("Check the genome format")
#   }

#   feature <- feature[,1:4]
#   colnames(feature) <- c("chrom", "start", "end", "name")
  
#   ## feature processing
#   feature <- data.table::as.data.table(feature)
#   data.table::setkey(feature, start, end)


#   ## region processing
#   region <- region[,1:4] 
#   slop_dat <- valr::bed_slop(
#     region[,1:4], 
#     genome = genome, 
#     both   = 1000000,
#     trim   = TRUE)
#   colnames(slop_dat) <- c("chrom", "slop_start", "slop_end", "name")

#   region <- region %>% 
#     left_join(slop_dat, by=c("chrom", "name"))

#   region_num <- nrow(region)
  
#   if (parallel==1){

#     result <- lapply(1:region_num, function(x){

#       tmp_region <- region[x,]

#       chr     <- pull(tmp_region, chrom)
#       min_num <- pull(tmp_region, slop_start)
#       max_num <- pull(tmp_region, slop_end)

#       tmp_dat <- feature[chrom %in% chr & end >= min_num & start <= max_num] %>%
#         valr::bed_closest(., tmp_region)
  
#       tmp_res <- -pull(tmp_dat, .dist)
#       names(tmp_res) <- pull(tmp_dat, name.x)
#       tmp_res
  
#     })

#   } else if (parallel > 1) {

#     ## if parallel
#     gc(verbose = FALSE)
#     doParallel::registerDoParallel(parallel)
#     if (region_num < parallel) {
#       split_n <- split(1:region_num, 1:region_num)
#     } else {
#       split_n <- split(1:region_num, cut(1:region_num, parallel))
#     }

#     result <- foreach(n = split_n, .packages = "magrittr", .combine=c) {

#       lapply(n, function(x){

#         tmp_region <- region[x,]
        
#         chr     <- pull(tmp_region, chrom)
#         min_num <- pull(tmp_region, slop_start)
#         max_num <- pull(tmp_region, slop_end)

#         tmp_dat <- feature[chrom == chr & end >= min_num & start <= max_num] %>%
#           valr::bed_closest(., tmp_region)
  
#         tmp_res <- -pull(tmp_dat, .dist)
#         names(tmp_res) <- pull(tmp_dat, name.x)
#         tmp_res
  
#       })

#     }

#   }

#   names(result) <- pull(region, name)
#   return(result)

  
# }


# FeatRegionObj <- function(
#   feature, region=region, genome=genome, parallel=1
# ) {



# }







#' Associate the factors to each element with specified distance, and compile the information to an object.
#' 
#' @param element Input the element data.frame or the path to bed file.
#' @param factor Input the factor data.frame or the path to the bed file.
#' @param dist Distance information for associating factors to each element.
#' @param strand If set to TRUE, it means that the input element contains strand information in column 6, 
#' and the analysis will take the strand information into consideration.
#' @param parallel If a number greater than 1 is assigned, the function will run in parallel.
#' 
#' @export
FactorElementCorObj <- function(
  element, factor=factor, dist=1000000, strand=FALSE,
  parallel=1) {
  
  if (any(!is.numeric(parallel))) {
    stop("Check the parallel input is numeric")
  }
  if (length(parallel) > 1) {
    stop("Check the parallel input is a single number")
  }
  if (!parallel > 1) {
    stop("Check the parallel input is larger than 1")
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

  if (parallel == 1) {

    if (isTRUE(strand)) {
        
        ## // make the element to plus and minus strand
        plus.element  <- element[element$strand=="+",]
        minus.element <- element[element$strand=="-",]
        
        plus.result <- lapply(1:nrow(plus.element), function(x) {
          table <- valr::bed_closest(factor, plus.element[x,]) %>% 
            filter(abs(.dist) <= dist) %>% select(factor_name.x, .dist)
          number <- table$.dist
          names(number) <- table$factor_name.x
          number
        })
        
        minus.result <- lapply(1:nrow(minus.element), function(x) {
          table <- valr::bed_closest(factor, minus.element[x,]) %>% 
            filter(abs(.dist) <= dist) %>% select(factor_name.x, .dist)
          number <- table$.dist * -1
          names(number) <- table$factor_name.x
          number
        })
        
        result <- c(plus.result, minus.result)
        
      } else if (!isTRUE(strand)) {
        
        result <- lapply(1:nrow(element), function(x) {
          table <- valr::bed_closest(factor, element[x,]) %>% 
            filter(abs(.dist) <= dist) %>% select(factor_name.x, .dist)
          number <- table$.dist
          names(number) <- table$factor_name.x
          number
        })
    }

  } else {

    if (isTRUE(strand)) {

      ## // make the element to plus and minus strand
      plus.element  <- element[element$strand=="+",]
      minus.element <- element[element$strand=="-",]

      ## Create the index for parallel
      if (nrow(plus.element) < parallel) {
        split_n <- split(1:nrow(plus.element), 1:nrow(plus.element))
      } else {
        split_n <- split(1:nrow(plus.element), cut(1:nrow(plus.element), parallel))
      }

      if (nrow(minus.element) < parallel) {
        split_n_minus <- split(1:nrow(minus.element), 1:nrow(minus.element))
      } else {
        split_n_minus <- split(1:nrow(minus.element), cut(1:nrow(minus.element), parallel))
      }

      # Run the parallel
      gc(verbose = FALSE)
      doParallel::registerDoParallel(parallel)

      plus.result <- foreach(n = split_n, .combine=c) %dopar% {
        lapply(n, function(x) {
          table <- valr::bed_closest(factor, plus.element[x,]) %>% 
            filter(abs(.dist) <= dist) %>% select(factor_name.x, .dist)
          number <- table$.dist
          names(number) <- table$factor_name.x
          number
        })
      }
      doParallel::stopImplicitCluster()

      gc(verbose = FALSE)
      minus.result <- foreach(n = split_n_minus, .combine=c) %dopar% {
        lapply(n, function(x) {
          table <- valr::bed_closest(factor, minus.element[x,]) %>% 
            filter(abs(.dist) <= dist) %>% select(factor_name.x, .dist)
          number <- table$.dist * -1
          names(number) <- table$factor_name.x
          number
        })
      }
      doParallel::stopImplicitCluster()

      result <- c(plus.result, minus.result)

    } else {

      ## Create the index for parallel
      if (nrow(element) < parallel) {
        split_n <- split(1:nrow(element), 1:nrow(element))
      } else {
        split_n <- split(1:nrow(element), cut(1:nrow(element), parallel))
      }

      # Run the parallel
      gc(verbose = FALSE)
      doParallel::registerDoParallel(parallel)

      result <- foreach(n = split_n, .combine=c) %dopar% {
        lapply(n, function(x) {
          table <- valr::bed_closest(factor, element[x,]) %>% 
            filter(abs(.dist) <= dist) %>% select(factor_name.x, .dist)
          number <- table$.dist
          names(number) <- table$factor_name.x
          number
        })
      }
      doParallel::stopImplicitCluster()

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
#' @param genome This parameter specifies the data or file path containing the names and sizes of chromosomes or contigs. 
#' Each name-size pair should be listed on a separate line and delimited by a tab.
#' @param incl The interval information to include the input regions.
#' @param excl The interval information to exclude the input regions.
#' @param seed 	A number used to initialize the random character generator in R, ensuring consistent results.
#' 
#' @export
ShufFactorElementCorObj <- function(
  element, factor=factor, dist=1000000, strand=FALSE,
  parallel=1, genome=genome, incl=NULL, excl=NULL, seed=1) {

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
    parallel      = parallel)

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

        result <- sum(abs(numbers) < dist)

      } else if (include=="upstream") {

        result <- sum(numbers <= 0 & abs(numbers) < dist)

      } else if (include=="downstream") {

        result <- sum(numbers >= 0 & abs(numbers) < dist)

      } else {
        stop("Check the include parameter")
      }

    } else {

      if (include=="all") {

        result <- sum(abs(numbers) <= dist & abs(numbers) > 0)

      } else if (include=="upstream") {

        result <- sum(numbers < 0 & abs(numbers) <= dist)

      } else if (include=="downstream") {

        result <- sum(numbers > 0 & abs(numbers) <= dist)

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

  return(result)

}

#' Compile the shuffle objects.
#'
#' @param dir The directory path to the shuffle files.
#' @param shuffle_name The prefix of the shuffle files.
#' @param ext_file The extension of the shuffle files.
#' @param dist Distance information for associating factors to each element.
#' @param intersect If set to TRUE, results will include the factor intersect with elements.
#' @param include Could be specified one of the: \cr
#' \cr
#' "all" - Include all factors associated with elements. \cr
#' \cr
#' "upstream" - Include all factors associated with elements at upstream. \cr
#' \cr
#' "downstream" -  Include all factors associated with elements at downstream.
#' @param shuffle_times The number of shuffle files.
#' @param parallel If a number greater than 1 is assigned, the function will run in parallel.
#'
#' @export
shuffleCompile <- function(
  dir, shuffle_name = "shuffle_", ext_file = ".rds", 
  dist = NULL, intersect = FALSE, include = "all", 
  shuffle_times = 100, parallel = 1) {

  if (!is.character(dir)) {
    stop("Check the dir input")
  }
  if (!dir.exists(dir)) {
    stop("Check the dir exsist")
  }

  if (!is.character(shuffle_name)) {
    stop("Check the shuffle_name input")
  }

  if (!is.character(ext_file)) {
    stop("Check the ext_file input")
  }

  if (!is.numeric(dist)) {
    stop("Check the sig_dist input")
  }

  if (!is.numeric(shuffle_times)) {
    stop("Check the shuffle_times input")
  }

  if (any(!is.numeric(parallel))) {
    stop("Check the parallel input is numeric")
  }
  if (length(parallel) > 1) {
    stop("Check the parallel input is a single number")
  }
  if (!parallel > 1) {
    stop("Check the parallel input is larger than 1")
  }

  # Check the shuffle files
  lapply(
    1:shuffle_times, 
    function(x) {
      file <- file.path(dir, paste0(shuffle_name, x, ext_file))
      if (!file.exists(file)) {
        stop("Check the shuffle files:", file, "exsist in path")
      }
    }
  )

  gc(verbose=FALSE)
  if (parallel == 1) {

    result <- lapply(
      1:shuffle_nums, 
      function(x) {
        file <- file.path(dir, paste0(shuffle_name, x, ext_file))
        readRDS(file) %>%
          CompileInfo(
            dist      = dist,
            intersect = FALSE, 
            include   = "all"
          )
      }
    )

  } else {

    if (shuffle_times < parallel) {
      split_n <- split(1:shuffle_times, 1:shuffle_times)
    } else {
      split_n <- split(1:shuffle_times, cut(1:shuffle_times, parallel))
    }

    gc(verbose = FALSE)
    doParallel::registerDoParallel(parallel)
    result <- foreach(n = split_n, .combine=c) %dopar% {
      lapply(
        n, 
        function(x) {
          file <- file.path(dir, paste0(shuffle_name, x, ext_file))
          readRDS(file) %>%
            CompileInfo(
              dist      = dist,
              intersect = FALSE, 
              include   = "all"
            )
        }
      )
    }
    doParallel::stopImplicitCluster()

  }

  return(result)
  
}



#' Compare the observed compilation information with the expected compilation information using a binomial distribution.
#' 
#' @param observe The observe compile information.
#' @param expect.data The expect compile information. Typically will be a list.
#' @param p.adjust If set to TRUE, the function will adjust the p-value using the FDR method.
#' @param parallel If a number greater than 1 is assigned, the function will run in parallel.
#'
#' @export
binomialPeakCompile <- function(
  observe, expect.data=expect.data, p.adjust=TRUE, parallel=1) {

  # / Check the input
  if (!is.numeric(observe)) {
    stop("Check the observe data")
  }

  if (any(!is.numeric(parallel))) {
    stop("Check the parallel input is numeric")
  }
  if (length(parallel) > 1) {
    stop("Check the parallel input is a single number")
  }
  if (!parallel > 1) {
    stop("Check the parallel input is larger than 1")
  }
  
  # / Check the observe and expect data
  lapply(expect.data, function(x){
    
    if (!is.numeric(x)) {
      stop("Check the expect data")
    }
    
    if (!identical(names(x), names(observe))) {
      stop("Check the names of expect and observe data")
    }
    
  })

  # / get logic list
  increase.logic <- lapply(expect.data, function(x){observe > x})
  decrease.logic <- lapply(expect.data, function(x){observe < x})

  # / get the number of increase and decrease from the logic list
  if (parallel == 1) {
      
    increase.number <- lapply(1:length(observe), function(x){
          
      lapply(increase.logic, function(y){y[[x]]}) %>% unlist() %>% sum()
    
    }) %>% unlist()
  
    decrease.number <- lapply(1:length(observe), function(x){
          
      lapply(decrease.logic, function(y){y[[x]]}) %>% unlist() %>% sum()
    
    }) %>% unlist()

  } else {

    if (length(observe) < parallel) {
      split_n <- split(1:length(observe), 1:length(observe))
    } else {
      split_n <- split(1:length(observe), cut(1:length(observe), parallel))
    }

    gc(verbose = FALSE)
    doParallel::registerDoParallel(parallel)

    increase.number <- foreach(n = split_n, .combine=c) %dopar% {
      lapply(
        n, 
        function(x){
        lapply(
          increase.logic, 
          function(y){
            y[[x]]
          }
        ) %>% unlist() %>% sum()
        }
      )
    } %>% unlist()

    decrease.number <- foreach(n = split_n, .combine=c) %dopar% {
      lapply(
        n, 
        function(x){
        lapply(
          decrease.logic, 
          function(y){
            y[[x]]
          }
        ) %>% unlist() %>% sum()
        }
      )
    } %>% unlist()

  }

  # / calculate the shuffle mean for getting the log2FC and calculate the p-value by binomial test
  if (parallel == 1) {

    shuffle.mean <- lapply(1:length(observe), function(x){
        
      lapply(expect.data, function(y){y[x]}) %>% unlist() %>% mean()
  
    }) %>% unlist()

    result <- lapply(1:length(observe), function(x){
        
        data.frame(
          name    = names(observe)[x],
          observe = observe[x],
          expect  = shuffle.mean[x],
          log2FC  = log2(observe[x]+1) - log2(shuffle.mean[x]+1),
          upper.p = binom.test(
            increase.number[x], length(expect.data), 0.5, alternative="greater")$p.value,
          lower.p = binom.test(
            decrease.number[x], length(expect.data), 0.5, alternative="greater")$p.value)
  
    }) %>% Reduce(rbind, .)

  } else {

    if (length(observe) < parallel) {
      split_n <- split(1:length(observe), 1:length(observe))
    } else {
      split_n <- split(1:length(observe), cut(1:length(observe), parallel))
    }

    gc(verbose = FALSE)
    doParallel::registerDoParallel(parallel)

    shuffle.mean <- foreach(n = split_n, .combine=c) %dopar% {
      lapply(
        n, 
        function(x){
          lapply(expect.data, function(y){y[x]}) %>% unlist() %>% mean()
        }
      ) %>% unlist()
    }

    result <- foreach(n = split_n, .combine=rbind) %dopar% {
      lapply(
        n, 
        function(x){
          data.frame(
            name    = names(observe)[x],
            observe = observe[x],
            expect  = shuffle.mean[x],
            log2FC  = log2(observe[x]+1) - log2(shuffle.mean[x]+1),
            upper.p = binom.test(
              increase.number[x], length(expect.data), 0.5, alternative="greater")$p.value,
            lower.p = binom.test(
              decrease.number[x], length(expect.data), 0.5, alternative="greater")$p.value)
        }
      ) %>% Reduce(rbind, .)
    }

  }

  result$pval <- ifelse(result$log2FC > 0, result$upper.p, result$lower.p)

  if (isTRUE(p.adjust)) {
    result$upper.FDR <- p.adjust(result$upper.p, method="fdr")
    result$lower.FDR <- p.adjust(result$lower.p, method="fdr")
    result$FDR <- ifelse(result$log2FC > 0, result$upper.FDR, result$lower.FDR)
  }

  return(result)

}



