

FactorElementCorObj <- function(
  element, factor=factor, dist=1000000, strand=FALSE, 
  enrichType="upstream", parallel=1, parallel.type="mclapply") {
  
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
    
    if (enrichType=="upstream") {
      
      plus.element  <- element[element$strand=="+",]
      minus.element <- element[element$strand=="-",]
      
      if (parallel > 1) {
        
        if (parallel.type=="mclapply") {
          
          plus.result <- parallel::mclapply(1:nrow(plus.element), function(x) {
            table <- valr::bed_closest(factor, plus.element[x,]) %>% 
              filter(abs(.dist) <= dist, .dist <= 0) %>% select(factor_name.x, .dist)
            number <- table$.dist
            names(number) <- table$factor_name.x
            number
          }, mc.cores=parallel)
          
          minus.result <- parallel::mclapply(1:nrow(minus.element), function(x) {
            table <- valr::bed_closest(factor, minus.element[x,]) %>% 
              filter(abs(.dist) <= dist, .dist >= 0) %>% select(factor_name.x, .dist)
            number <- table$.dist
            names(number) <- table$factor_name.x
            number
          }, mc.cores=parallel)
          
        } else if (parallel.type=="bplapply") {
          
          BiocParallel::register(BiocParallel::MulticoreParam(workers = parallel))
          plus.result <- BiocParallel::bplapply(1:nrow(plus.element), function(x) {
            table <- valr::bed_closest(factor, plus.element[x,]) %>% 
              filter(abs(.dist) <= dist, .dist <= 0) %>% select(factor_name.x, .dist)
            number <- table$.dist
            names(number) <- table$factor_name
            number
          })
          
          minus.result <- BiocParallel::bplapply(1:nrow(minus.element), function(x) {
            table <- valr::bed_closest(factor, minus.element[x,]) %>% 
              filter(abs(.dist) <= dist, .dist >= 0) %>% select(factor_name.x, .dist)
            number <- table$.dist
            names(number) <- table$factor_name
            number
          })
          BiocParallel::register(BiocParallel::SerialParam())
          
        } 
        
      } else {
        
        plus.result <- lapply(1:nrow(plus.element), function(x) {
          table <- valr::bed_closest(factor, plus.element[x,]) %>% 
            filter(abs(.dist) <= dist, .dist <= 0) %>% select(factor_name.x, .dist)
          number <- table$.dist
          names(number) <- table$factor_name
          number
        })
        
        minus.result <- lapply(1:nrow(minus.element), function(x) {
          table <- valr::bed_closest(factor, minus.element[x,]) %>% 
            filter(abs(.dist) <= dist, .dist >= 0) %>% select(factor_name.x, .dist)
          number <- table$.dist
          names(number) <- table$factor_name
          number
        })
        
      }
      
      result <- c(plus.result, minus.result)
      
    }
    
    
  } else {
    
    if (parallel > 1) {
      
      if (parallel.type=="mclapply") {
        
        result <- parallel::mclapply(1:nrow(element), function(x) {
          table <- valr::bed_closest(factor, element[x,]) %>% 
            filter(abs(.dist) <= dist) %>% select(factor_name.x, .dist) %>% data.frame()
          number <- table$.dist
          names(number) <- table$factor_name.x
          number
        }, mc.cores=parallel)
        
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

# ExpFactorElementCor <- function(
#   element, factor=factor, dist=1000000, strand=FALSE, enrichType="upstream", parrallel=FALSE, 
#   seed=1, genome=genome, incl=NULL, excl=NULL) {

  


# }



# factorEnirchElement <- function(
#     element, factor=factor,)



# factorElementCor <- function(
#   element, factor=factor, dist=1000000, strand=FALSE, enrichType="upstream", parrallel=TRUE) {
  
#     ## Check factor format and input
#   if (is.character(factor)) {
    
#     if (!file.exists(factor)) {
#       stop("Check the factor file exsist in path")
#     }
    
#     factor <- valr::read_bed(factor, n_fields=4)[,1:4] 
    
#   } else if (is.data.frame(factor)) {
    
#     factor <- factor[,1:4]
    
#   } else if (class(factor)[[1]]=="GRanges") {
    
#     stop("Element format is GRanges, please convert it to data.frame with bed format")
    
#   } else {
#     stop("Check the factor format")
#   }
  
#   colnames(factor) <- c("chrom", "start", "end", "factor_name")
  
#   ## Check element format and input
#   if (is.character(element)) {
    
#     if (!file.exists(element)) {
#       stop("Check the element file exsist in path")
#     }
    
#     if (isTRUE(strand)) {
      
#       element <- valr::read_bed(element, n_fields=6)[,1:6]
      
#     } else {
      
#       element <- valr::read_bed(element, n_fields=4)[,1:4]
      
#     }
    
#   } else if (is.data.frame(element)) {
    
#     if (isTRUE(strand)) {
      
#       element <- element[,1:6]
      
#     } else {
      
#       element <- element[,1:4]
      
#     }
    
#   } else if (class(element)[[1]]=="GRanges") {
    
#     stop("Element format is GRanges, please convert it to data.frame with bed format")
    
#   } else {
#     stop("Check the element format")
#   }
  
#   if (isTRUE(strand)) {
    
#     colnames(element) <- c("chrom", "start", "end", "element_name", "score", "strand")
    
#   } else {
    
#     colnames(element) <- c("chrom", "start", "end", "element_name")
    
#   }

#   if (isTRUE(strand)) {

#     if (enrichType=="upstream") {

#       plus.element  <- element[element$strand=="+",]
#       minus.element <- element[element$strand=="-",]

#       if (isTRUE(parrallel)) {

#         plus.result <- BiocParallel::bplapply(1:nrow(plus.element), function(x) {
#           table <- valr::bed_closest(factor, plus.element[x,]) %>% 
#             filter(abs(.dist) <= dist, .dist <= 0) %>% select(factor_name, .dist)
#           number <- table$.dist
#           names(number) <- table$factor_name
#           number
#         })

#         minus.result <- BiocParallel::bplapply(1:nrow(minus.element), function(x) {
#           table <- valr::bed_closest(factor, minus.element[x,]) %>% 
#             filter(abs(.dist) <= dist, .dist >= 0) %>% select(factor_name, .dist)
#           number <- table$.dist
#           names(number) <- table$factor_name
#           number
#         })

#       } else {

#         plus.result <- lapply(1:nrow(plus.element), function(x) {
#           table <- valr::bed_closest(factor, plus.element[x,]) %>% 
#             filter(abs(.dist) <= dist, .dist <= 0) %>% select(factor_name, .dist)
#           number <- table$.dist
#           names(number) <- table$factor_name
#           number
#         })

#         minus.result <- lapply(1:nrow(minus.element), function(x) {
#           table <- valr::bed_closest(factor, minus.element[x,]) %>% 
#             filter(abs(.dist) <= dist, .dist >= 0) %>% select(factor_name, .dist)
#           number <- table$.dist
#           names(number) <- table$factor_name
#           number
#         })

#       }

#       result <- c(plus.result, minus.result)

#     }


#   } else {

#     if (isTRUE(parrallel)) {

#       result <- BiocParallel::bplapply(1:nrow(element), function(x) {
#         table <- valr::bed_closest(factor, element[x,]) %>% 
#           filter(abs(.dist) <= dist) %>% select(name.x, .dist)
#         number <- table$.dist
#         names(number) <- table$name.x
#         number
#       })

#     } else {

#       result <- lapply(1:nrow(element), function(x) {
#         table <- valr::bed_closest(factor, element[x,]) %>% 
#           filter(abs(.dist) <= dist) %>% select(name.x, .dist)
#         number <- table$.dist
#         names(number) <- table$name.x
#         number
#       })

#     }
#   }

#   return(result)
  
# }