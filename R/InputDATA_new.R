
#' Input factors and associate to elements. Output data.frame conatined the information of factors to nearest elements with the distance.
#' 
#' @param factor input the factor data.frame or the path to bed file.
#' @param element input the element data.frame or the path to bed file.
#' @param strand if assign as TRUE means input ele,ent contain the strand information at column 6, and will consider the strand information in the analysis.
#' @param tag if assign as TRUE means output the tag information in the analysis.
#' @param outloc The location of output file.
#' @param genome The genome information, could be the path to genome file or data frame contained the genome information.
#' @param incl The interval information contained included regions.
#' @param excl The interval information contained excluded regions.
#' @param seed The seed number for shuffle.
#' 
#' @noRd
FactorShufCorrelate_new <- function(
  factor, element=element, strand=FALSE, tag=FALSE, outloc=NULL, 
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

  result <- FactorElementCorrelate(
    element, element=shuffle, strand=strand, tag=tag, outloc=outloc)
  return(result)

}


#' Create the object contained the information of observe and expect.
#' 
#' @param factor input the factor data frame or the path to bed file.
#' @param element input the element data frame or the path to bed file.
#' @param strand if assign as TRUE means input ele,ent contain the strand information at column 6, and will consider the strand information in the analysis.
#' @param tag if assign as TRUE means output the tag information in the analysis.
#' @param outloc The location of output file.
#' @param genome The genome information, could be the path to genome file or data frame contained the genome information.
#' @param incl The interval information contained included regions.
#' @param excl The interval information contained excluded regions.
#' @param random.n Times of shuffle.
#' @param intersect If assign as TRUE, result will contained intersect number.
#' @param condition Range of distance to nearest factor. Two number saperate by "-".
#' @param parallel If assign number > 1, the function will run in parallel
#' @param parallel.type Could be specify one of: \cr
#' \cr
#' "mclapply" - Use mclapply to run in parallel\cr
#' \cr
#' "bplapply" - Use BiocParallel to run in parallel
#' 
#' @noRd
ObsExpObj_new <- function(
  factor, element=element, strand=FALSE, tag=FALSE, outloc=NULL, 
  genome=genome, incl=NULL, excl=NULL, random.n=10000, intersect=TRUE,
  condition=c("0-3000", "3000-10000", "10000-20000", "20000-30000", "40000-50000"),
  parallel=1, parallel.type="mclapply") {

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

  ## observe result
  observe <- FactorElementCorrelate(
    factor = element, element = factor, strand = strand, tag = tag) %>%
    CountCorrelation(intersect = intersect, condition = condition)

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

  }

  if (parallel==1) {

    expect <- lapply(1:random.n, function(x)
      FactorShufCorrelate(
        factor = element, element = factor, strand = strand, tag = tag, outloc = outloc, 
        genome = genome, incl = incl, seed = x) %>% 
      CountCorrelation(intersect = intersect, condition = condition))

  } else if (parallel>1) {
    
    gc(verbose = FALSE)

    if (parallel.type=="mclapply") {

      expect <- parallel::mclapply(1:random.n, function(x) {
        FactorShufCorrelate(
          factor = element, element = factor, strand = strand, tag = tag, outloc = outloc, 
          genome = genome, incl = incl, seed = x) %>% 
        CountCorrelation(intersect = intersect, condition = condition)
      } ,mc.cores = parallel)

    } else if (parallel.type=="bplapply") {

      BiocParallel::register(BiocParallel::MulticoreParam(workers = parallel))
      expect <- BiocParallel::bplapply(1:random.n, function(x) {
        FactorShufCorrelate(
          factor = element, element = factor, strand = strand, tag = tag, outloc = outloc, 
          genome = genome, incl = incl, seed = x) %>% 
        CountCorrelation(intersect = intersect, condition = condition)
      })
      BiocParallel::register(BiocParallel::SerialParam())

    }
    
    

  } else {

    stop("Check the parallel parameter must be numeric")
  
  }

  result <- list(observe = observe, expect = expect)

  if (!is.null(outloc)) {

    if (tools::file_ext(outloc)=="gz") {

      saveRDS(result, outloc, compress="gzip")

    } else {

      saveRDS(result, outloc)

    }

  }

  return(result)

}

#' Create the object contained the information of observe and expect by each bin.
#' 
#' @param factor input the factor data frame or the path to bed file.
#' @param element input the element data frame or the path to bed file.
#' @param strand if assign as TRUE means input ele,ent contain the strand information at column 6, and will consider the strand information in the analysis.
#' @param tag if assign as TRUE means output the tag information in the analysis.
#' @param outloc The location of output file.
#' @param genome The genome information, could be the path to genome file or data frame contained the genome information.
#' @param incl The interval information contained included regions.
#' @param excl The interval information contained excluded regions.
#' @param random.n Times of shuffle.
#' @param intersect If assign as TRUE, result will contained intersect number.
#' @param bin Interval size of each window.
#' @param min Minimum number of conditions
#' @param max Maximum number of conditions
#' @param count.type Could be specify one of: \cr
#' \cr
#' "continue" - All conditions will continues from min+interval*(max/bin(n)-1) to min+interval(max/bin(n))\cr
#' \cr
#' "within"   - All conditions will start from 0 to each intervals {0+interval(max/bin(n))}
#' @param parallel If assign number > 1, the function will run in parallel
#' @param parallel.type Could be specify one of: \cr
#' \cr
#' "mclapply" - Use mclapply to run in parallel\cr
#' \cr
#' "bplapply" - Use BiocParallel to run in parallel 
#' 
#' @noRd
ObsExpObjEachBin_new <- function(
  factor, element=element, strand=FALSE, tag=FALSE, outloc=NULL, 
  genome=genome, incl=NULL, excl=NULL, random.n=10000, intersect=TRUE,
  bin=1000, min=0, max=1000000, count.type="within", parallel=FALSE, parallel.type="mclapply") {

  condition <- BinsDefine(bin=bin, min=min, max=max, type=count.type)
  result    <- ObsExpObj_new(
    factor, element=element, strand=strand, tag=tag, outloc=outloc, 
    genome=genome, incl=incl, excl=excl, random.n=random.n, intersect=intersect,
    condition=condition, parallel=parallel, parallel.type=parallel.type)
  return(result)

}
