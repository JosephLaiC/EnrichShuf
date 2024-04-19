#' Generate shuffled features
#'
#' @param factor Input the factor data.frame or the path to the bed file.
#' @param genome This parameter specifies the data or file path containing the names and sizes of chromosomes or contigs. 
#' Each name-size pair should be listed on a separate line and delimited by a tab.
#' @param incl The interval information to include the input regions.
#' @param excl The interval information to exclude the input regions.
#' @param seed A number used to initialize the random character generator in R, ensuring consistent results.
#' 
#' @export 
shuffleGenerate <- function(
  factor, genome=genome, incl=NULL, excl=NULL, seed=1) {

  # / check input factor
  if (is.character(factor)) {
    
    if (!file.exists(factor)) {
      stop("Check the factor file exsist in path")
    }
    
    factor <- valr::read_bed(factor, n_fields=4)[,1:4] 
    
  } else if (is.data.frame(factor)) {
    
    factor <- factor[,1:4]
    
  } else if (class(factor)[[1]]=="GRanges") {
    
    stop("Factor format is GRanges, please convert it to data.frame with bed format")
    
  } else {
    stop("Check the factor format")
  }

  colnames(factor) <- c("chrom", "start", "end", "factor_name")

  # / check genome
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

  # / Check incl and excl
  if (is.character(incl)) {
    incl <- valr::read_bed(incl, n_fields=3)
  } else if (is.data.frame(incl)) {
    incl <- incl[,1:3]
  }

  if (is.character(excl)) {
    excl <- valr::read_bed(excl, n_fields=3)
  } else if (is.data.frame(excl)) {
    excl <- excl[,1:3]
  }

  if (all(!is.null(incl), !is.null(excl))) {
    
    stop("Check the incl and excl, only could specify one")

  } else if (!is.null(excl)) {

    incl    <- data.frame(
      chrom = data.frame(genome)[,1],
      start = 0,
      end   = data.frame(genome)[,2]) %>% valr::bed_subtract(excl)
    result <- valr::bed_shuffle(factor, genome, seed=seed, incl=incl)

  } else if (!is.null(incl)) {

    result <- valr::bed_shuffle(factor, genome, seed=seed, incl=incl)

  } else {

    result <- valr::bed_shuffle(factor, genome, seed=seed)

  }

  return(result)

}

#' Input the data processed to calculate the observed and expected distribution of each chromosome.
#' 
#' @param factor Input the factor data.frame or the path to the bed file.
#' @param genome This parameter specifies the data or file path containing the names and sizes of chromosomes or contigs. 
#' Each name-size pair should be listed on a separate line and delimited by a tab.
#' @param incl The interval information to include the input regions.
#' @param excl The interval information to exclude the input regions.
#' @param seed A number used to initialize the random character generator in R, ensuring consistent results.
#' @param random.n Times of randomization.
#' @param parallel If a number greater than 1 is assigned, the function will run in parallel.
#' 
#' @export 
ObsExpChrOb <- function(
  factor, outloc=NULL, 
  genome=genome, incl=NULL, excl=NULL, 
  random.n=10000, parallel=1) {

  # / check input factor
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

  # / check genome
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

  # / observe result
  observe <- factor %>%
    group_by(chrom) %>%
    summarise(n = n()) %>%
    { setNames(.$n, .$chrom) }

  # / expect result
  if (parallel == 1) {
    expect <- lapply(
      1:random.n, 
      function(x) {
        valr::bed_shuffle(factor, genomesize, seed = x) %>%
          group_by(chrom) %>%
          summarise(n = n()) %>%
          { setNames(.$n, .$chrom) }
      }
    )
  } else {

    gc(verbose = FALSE)
    doParallel::registerDoParallel(parallel)
    if (random.n < parallel) {
      split_n <- split(1:random.n, 1:random.n)
    } else {
      split_n <- split(1:random.n, cut(1:random.n, parallel))
    }

    expected_dat <- foreach(
      n = split_n, .packages = "magrittr", .combine=c
    ) %dopar% {

      lapply(
        n, 
        function(x) {
          valr::bed_shuffle(factor, genomesize, seed = x) %>%
            group_by(chrom) %>%
            summarise(n = n()) %>%
            { setNames(.$n, .$chrom) }
        }
      )

    }
    doParallel::stopImplicitCluster()
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