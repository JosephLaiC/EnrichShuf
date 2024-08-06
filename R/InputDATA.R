#' Input factors and associate them with elements. The output data.frame contains information 
#' about factors and their nearest elements with the distance.
#' 
#' @param factor Input the factor data.frame or the path to the bed file.
#' @param element Input the element data.frame or the path to bed file.
#' @param strand If set to TRUE, it means that the input element contains strand information 
#' in column 6, and the analysis will take the strand information into consideration.
#' @param tag If a character is assigned, the tag information will be outputted in the final table column.
#' @param outloc The location of the output file.
#' 
#' @export 
FactorElementCorrelate <- function(
    factor, element=element, strand=FALSE, tag=FALSE, outloc=NULL) {
  
  ## strand=TRUE only applied in element
  
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
  
  if (!isTRUE(strand)) {
    
    table <- valr::bed_closest(factor, element) %>%
      select(factor_name.x, element_name.y, .dist, .overlap) %>%
      setNames(c("factor_name", "element_name", "distance", "overlap")) %>%
      mutate(abs_distance = abs(distance)) %>%
      data.table::as.data.table()
    
    data.table::setorder(table, factor_name, -overlap, abs_distance, element_name)
    
  } else {
    
    table <- valr::bed_closest(factor, element) %>%
      select(factor_name.x, element_name.y, .dist, .overlap, strand.y) %>%
      setNames(c("factor_name", "element_name", "distance", "overlap", "strand")) %>%
      mutate(abs_distance = abs(distance)) %>%
      data.table::as.data.table()
    
    table[, distance := data.table::fifelse(strand == "-", -distance, distance)]
    data.table::setorder(table, factor_name, -overlap, abs_distance, element_name)
    
  }
  
  result <- table[, .SD[1], by = factor_name]
  # result <- result %>%
  #   dplyr::as_tibble() %>%
  #   select(factor_name, element_name, distance) %>%
  #   setNames(c("name", "tag", "distance"))
    
  
  return(result)
  
}

#' The input factors are shuffled and associated with elements. 
#' The output data.frame contains information about the random shuffle factors and their nearest elements with the distance.
#' 
#' @param factor Input the factor data.frame or the path to the bed file.
#' @param element Input the element data.frame or the path to bed file.
#' @param strand 	If set to TRUE, it means that the input element contains strand information in column 6, 
#' and the analysis will take the strand information into consideration.
#' @param tag If a character is assigned, the tag information will be outputted in the final table column.
#' @param outloc The location of the output file.
#' @param genome This parameter specifies the data or file path containing the names and sizes of chromosomes or contigs. 
#' Each name-size pair should be listed on a separate line and delimited by a tab.
#' @param incl The interval information to include the input regions.
#' @param excl The interval information to exclude the input regions.
#' @param seed A number used to initialize the random character generator in R, ensuring consistent results.
#' 
#' @export 
FactorShufCorrelate <- function(
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

  # Check incl and excl
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
    shuffle <- valr::bed_shuffle(factor, genome, seed=seed, incl=incl)

  } else if (!is.null(incl)) {

    shuffle <- valr::bed_shuffle(factor, genome, seed=seed, incl=incl)

  } else {

    shuffle <- valr::bed_shuffle(factor, genome, seed=seed)

  }

  result <- FactorElementCorrelate(
    shuffle, element=element, strand=strand, tag=tag, outloc=outloc)
  return(result)

}


#' Input the data processed by `RegionAnnoFromData`` and tally the occurrences for each condition 
#' (e.g., intersections or specific distance ranges).
#'
#' @param data Input data.frame or compressed/uncompressed file. 
#' Data should be in the format outputted by RegionAnnoFromData . It contain the information of: \cr
#' \cr
#' column 1: Name of the interval for the input factor. \cr
#' column 2: Name of the interval for the input element. \cr
#' column 3: Relative position information of factors with respect to the nearest element. \cr
#' column 4: Distance from the factors to the nearest element.
#' @param intersect If set to TRUE, the output result will place the number of intersecting peaks 
#' (where column 4's "distance" is equal to 0) in the first row.
#' @param condition A list of two numbers separated by a hyphen ("-"). Input the annotated data and tally the 
#' occurrences for each condition (e.g., intersections or specific distance ranges).
#' 
#' @export
CountCorrelation <- function(
    data, intersect=TRUE,
    condition=c("0-3000", "3000-10000", "10000-20000", "20000-30000", "30000-40000", "40000-50000")){

  #library(readr); library(stringr)

  if (is.data.frame(data)){
    data <- data
  } else if (is.character(data)){
    data <- data.frame(readr::read_tsv(data, show_col_types=FALSE))
  } else {
    stop("Check the input data format, could be dataframe or file path")
  }

  result <- NULL

  if (isTRUE(intersect)){
    number <- length(which(data[,4] == 0))
    names(number) <- "intersect"
    result <- c(result, number)
  }

  result <- c(result, unlist(lapply(condition, function(x){
    
    condition.list <- stringr::str_split(x, "-")[[1]] %>% as.numeric()
    
    if (any(!length(condition.list)==2, is.na(condition.list))){
      
      message("error with condition of ", x, " , ignore..")
      return(NA)

    } else {
      
      MIN <- min(condition.list)
      MAX <- max(condition.list)
      number <- length(which(data[,4] > MIN& data[,4] <= MAX))
      names(number) <- paste(MIN, MAX, sep="_")
      number

    }
  }
  )))

  return(result)

}


#' Given the input of bin size, minimum, and maximum values, generate intervals from the minimum to the maximum in specified bin increments. 
#' This string could be used in counting the occurrences of factors with associated elements.
#'
#' @param bin Interval size of each window..
#' @param min Minimum number of conditions
#' @param max Maximum number of conditions
#' @param type Could be specify: \cr
#' \cr
#' "continue" - Output a string of intervals with the range of distances between 0 and a series of consecutive numbers. \cr
#' \cr
#' "eachBin"   - Output a string of intervals based on specified bin sizes, dividing the range of distances into equal segments.
#'
#' @export
BinsDefine <- function(bin=NULL, min=NULL, max=NULL, type=NULL){

  if (any(is.null(bin), is.null(min), is.null(max), is.null(type))){
    stop("Please specify all parameters, bin,min,max,type.")
  }

  if (!all(is.numeric(bin), is.numeric(min), is.numeric(max), is.character(type))){
    stop("Check the format of input parameter.")
  }

  num <- seq(min, max, bin)

  if (type=="continue"){
    condition <- NULL
    for (i in 1:(length(num)-1)){
      condition <- c(condition, paste(num[i], num[i+1], sep="-"))
    }
  } else if (type=="within"){
    condition <- paste(num[1], num[-1], sep = "-")
  } else {
    stop("Please specify the type as continue or whithin")
  }

  return(condition)

}

#' Input the correlated (annotation) data
#'
#' Input the correlated (annotation) data and count the number of each condition by bins. (e.g. intersect or with specific distances)
#'
#' @param data Input data.frame or compressed/uncompressed file. 
#' Data should be in the format outputted by RegionAnnoFromData . It contain the information of: \cr
#' \cr
#' column 1: Name of the interval for the input factor. \cr
#' column 2: Name of the interval for the input element. \cr
#' column 3: Relative position information of factors with respect to the nearest element. \cr
#' column 4: Distance from the factors to the nearest element.
#' @param intersect If set to TRUE, the output result will place the number of intersecting peaks 
#' (where column 4's "distance" is equal to 0) in the first row.
#' @param bin Interval size of each window..
#' @param min Minimum number of conditions
#' @param max Maximum number of conditions
#' @param type Could be specify: \cr
#' \cr
#' "continue" - Output a string of intervals with the range of distances between 0 and a series of consecutive numbers. \cr
#' \cr
#' "eachBin"   - Output a string of intervals based on specified bin sizes, dividing the range of distances into equal segments.
#'
#' @export
CountCorrelationByBin <- function(
    data, intersect=TRUE, bin=1000, min=0, max=1000000, count.type="within"){

  condition <- BinsDefine(bin=bin, min=min, max=max, type=count.type)
  result    <- CountCorrelation(data, intersect=intersect, condition=condition)
  return(result)

}

#' Input the data processed by `RegionAnnoFromData` and tally the occurrences for each condition 
#' (e.g., intersections or specific distance ranges) build by `BinsDefine`.
#' 
#' @param factor Input the factor data.frame or the path to the bed file.
#' @param element Input the element data.frame or the path to bed file.
#' @param strand 	If set to TRUE, it means that the input element contains strand information in column 6, 
#' and the analysis will take the strand information into consideration.
#' @param tag If a character is assigned, the tag information will be outputted in the final table column.
#' @param outloc The location of the output file.
#' @param genome This parameter specifies the data or file path containing the names and sizes of chromosomes or contigs. 
#' Each name-size pair should be listed on a separate line and delimited by a tab.
#' @param incl The interval information to include the input regions.
#' @param excl The interval information to exclude the input regions.
#' @param random.n Times of randomization.
#' @param intersect If set to TRUE, the output result will place the number of intersecting peaks 
#' (where column 4's "distance" is equal to 0) in the first row.
#' @param condition 	A list of two numbers separated by a hyphen ("-"). Input the annotated data and tally the occurrences for each condition 
#' (e.g., intersections or specific distance ranges).
#' @param parallel If a number greater than 1 is assigned, the function will run in parallel.
#' 
#' @importFrom foreach %dopar% foreach
#' @importFrom magrittr %>%
#' 
#' @export
ObsExpObj <- function(
  factor, element=element, strand=FALSE, tag=FALSE, outloc=NULL, 
  genome=genome, incl=NULL, excl=NULL, random.n=10000, intersect=TRUE,
  condition=c("0-3000", "3000-10000", "10000-20000", "20000-30000", "30000-40000", "40000-50000"),
  parallel=1) {

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
    factor = factor, element = element, strand = strand, tag = tag) %>%
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

  # Check incl and excl
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

  }

  if (parallel==1) {

    expect <- lapply(1:random.n, function(x)
      FactorShufCorrelate(
        factor = factor, element = element, strand = strand, tag = tag, outloc = NULL, 
        genome = genome, incl = incl, seed = x) %>% 
      CountCorrelation(intersect = intersect, condition = condition))

  } else if (parallel > 1) {
    
    gc(verbose = FALSE)
    doParallel::registerDoParallel(parallel)
    if (random.n < parallel) {
      split_n <- split(1:random.n, 1:random.n)
    } else {
      split_n <- split(1:random.n, cut(1:random.n, parallel))
    }


    expect <- foreach(n = split_n, .packages = "magrittr", .combine=c) %dopar% {
      lapply(n, function(x) {
        FactorShufCorrelate(
          factor = factor, element = element, strand = strand, 
          tag = tag, outloc = outloc, genome = genome, incl = incl, seed = x) %>% 
          CountCorrelation(intersect = intersect, condition = condition)
      })
    }
    doParallel::stopImplicitCluster()

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


#' 
#' Create the object contained the information of observe and expect by each bin.
#' 
#' @param factor Input the factor data.frame or the path to the bed file.
#' @param element Input the element data.frame or the path to bed file.
#' @param strand 	If set to TRUE, it means that the input element contains strand information in column 6, 
#' and the analysis will take the strand information into consideration.
#' @param tag If a character is assigned, the tag information will be outputted in the final table column.
#' @param outloc The location of the output file.
#' @param genome This parameter specifies the data or file path containing the names and sizes of chromosomes or contigs. 
#' Each name-size pair should be listed on a separate line and delimited by a tab.
#' @param incl The interval information to include the input regions.
#' @param excl The interval information to exclude the input regions.
#' @param random.n Times of randomization.
#' @param intersect If set to TRUE, the output result will place the number of intersecting peaks 
#' (where column 4's "distance" is equal to 0) in the first row.
#' @param bin Interval size of each window..
#' @param min Minimum number of conditions
#' @param max Maximum number of conditions
#' @param type Could be specify: \cr
#' \cr
#' "continue" - Output a string of intervals with the range of distances between 0 and a series of consecutive numbers. \cr
#' \cr
#' "eachBin"   - Output a string of intervals based on specified bin sizes, dividing the range of distances into equal segments.
#' @param parallel If a number greater than 1 is assigned, the function will run in parallel.
#' 
#' @export
ObsExpObjEachBin <- function(
  factor, element=element, strand=FALSE, tag=FALSE, outloc=NULL, 
  genome=genome, incl=NULL, excl=NULL, random.n=10000, intersect=TRUE,
  bin=1000, min=0, max=1000000, count.type="within", parallel=FALSE) {

  condition <- BinsDefine(bin=bin, min=min, max=max, type=count.type)
  result    <- ObsExpObj(
    factor, element=element, strand=strand, tag=tag, outloc=outloc, 
    genome=genome, incl=incl, excl=excl, random.n=random.n, intersect=intersect,
    condition=condition, parallel=parallel)
  return(result)

}

