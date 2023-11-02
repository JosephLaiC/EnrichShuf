#' Input factors and associate to elements. Output data.frame conatined the information of factors to nearest elements with the distance.
#' 
#' @param factor input the factor data.frame or the path to bed file.
#' @param element input the element data.frame or the path to bed file.
#' @param strand if assign as TRUE means input ele,ent contain the strand information at column 6, and will consider the strand information in the analysis.
#' @param tag if assign as TRUE means output the tag information in the analysis.
#' @param outloc The location of output file.
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
  
  table <- valr::bed_closest(factor, element)[
    ,c("factor_name.x", "element_name.y", ".dist", ".overlap")] %>% data.frame()
  colnames(table) <- c("factor_name", "element_name", "distance", "overlap")
  
  intersect <- table[which(table$distance==0),] %>% 
    { .[order(.$overlap, decreasing=TRUE),] }
  
  multi.overlap <- table(intersect$factor_name) %>% 
    { .[which(.>1)] } %>% names()
  
  if (length(multi.overlap) > 0) {
    
    multi.info <- lapply(multi.overlap, function(x) 
      intersect[which(intersect$factor_name==x),] %>% { .[1,] }) %>%
      Reduce(rbind, .)
    
    intersect.info <- rbind(
      intersect[!intersect$factor_name%in%multi.overlap,], multi.info) 
    
  } else {
    
    intersect.info <- intersect
    
  }
  
  associate <- table[!table$factor_name%in%intersect.info$factor_name,]
  multi.associate <- table(associate$factor_name) %>% 
    { .[which(.>1)] } %>% names()
  
  if (length(multi.associate) > 0) {
    
    multi.info <- lapply(multi.associate, function(x) 
      associate[which(associate$factor_name==x),] %>% { .[1,] }) %>%
      Reduce(rbind, .)
    
    associate.info <- rbind(
      associate[!associate$factor_name%in%multi.associate,], multi.info) 
    
  } else {
    
    associate.info <- associate
    
  }
  
  result <- rbind(intersect.info, associate.info)
  
  result$annotation <- NA
  
  if (isTRUE(strand)) {
    
    element <- data.frame(element)
    
    if (!all(unique(element$strand)%in%c("+", "-"))) {
      stop("Check the element strand information")
    }
    
    forward <- element[which(element$strand=="+"),"element_name"] %>%
      { result[result$element_name%in%.,] }
    reverse <- element[which(element$strand=="+"),"element_name"] %>%
      { result[result$element_name%in%.,] }
    
    forward[which(forward$distance==0),"annotation"] <- "overlap"
    forward[which(forward$distance <0),"annotation"] <- "upstream"
    forward[which(forward$distance >0),"annotation"] <- "downstream"  
    
    reverse[which(reverse$distance==0),"annotation"] <- "overlap"
    reverse[which(reverse$distance <0),"annotation"] <- "downstream"
    reverse[which(reverse$distance >0),"annotation"] <- "upstream" 
    
    result <- rbind(forward, reverse)
    
  } else {
    
    result[which(result$distance==0),"annotation"] <- "overlap"
    result[which(result$distance <0),"annotation"] <- "upstream"
    result[which(result$distance >0),"annotation"] <- "downstream" 
    
  }
  
  result <- result[,c("factor_name", "element_name", "annotation", "distance")]
  result$distance <- abs(result$distance)
  
  
  if (is.character(tag)) {
    
    colnames(result) <- c("name", paste(tag, c("name", "annotation", "distance"), sep="_")) 
    
  } else {
    
    colnames(result) <- c("name", "tag", "annotation", "distance")
    
  }
  
  if (!is.null(outloc)) {

    if (!is.character(outloc)) {
      stop("Check the output location")
    }

    write.table(
      result, outloc, sep="\t", quote=FALSE, row.names=FALSE)

  }

  return(result)
  
}

# Under developed
# ExpShuffleFactor <- function(
#   factor, tag=FALSE, outloc=NULL, 
#   genome=genome, incl=NULL, excl=NULL, seed=1) {

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

# }

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
    shuffle, element=element, strand=strand, tag=tag, outloc=outloc)
  return(result)

}


#' Input the correlated (annotation) data
#'
#' Input the correlated (annotation) data and count the number of each condition. (e.g. intersect or with specific distances)
#'
#' @param data Input data.frame or compressed/uncompressed data. Data should output from RegionAnnoFromData, contain the information of
#' name (1st column); $tag_region (2nd column); $tag_type (3rd column); $tag_dist (4th column)
#' @param intersect If assign TRUE, output result will contained the number intersect peaks ($tag_dist==0) at first row
#' @param condition list of two number saperate with "-", will be used to count the number of factor A within the distance of numbers.
#'
#' @export
CountCorrelation <- function(
    data, intersect=TRUE,
    condition=c("0-3000", "3000-10000", "10000-20000", "20000-30000", "40000-50000")){

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


#' Get the intervals from min to max by bins.
#'
#' @param bin Interval size of each window.
#' @param min Minimum number of conditions
#' @param max Maximum number of conditions
#' @param type Could be specify: \cr
#' \cr
#' "continue" - All conditions will continues from min+interval*(max/bin(n)-1) to min+interval(max/bin(n))\cr
#' \cr
#' "within"   - All conditions will start from 0 to each intervals {0+interval(max/bin(n))}
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
#' @param data Input data.frame or compressed/uncompressed data. Data should output from RegionAnnoFromData, contain the information of
#' name (1st column); $tag_region (2nd column); $tag_type (3rd column); $tag_dist (4th column)
#' @param intersect If assign TRUE, output result will contained the number intersect peaks ($tag_dist==0) at first row
#' @param bin Interval size of each window.
#' @param min Minimum number of conditions
#' @param max Maximum number of conditions
#' @param type Could be specify one of: \cr
#' \cr
#' "continue" - All conditions will continues from min+interval*(max/bin(n)-1) to min+interval(max/bin(n))\cr
#' \cr
#' "within"   - All conditions will start from 0 to each intervals {0+interval(max/bin(n))}
#'
#' @export
CountCorrelationByBin <- function(
    data, intersect=TRUE, bin=1000, min=0, max=1000000, count.type="within"){

  condition <- BinsDefine(bin=bin, min=min, max=max, type=count.type)
  result    <- CountCorrelation_factor(data, intersect=intersect, condition=condition)
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
#' @export
ObsExpObj <- function(
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
        factor = factor, element = element, strand = strand, tag = tag, outloc = outloc, 
        genome = genome, incl = incl, seed = x) %>% 
      CountCorrelation(intersect = intersect, condition = condition))

  } else if (parallel>1) {
    
    gc(verbose = FALSE)

    if (parallel.type=="mclapply") {

      expect <- parallel::mclapply(1:random.n, function(x) {
        FactorShufCorrelate(
          factor = factor, element = element, strand = strand, tag = tag, outloc = outloc, 
          genome = genome, incl = incl, seed = x) %>% 
        CountCorrelation(intersect = intersect, condition = condition)
      } ,mc.cores = parallel)

    } else if (parallel.type=="bplapply") {

      BiocParallel::register(BiocParallel::MulticoreParam(workers = parallel))
      expect <- BiocParallel::bplapply(1:random.n, function(x) {
        FactorShufCorrelate(
          factor = factor, element = element, strand = strand, tag = tag, outloc = outloc, 
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
#' @export
ObsExpObjEachBin <- function(
  factor, element=element, strand=FALSE, tag=FALSE, outloc=NULL, 
  genome=genome, incl=NULL, excl=NULL, random.n=10000, intersect=TRUE,
  bin=1000, min=0, max=1000000, count.type="within", parallel=FALSE, parallel.type="mclapply") {

  condition <- BinsDefine(bin=bin, min=min, max=max, type=count.type)
  result    <- ObsExpObj(
    factor, element=element, strand=strand, tag=tag, outloc=outloc, 
    genome=genome, incl=incl, excl=excl, random.n=random.n, intersect=intersect,
    condition=condition, parallel=parallel, parallel.type=parallel.type)
  return(result)

}

