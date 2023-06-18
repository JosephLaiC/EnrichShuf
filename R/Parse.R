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
    table  <- data.frame(
      number = length(which(data[,4] == 0)),
      type   = "intersect")
    result <- rbind(result, table)
  }

  for (i in 1:length(condition)){

    condition.list <- stringr::str_split(condition[i], "-")[[1]] %>% as.numeric()

    if (any(!length(condition.list)==2, is.na(condition.list))){

      message("error with condition of ", condition[i], " , ignore..")
      next

    } else {

      MIN <- min(condition.list)
      MAX <- max(condition.list)

      table  <- data.frame(
        number = length(which(data[,4] > MIN& data[,4] <= MAX)),
        type   = paste(MIN, MAX, sep="_"))
      result <- rbind(result, table)

    }

  }; result$type <- factor(result$type, levels=result$type)

  return(result)

}



CountCorrelation_factor <- function(
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

  # for (i in 1:length(condition)){

  #   condition.list <- stringr::str_split(condition[i], "-")[[1]] %>% as.numeric()

  #   if (any(!length(condition.list)==2, is.na(condition.list))){

  #     message("error with condition of ", condition[i], " , ignore..")
  #     next

  #   } else {

  #     MIN <- min(condition.list)
  #     MAX <- max(condition.list)

  #     number <- length(which(data[,4] > MIN& data[,4] <= MAX))
  #     names(number) <- paste(MIN, MAX, sep="_")
  #     result <- c(result, number)

  #   }

  # }; result$type <- factor(result$type, levels=result$type)

  return(result)

}


#' Input shuffle peaks correlated (annotation) information as a list object.
#'
#' @param dir Directory of shuffle peak sets
#' @param parallel If set TRUE, will compile the shuffle peak sets with multiple threads
#' @param shuffle.n Number of peaks input into list. Reflect the flag (1~n) of shuffle_$flag.txt.gz
#' @param shuffle.prefix Name of prefix of peak sets.
#' @param file.ext Extension name of peak sets.
#' @param intersect If assign TRUE, output result will contained the number intersect peaks ($tag_dist==0) at first row
#' @param condition list of two number saperate with "-", will be used to count the number of factor A within the distance of numbers.
#'
#'
#' @export
compileCorrelation <- function(
    dir, parallel=FALSE, shuffle.n=10000, shuffle.prefix="shuffle_", file.ext=".txt.gz",
    intersect=TRUE, condition=c("0-3000", "3000-10000", "10000-20000", "20000-30000", "40000-50000")){

  #library(dplyr)

  file.list <- paste0(shuffle.prefix, rep(1:shuffle.n), file.ext)

  if (!all(file.list%in%list.files(dir))){
    stop("Check the intput files in expected (shuffle peaksets) directory")
  }

  if (isTRUE(parallel)){

    expect.list <- BiocParallel::bplapply(
      file.path(dir, file.list), CountCorrelation, intersect=intersect,
      condition=condition)


  } else {

    expect.list <- NULL
    for (i in 1:length(file.list)){
      expect.list <- CountCorrelation(
        file.path(dir, file.list[i]), condition=condition)
    }

  }

  return(expect.list)

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
    data, intersect=TRUE, bin=1000, min=0, max=1000000, count.type="continue"){

  condition <- BinsDefine(bin=bin, min=min, max=max, type=count.type)
  table     <- CountCorrelation(data, intersect=intersect, condition=condition)
  return(table)

}

CountCorrelationByBin_factor <- function(
    data, intersect=TRUE, bin=1000, min=0, max=1000000, count.type="continue"){

  condition <- BinsDefine(bin=bin, min=min, max=max, type=count.type)
  table     <- CountCorrelation_factor(data, intersect=intersect, condition=condition)
  return(table)

}





#' Input shuffle peaks correlated (annotation) information as a list object.
#'
#' @param dir Directory of shuffle peak sets
#' @param parallel If set TRUE, will compile the shuffle peak sets with multiple threads
#' @param shuffle.n Number of peaks input into list. Reflect the flag (1~n) of shuffle_$flag.txt.gz
#' @param shuffle.prefix Name of prefix of peak sets.
#' @param file.ext Extension name of peak sets.
#' @param intersect If assign TRUE, output result will contained the number intersect peaks ($tag_dist==0) at first row
#' @param condition list of two number saperate with "-", will be used to count the number of factor A within the distance of numbers.
#'
#' @export
compileCorrelationByBin <- function(
    dir, parallel=FALSE, shuffle.n=10000, shuffle.prefix="shuffle_", file.ext=".txt.gz",
    intersect=TRUE, bin=1000, min=0, max=1000000, count.type="continue"){

  condition <- BinsDefine(bin=bin, min=min, max=max, type=count.type)
  result    <- compileCorrelation(
    dir, parallel=parallel, shuffle.n=shuffle.n, shuffle.prefix=shuffle.prefix,
    file.ext=file.ext, intersect=intersect, condition=condition)
  return(result)
}
