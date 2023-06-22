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
