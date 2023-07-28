#' Re parse and perform statistic the shuffle peak sets
#' @param data Input suject
#' @param data.type If assign "old", input data type is old version
#' @param from If assign "within", data type of input data
#' @param to Condition of output data
#' 
#' @export
ReconstructExpect <- function(
    data, from="within", data.type="old",
    to=c("intersect", "0-3000", "3000-10000", "10000-20000", "20000-30000", "40000-50000")) {
  
  if (data.type=="old"){
    
    stat.info <- data$statistic
    expect.numbers <- data$expectNum
    
    intersect.info   <- stat.info[which(stat.info$type=="intersect"),]
    intersect.expect <- expect.numbers[which(stat.info$type=="intersect")]
    
    within.info   <- stat.info[which(stat.info$type!="intersect"),]
    within.expect <- expect.numbers[which(stat.info$type!="intersect")]
    
    number.table <- stringr::str_split_fixed(within.info$type, "_", 2)
    colnames(number.table) <- c("min", "max")
    within.info <- cbind(within.info, number.table)
    
    condition.numbers <- stringr::str_split_fixed(to[!to%in%"intersect"], "-", 2)
    
    if (isTRUE("intersect"%in%to)) {
      
      observe.number <- intersect.info[,"observe.num"]
      expect.list    <- intersect.expect[[1]]
      
      SD   <- sd(expect.list)
      MEAN <- mean(expect.list)
      
      STAT.INFO <- data.frame(
        type        = "intersect",
        observe.num = observe.number,
        expect.mean = MEAN,
        log2FC      = log2(observe.number/MEAN), 
        upper.pVal  = pnorm(observe.number, mean=MEAN, sd=SD, lower.tail=FALSE),
        lower.pVal  = pnorm(observe.number, mean=MEAN, sd=SD, lower.tail=TRUE))
      
      EXPECT.NUM.LIST <- intersect.expect
      
    } else {
      STAT.INFO <- NULL; EXPECT.NUM.LIST <- NULL
    }
    
    STAT.INFO <- NULL; EXPECT.NUM.LIST <- NULL
    for (x in 1:nrow(condition.numbers)) {
      
      min <- as.numeric(condition.numbers[x,1])
      max <- as.numeric(condition.numbers[x,2])
      
      if (!min < max) {
        stop("Check the input condition numbers")
      }
      
      if (min==0) {
        
        observe.number <- within.info[which(within.info$max==max),"observe.num"]
        expect.list    <- within.expect[which(within.info$max==max)][[1]]
        
      } else {
        
        observe.min <- within.info[which(within.info$max==min),"observe.num"]
        observe.max <- within.info[which(within.info$max==max),"observe.num"]
        
        expect.min <- within.expect[which(within.info$max==min)][[1]]
        expect.max <- within.expect[which(within.info$max==max)][[1]]
        
        observe.number <- observe.max - observe.min
        expect.list    <- expect.max - expect.min
        
      }
      
      SD   <- sd(expect.list)
      MEAN <- mean(expect.list)
      
      RESULT <- data.frame(
        type        = paste0(min, "_", max),
        observe.num = observe.number,
        expect.mean = MEAN,
        log2FC      = log2(observe.number/MEAN), 
        upper.pVal  = pnorm(observe.number, mean=MEAN, sd=SD, lower.tail=FALSE),
        lower.pVal  = pnorm(observe.number, mean=MEAN, sd=SD, lower.tail=TRUE))
      
      STAT.INFO <- rbind(STAT.INFO, RESULT)
      EXPECT.NUM.LIST <- c(EXPECT.NUM.LIST, list(expect.list))
      
    }
    
    DATA <- list(statistic = STAT.INFO, expectNum = EXPECT.NUM.LIST)
    
    return(DATA)
  }
  
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
