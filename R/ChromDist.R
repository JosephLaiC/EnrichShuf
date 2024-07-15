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
        valr::bed_shuffle(factor, genome, seed = x) %>%
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

    expect <- foreach(
      n = split_n, .packages = "magrittr", .combine=c
    ) %dopar% {

      lapply(
        n, 
        function(x) {
          valr::bed_shuffle(factor, genome, seed = x) %>%
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

#' Input the data processed to calculate the observed and expected values for each chromosome.
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
#' @param condition A list of two numbers separated by a hyphen ("-"). Input the annotated data and tally the occurrences for each condition 
#' (e.g., intersections or specific distance ranges).
#' @param parallel If a number greater than 1 is assigned, the function will run in parallel.
#' @param chrom_list A vector of chromosome names to be analyzed. If NULL, all chromosomes present in the genome parameter will be analyzed. This allows for selective analysis of specific chromosomes or contigs.
#' 
#' @importFrom foreach %dopar% foreach
#' @importFrom magrittr %>%
#' 
#' @export
ObsExpObjEachChrom <- function(
  factor, element=element, strand=FALSE, tag=FALSE, outloc=NULL, 
  genome=genome, incl=NULL, excl=NULL, random.n=10000, intersect=TRUE,
  condition=c("0-3000", "3000-10000", "10000-20000", "20000-30000", "30000-40000", "40000-50000"),
  parallel=1, chrom_list=NULL) {

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


  if (is.null(chrom_list)) {
    chrom_list <- pull(genome, chrom)
  }

  if (is.null(incl)) {
    
    result <- lapply(
      chrom_list,
      function(x){
        
        ObsExpObj(
          factor    = factor[factor$chrom==x,], 
          element   = element[element$chrom==x,],
          strand    = strand,
          tag       = tag,
          genome    = genome[genome$chrom==x,],
          random.n  = random.n,
          intersect = intersect,
          condition = condition,
          parallel  = parallel 
        )
        
      }
    ) %>% setNames(chrom_list)

  } else {

    result <- lapply(
      chrom_list,
      function(x){

        tmp_incl <- incl[incl$chrom==x,]

        if (nrow(tmp_incl)==0) {

          ObsExpObj(
            factor    = factor[factor$chrom==x,], 
            element   = element[element$chrom==x,],
            strand    = strand,
            tag       = tag,
            genome    = genome[genome$chrom==x,],
            random.n  = random.n,
            intersect = intersect,
            condition = condition,
            parallel  = parallel 
          )

        } else {

          ObsExpObj(
            factor    = factor[factor$chrom==x,], 
            element   = element[element$chrom==x,],
            strand    = strand,
            tag       = tag,
            genome    = genome[genome$chrom==x,],
            incl      = tmp_incl,
            random.n  = random.n,
            intersect = intersect,
            condition = condition,
            parallel  = parallel 
          )

        }
        
      }
    ) %>% setNames(chrom_list)

  }

  return(result)

} 

#' Input the data processed to calculate the correlation between features and regions of each chromosome for each bin.
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
#' @param chrom_list A vector of chromosome names to be analyzed. If NULL, all chromosomes present in the genome parameter will be analyzed. This allows for selective analysis of specific chromosomes or contigs.
#'
#' @export 
ObsExpObjEachBinEachChrom <- function(
  factor, element=element, strand=FALSE, tag=FALSE, outloc=NULL, 
  genome=genome, incl=NULL, excl=NULL, random.n=10000, intersect=TRUE,
  bin=1000, min=0, max=1000000, count.type="within", parallel=FALSE, chrom_list=NULL
) {

  condition <- BinsDefine(bin=bin, min=min, max=max, type=count.type)
  result    <- ObsExpObjEachChrom(
    factor, element=element, strand=strand, tag=tag, outloc=outloc, 
    genome=genome, incl=incl, excl=excl, random.n=random.n, intersect=intersect,
    condition=condition, parallel=parallel, chrom_list=chrom_list)
  return(result)

}

#' Input the data processed to calculate the correlation between features and regions for each chromosome.
#'
#' @param data Input object 
#' @param parallel If a number greater than 1 is assigned, the function will run in parallel.
#'
#' @export
ChromSTAT <- function(data, parallel=1) {

  result <- lapply(
    data,
    function(x) {

      ObsExpSTAT(x, parallel=parallel)

    }
  )

  return(result)

}

#' Function for the combination of the results across all chromosome.
#'
#' @param data Input object
#'
#' @export
CombineChrom <- function(data) {

  random.n <- length(data[[1]]$expect)

  result <- list(
    observe = {
      lapply(
        data, 
        function(x) {
          x$observe
        }
      ) %>% Reduce(`+`, .)
    },
    expect  = {
      lapply(
        1:random.n, 
        function(n){
            
          lapply(
            data, 
            function(x) {
              x$expect[[n]]
            }
          ) %>% Reduce(`+`, .) 
            
        }
      )
    }
  )

  return(result)

}





