#' Input factors and associate to elements. Output data.frame conatined the information of factors to nearest elements with the distance.
#' 
#' @param factor input the factor data.frame or the path to bed file.
#' @param element input the element data.frame or the path to bed file.
#' @param strand if assign as TRUE means input ele,ent contain the strand information at column 6, and will consider the strand information in the analysis.
#' @param tag if assign as TRUE means output the tag information in the analysis.
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
    
    factor <- valr::bed_read(factor, n_fields=4)[,1:4] 
    
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
      
      element <- valr::bed_read(element, n_fields=6)[,1:6]
      
    } else {
      
      element <- valr::bed_read(element, n_fields=4)[,1:4]
      
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
  
  intersect <- table[table$distance==0,] %>% 
    { .[order(.$overlap, decreasing=TRUE),] }
  
  multi.overlap <- table(intersect$factor_name) %>% 
    { .[.>1] } %>% names()
  
  if (length(multi.overlap) > 0) {
    
    multi.info <- lapply(multi.overlap, function(x) 
      intersect[intersect$factor_name==x,] %>% { .[1,] }) %>%
      Reduce(rbind, .)
    
    intersect.info <- rbind(
      intersect[!intersect$factor_name%in%multi.overlap,], multi.info) 
    
  } else {
    
    intersect.info <- intersect
    
  }
  
  associate <- table[!table$factor_name%in%intersect.info$factor_name,]
  multi.associate <- table(associate$factor_name) %>% 
    { .[.>1] } %>% names()
  
  if (length(multi.associate) > 0) {
    
    multi.info <- lapply(multi.associate, function(x) 
      associate[associate$factor_name==x,] %>% { .[1,] }) %>%
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
    
    forward <- element[element$strand=="+","element_name"] %>%
      { result[result$element_name%in%.,] }
    reverse <- element[element$strand=="+","element_name"] %>%
      { result[result$element_name%in%.,] }
    
    forward[forward$distance==0,"annotation"] <- "overlap"
    forward[forward$distance <0,"annotation"] <- "upstream"
    forward[forward$distance >0,"annotation"] <- "downstream"  
    
    reverse[reverse$distance==0,"annotation"] <- "overlap"
    reverse[reverse$distance <0,"annotation"] <- "downstream"
    reverse[reverse$distance >0,"annotation"] <- "upstream" 
    
    result <- rbind(forward, reverse)
    
  } else {
    
    result[result$distance==0,"annotation"] <- "overlap"
    result[result$distance <0,"annotation"] <- "upstream"
    result[result$distance >0,"annotation"] <- "downstream" 
    
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

FactorShufCorrelate <- function(
  factor, element=element, strand=FALSE, tag=FALSE, genome=genome, incl=NULL, excl=NULL, seed=1) {

  if (is.character(factor)) {
    
    if (!file.exists(factor)) {
      stop("Check the factor file exsist in path")
    }
    
    factor <- valr::bed_read(factor, n_fields=4)[,1:4] 
    
  } else if (is.data.frame(factor)) {
    
    factor <- factor[,1:4]
    
  } else if (class(factor)[[1]]=="GRanges") {
    
    stop("Element format is GRanges, please convert it to data.frame with bed format")
    
  } else {
    stop("Check the factor format")
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

  

}


