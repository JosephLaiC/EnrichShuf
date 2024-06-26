% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/InputDATA.R
\name{ObsExpObj}
\alias{ObsExpObj}
\title{Input the data processed by `RegionAnnoFromData` and tally the occurrences for each condition 
(e.g., intersections or specific distance ranges) build by `BinsDefine`.}
\usage{
ObsExpObj(
  factor,
  element = element,
  strand = FALSE,
  tag = FALSE,
  outloc = NULL,
  genome = genome,
  incl = NULL,
  excl = NULL,
  random.n = 10000,
  intersect = TRUE,
  condition = c("0-3000", "3000-10000", "10000-20000", "20000-30000", "30000-40000",
    "40000-50000"),
  parallel = 1
)
}
\arguments{
\item{factor}{Input the factor data.frame or the path to the bed file.}

\item{element}{Input the element data.frame or the path to bed file.}

\item{strand}{If set to TRUE, it means that the input element contains strand information in column 6, 
and the analysis will take the strand information into consideration.}

\item{tag}{If a character is assigned, the tag information will be outputted in the final table column.}

\item{outloc}{The location of the output file.}

\item{genome}{This parameter specifies the data or file path containing the names and sizes of chromosomes or contigs. 
Each name-size pair should be listed on a separate line and delimited by a tab.}

\item{incl}{The interval information to include the input regions.}

\item{excl}{The interval information to exclude the input regions.}

\item{random.n}{Times of randomization.}

\item{intersect}{If set to TRUE, the output result will place the number of intersecting peaks 
(where column 4's "distance" is equal to 0) in the first row.}

\item{condition}{A list of two numbers separated by a hyphen ("-"). Input the annotated data and tally the occurrences for each condition 
(e.g., intersections or specific distance ranges).}

\item{parallel}{If a number greater than 1 is assigned, the function will run in parallel.}
}
\description{
Input the data processed by `RegionAnnoFromData` and tally the occurrences for each condition 
(e.g., intersections or specific distance ranges) build by `BinsDefine`.
}
