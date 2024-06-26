% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ChromDist.R
\name{ObsExpChrOb}
\alias{ObsExpChrOb}
\title{Input the data processed to calculate the observed and expected distribution of each chromosome.}
\usage{
ObsExpChrOb(
  factor,
  outloc = NULL,
  genome = genome,
  incl = NULL,
  excl = NULL,
  random.n = 10000,
  parallel = 1
)
}
\arguments{
\item{factor}{Input the factor data.frame or the path to the bed file.}

\item{genome}{This parameter specifies the data or file path containing the names and sizes of chromosomes or contigs. 
Each name-size pair should be listed on a separate line and delimited by a tab.}

\item{incl}{The interval information to include the input regions.}

\item{excl}{The interval information to exclude the input regions.}

\item{random.n}{Times of randomization.}

\item{parallel}{If a number greater than 1 is assigned, the function will run in parallel.}

\item{seed}{A number used to initialize the random character generator in R, ensuring consistent results.}
}
\description{
Input the data processed to calculate the observed and expected distribution of each chromosome.
}
