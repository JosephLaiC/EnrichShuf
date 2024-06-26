% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CompilePeak.R
\name{CountNumber}
\alias{CountNumber}
\title{Count the number of elements associated with factors within the specified distance.}
\usage{
CountNumber(numbers, dist = 1e+06, intersect = FALSE, include = "all")
}
\arguments{
\item{numbers}{The data included factors associated with each element, along with their distances.}

\item{dist}{Distance to include associating factors to each element.}

\item{intersect}{If set to TRUE, results will include the factor intersect with elements.}

\item{include}{Could be specified one of the: \cr
\cr
"all" - Include all factors associated with elements. \cr
\cr
"upstream" - Include all factors associated with elements at upstream. \cr
\cr
"downstream" -  Include all factors associated with elements at downstream.}
}
\description{
Count the number of elements associated with factors within the specified distance.
}
