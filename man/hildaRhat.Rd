% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/hilda_diagnosis.R
\name{hildaRhat}
\alias{hildaRhat}
\title{Output the maximum potential scale reduction statistic of all parameters
estimated}
\usage{
hildaRhat(jagsOutput)
}
\arguments{
\item{jagsOutput}{the output jags file generated by the jags function from
the R2jags package.}
}
\value{
a number for the Rhat statistic.
}
\description{
Output the maximum potential scale reduction statistic of all parameters
estimated
}
\examples{

inputFile <- system.file("extdata/hildaLocal.rdata", package="HiLDA")
hildaLocal <- readRDS(inputFile)
hildaRhat(hildaLocal)

}
