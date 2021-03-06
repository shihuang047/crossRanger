% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ranger_RF_util.R
\name{get.predict.probability.from.forest}
\alias{get.predict.probability.from.forest}
\title{get.predict.probability.from.forest}
\usage{
get.predict.probability.from.forest(model, newx)
}
\arguments{
\item{model}{A object of class \code{rf.out.of.bag} with elements.}

\item{newx}{The external test data.}
}
\value{
\code{probs}
}
\description{
to get votes from the forest in the external test. Last update: Apr. 7, 2019
}
\author{
Shi Huang
}
