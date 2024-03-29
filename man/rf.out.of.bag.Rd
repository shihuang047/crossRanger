% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ranger_RF_util.R
\name{rf.out.of.bag}
\alias{rf.out.of.bag}
\title{rf.out.of.bag}
\usage{
rf.out.of.bag(
  x,
  y,
  ntree = 500,
  verbose = FALSE,
  sparse = FALSE,
  imp_pvalues = FALSE
)
}
\arguments{
\item{x}{Training data: data.matrix or data.frame.}

\item{y}{A response vector. If a factor, classification is assumed, otherwise regression is assumed.}

\item{ntree}{The number of trees.}

\item{verbose}{A boolean value indicates if showing computation status and estimated runtime.}

\item{sparse}{A boolean value indicates if the input matrix transformed into sparse matrix for rf modeling.}

\item{imp_pvalues}{If compute both importance score and pvalue for each feature.}
}
\value{
Object of class \code{rf.out.of.bag} with elements including the ranger object and critical metrics for model evaluation.
}
\description{
It runs standard random forests with out-of-bag error estimation for both classification and regression using \code{ranger}.
This is merely a wrapper that extracts relevant info from \code{ranger} output.
}
\details{
Ranger is a fast implementation of random forests (Breiman 2001) or recursive partitioning,
particularly suited for high dimensional data.
Classification, regression, and survival forests are supported.
}
\examples{
set.seed(123)
x <- data.frame(rbind(t(rmultinom(7, 75, c(.201,.5,.02,.18,.099))),
            t(rmultinom(8, 75, c(.201,.4,.12,.18,.099))),
            t(rmultinom(15, 75, c(.011,.3,.22,.18,.289))),
            t(rmultinom(15, 75, c(.091,.2,.32,.18,.209))),
            t(rmultinom(15, 75, c(.001,.1,.42,.18,.299)))))
y<-factor(c(rep("A", 15), rep("B", 15), rep("C", 15), rep("D", 15)))
y<-factor(c(rep("A", 20), rep("B", 20), rep("C", 20)))
rf.out.of.bag(x, y, imp_pvalues=FALSE)
rf.out.of.bag(x, y, imp_pvalues=TRUE)
x_ <- data.frame(rbind(t(rmultinom(7, 7500, rep(c(.201,.5,.02,.18,.099), 1000))),
            t(rmultinom(8, 750, rep(c(.201,.4,.12,.18,.099), 1000))),
            t(rmultinom(15, 750, rep(c(.011,.3,.22,.18,.289), 1000))),
            t(rmultinom(15, 750, rep(c(.091,.2,.32,.18,.209), 1000))),
            t(rmultinom(15, 750, rep(c(.001,.1,.42,.18,.299), 1000)))))
y_<-factor(c(rep("A", 15), rep("B", 15), rep("C", 15), rep("D", 15)))
rf.out.of.bag(x_, y_, imp_pvalues=FALSE)
rf.out.of.bag(x, y_, imp_pvalues=TRUE)
y0<-factor(c(rep("old", 30), rep("young", 30)))
rf.out.of.bag(x, y0, imp_pvalues=FALSE)
rf.out.of.bag(x, y0, imp_pvalues=TRUE)
y<- 1:60
rf.out.of.bag(x, y, imp_pvalues=FALSE)
}
\references{
Wright, M. N. & Ziegler, A. (2017). ranger: A fast implementation of random forests for
high dimensional data in C++ and R. Journal of Statistical Software 77:1-17.
}
\seealso{
ranger
}
\author{
Shi Huang
}
