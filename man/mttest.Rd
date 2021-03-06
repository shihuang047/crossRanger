% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/BetweenGroup.test.R
\name{mttest}
\alias{mttest}
\title{mttest}
\usage{
mttest(x, y, p.adj.method = "bonferroni", paired = FALSE)
}
\arguments{
\item{x}{A data.martrix or data.frame including multiple numeric vectors.}

\item{y}{A factor with two or more levels.}

\item{p.adj.method}{A string indicating the p-value correction method.}

\item{paired}{A logical indicating if paired between-group comparison is desired.}
}
\description{
Perform the univariate test for all features in the compositional microbiome data.
}
\examples{
y <-factor(c(rep("A", 30), rep("B", 30)))
y <-factor(c(rep("A", 15), rep("B", 15), rep("C", 15), rep("D", 15)))
x <- data.frame(rbind(t(rmultinom(7, 75, c(.201,.5,.02,.18,.099))),
            t(rmultinom(8, 75, c(.201,.4,.12,.18,.099))),
            t(rmultinom(15, 75, c(.011,.3,.22,.18,.289))),
            t(rmultinom(15, 75, c(.091,.2,.32,.18,.209))),
            t(rmultinom(15, 75, c(.001,.1,.42,.18,.299)))))
mttest(x, y, clr_transform=FALSE)
mttest(x, y, clr_transform=TRUE)
}
\author{
Shi Huang
}
