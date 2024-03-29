% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ranger_crossRF_util.R
\name{rf_reg.cross_appl}
\alias{rf_reg.cross_appl}
\title{rf_reg.cross_appl}
\usage{
rf_reg.cross_appl(rf_list, x_list, y_list)
}
\arguments{
\item{rf_list}{A list of rf.model objects from \code{rf.out.of.bag}.}

\item{x_list}{A list of training datasets usually in the format of data.frame.}

\item{y_list}{A list of responsive vector for regression in the training datasets.}
}
\value{
...
}
\description{
Based on pre-computed rf models regressing \code{c_category} in each the sub-datasets splited by the \code{s_category},
perform cross-datasets application of the rf models. The inputs are precalculated
rf regression models, x_list and y_list.
}
\examples{
df <- data.frame(rbind(t(rmultinom(14, 14*5, c(.21,.6,.12,.38,.099))),
            t(rmultinom(16, 16*5, c(.001,.6,.42,.58,.299))),
            t(rmultinom(30, 30*5, c(.011,.6,.22,.28,.289))),
            t(rmultinom(30, 30*5, c(.091,.6,.32,.18,.209))),
            t(rmultinom(30, 30*5, c(.001,.6,.42,.58,.299)))))
df0 <- data.frame(t(rmultinom(120, 600,c(.001,.6,.2,.3,.299))))
metadata<-data.frame(f_s=factor(c(rep("A", 60), rep("B", 60))),
                     f_s1=factor(c(rep(TRUE, 60), rep(FALSE, 60))),
                     f_c=factor(c(rep("C", 30), rep("H", 30), rep("D", 30), rep("P", 30))),
                     age=c(1:60, 2:61)
                     )

table(metadata[, 'f_c'])
reg_res<-rf_reg.by_datasets(df, metadata, s_category='f_c', c_category='age')
rf_reg.cross_appl(reg_res, x_list=reg_res$x_list, y_list=reg_res$y_list)
reg_res<-rf_reg.by_datasets(df, metadata, nfolds=5, s_category='f_c', c_category='age')
rf_reg.cross_appl(reg_res, x_list=reg_res$x_list, y_list=reg_res$y_list)
}
\seealso{
ranger
}
\author{
Shi Huang
}
