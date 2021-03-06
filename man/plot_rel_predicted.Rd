% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ranger_reg_plot.R
\name{plot_rel_predicted}
\alias{plot_rel_predicted}
\title{plot_rel_predicted}
\usage{
plot_rel_predicted(
  relTrain_data,
  prefix = "train",
  target_field = "value",
  outdir = NULL
)
}
\arguments{
\item{relTrain_data}{The output dataframe of \code{calc_rel_predicted}.}

\item{prefix}{The prefix of a train dataset.}

\item{target_field}{A string indicating the target field in the metadata for regression.}

\item{outdir}{The output directory.}
}
\description{
Calculate the relative predicted values to the spline fit in the training data.
}
\examples{
set.seed(123)
train_x <- data.frame(rbind(t(rmultinom(7, 75, c(.201,.5,.02,.18,.099))),
            t(rmultinom(8, 75, c(.201,.4,.12,.18,.099))),
            t(rmultinom(15, 75, c(.011,.3,.22,.18,.289))),
            t(rmultinom(15, 75, c(.091,.2,.32,.18,.209))),
            t(rmultinom(15, 75, c(.001,.1,.42,.18,.299)))))
train_y<- 1:60
test_x <- data.frame(rbind(t(rmultinom(7, 75, c(.201,.5,.02,.18,.099))),
            t(rmultinom(8, 75, c(.201,.4,.12,.18,.099))),
            t(rmultinom(15, 75, c(.011,.3,.22,.18,.289))),
            t(rmultinom(15, 75, c(.091,.2,.32,.18,.209)))))
test_y<- 1:45
train_rf_model<-rf.out.of.bag(train_x, train_y)
predicted_test_y<-predict(train_rf_model$rf.model, test_x)$predictions
relTrain_data<-calc_rel_predicted(train_y, train_rf_model$predicted, test_y, predicted_test_y)
plot_rel_predicted(relTrain_data)
}
\author{
Shi Huang
}
