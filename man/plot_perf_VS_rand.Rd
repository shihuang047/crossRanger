% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ranger_reg_plot.R
\name{plot_perf_VS_rand}
\alias{plot_perf_VS_rand}
\title{plot_perf_VS_rand}
\usage{
plot_perf_VS_rand(
  x,
  y,
  predicted_y,
  prefix = "train",
  target_field,
  nfolds,
  metric = "MAE",
  permutation = 100,
  n_features = NA,
  outdir = NULL
)
}
\arguments{
\item{x}{The train data.}

\item{y}{The numeric labeling data.}

\item{predicted_y}{The predicted values for y.}

\item{prefix}{The prefix for the dataset in the training or testing.}

\item{target_field}{A string indicating the target field of the metadata for regression.}

\item{nfolds}{The number of folds in the cross validation. If nfolds > length(y)
or nfolds==-1, uses leave-one-out cross-validation. If nfolds was a factor, it means customized folds
(e.g., leave-one-group-out cv) were set for CV.}

\item{metric}{The regression performance metric applied, including MAE, RMSE, MSE, R_squared, Adj_R_squared, or Separman_rho.}

\item{permutation}{The permutation times for a random guess of regression performance.}

\item{n_features}{The number of features in the training data.}

\item{outdir}{The output directory.}
}
\description{
This outputs a histogram and a p-value showing if the performance of a real regression model
significantly better than null models.
}
\examples{
set.seed(123)
x <- data.frame(rbind(t(rmultinom(7, 75, c(.201,.5,.02,.18,.099))),
            t(rmultinom(8, 75, c(.201,.4,.12,.18,.099))),
            t(rmultinom(15, 75, c(.011,.3,.22,.18,.289))),
            t(rmultinom(15, 75, c(.091,.2,.32,.18,.209))),
            t(rmultinom(15, 75, c(.001,.1,.42,.18,.299)))))
y<- 1:60
rf_model<-rf.out.of.bag(x, y)
p<-plot_perf_VS_rand(x=x, y=y, predicted_y=rf_model$predicted, prefix="train", nfolds=5,
permutation=100, metric="MAE", target_field="age", n_features=5)
p
}
\author{
Shi Huang
}
