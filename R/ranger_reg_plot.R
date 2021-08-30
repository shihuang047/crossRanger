#' @importFrom stats cor wilcox.test residuals smooth.spline fitted median
#' @importFrom plyr ldply
#' @importFrom rlang .data
#' @import ggplot2
#' @importFrom gridExtra arrangeGrob
#' @importFrom reshape2 melt
#' @importFrom utils write.table


#' @title plot_obs_VS_pred
#' @description Plot a scatterplot of observed and predicted values from a ranger model.
#' @param y The numeric values for labeling data.
#' @param predicted_y The predicted values for y.
#' @param metric A regression performance metric (i.e., MAE, RMSE, MSE, R_squared, Adj_R_squared, or Separman_rho) showed in the scatter plot.
#' @param SampleIDs The sample ids in the data table.
#' @param prefix The prefix for the dataset in the training or testing.
#' @param target_field A string indicating the target field of metadata for regression.
#' @param span Controls the amount of smoothing for the default loess smoother in the scatterplot.
#' Smaller numbers produce wigglier lines, larger numbers produce smoother lines.
#' @param outdir The output directory.
#' @examples
#' set.seed(123)
#' x <- data.frame(rbind(t(rmultinom(7, 75, c(.201,.5,.02,.18,.099))),
#'             t(rmultinom(8, 75, c(.201,.4,.12,.18,.099))),
#'             t(rmultinom(15, 75, c(.011,.3,.22,.18,.289))),
#'             t(rmultinom(15, 75, c(.091,.2,.32,.18,.209))),
#'             t(rmultinom(15, 75, c(.001,.1,.42,.18,.299)))))
#' y<- 1:60
#' rf_model<-rf.out.of.bag(x, y)
#' plot_obs_VS_pred(y, predicted_y=rf_model$predicted, prefix="train", target_field="age", outdir=NULL)
#' @author Shi Huang
#' @export
plot_obs_VS_pred <- function(y, predicted_y, SampleIDs=NULL, prefix="train", target_field="value", metric="MAE", span=1, outdir=NULL){
  df<-data.frame(y, predicted_y)
  if(nrow(df)==length(SampleIDs)){
    df<-data.frame(SampleIDs, df)
  }else if(length(SampleIDs)>0 & nrow(df)!=length(SampleIDs)){
    stop("Please make sure that sample IDs match with y or predicted y.")
  }
  perf_value<-get.reg.performance(predicted_y, y)[[metric]]
  label = paste(metric, ": ", as.character(round(perf_value, 2)), sep="")
  p<-ggplot(df, aes(x=.data$y, y=.data$predicted_y))+
    ylab(paste("Predicted ",target_field,sep=""))+
    xlab(paste("Observed ",target_field,sep=""))+
    geom_point(alpha=0.1)+
    geom_smooth(method="loess",span=span)+
    annotate(geom="text", x=Inf, y=Inf, label=label, color="grey60", vjust=2, hjust=2)+
    theme_bw()
  #coord_flip()+
  if(!is.null(outdir)){
    ggsave(filename=paste(outdir, prefix, ".", target_field, ".obs_vs_pred.scatterplot.pdf",sep=""), plot=p, height=4, width=4)
    sink(paste(outdir, prefix, ".", target_field, ".obs_vs_pred.results.xls",sep=""));cat("\t")
    write.table(df, quote=FALSE,sep="\t", row.names = FALSE);sink(NULL)
  }
  p
}
#' @title plot_residuals
#' @description Plot the residuals of observed and predicted values from a ranger model.
#' @param y The numeric values for labeling data.
#' @param predicted_y The predicted values for y.
#' @param SampleIDs The sample ids in the data table.
#' @param prefix The prefix for the dataset in the training or testing.
#' @param target_field A string indicating the target field of the metadata for regression.
#' @param outdir The output directory.
#' @examples
#' set.seed(123)
#' x <- data.frame(rbind(t(rmultinom(7, 75, c(.201,.5,.02,.18,.099))),
#'             t(rmultinom(8, 75, c(.201,.4,.12,.18,.099))),
#'             t(rmultinom(15, 75, c(.011,.3,.22,.18,.289))),
#'             t(rmultinom(15, 75, c(.091,.2,.32,.18,.209))),
#'             t(rmultinom(15, 75, c(.001,.1,.42,.18,.299)))))
#' y<- 1:60
#' rf_model<-rf.out.of.bag(x, y)
#' plot_residuals(y, predicted_y=rf_model$predicted, prefix="train", target_field="age")
#' @author Shi Huang
#' @export
plot_residuals <- function(y, predicted_y, SampleIDs=NULL, prefix="train", target_field="value", outdir=NULL){
  df<-data.frame(y, predicted_y, rsdl=y-predicted_y)
  if(nrow(df)==length(SampleIDs)){
    df<-data.frame(SampleIDs, df)
  }else if(length(SampleIDs)>0 & nrow(df)!=length(SampleIDs)){
    stop("Please make sure that sample IDs match with y or predicted y.")
  }
  p<-ggplot(df, aes(x=.data$y, y=.data$rsdl))+
    ylab(paste("Residuals of prediceted ",target_field,sep=""))+
    xlab(paste("Observed ",target_field,sep=""))+
    geom_point(alpha=0.1)+
    geom_hline(yintercept=0)+
    #geom_smooth(method="loess",span=span)+
    theme_bw()
  if(!is.null(outdir)){
  ggsave(filename=paste(outdir, prefix, ".", target_field, ".obs_vs_residuals_of_pred.scatterplot.pdf",sep=""), plot=p, height=4, width=4)
  sink(paste(outdir, prefix, ".", target_field, ".obs_vs_pred.results.xls",sep=""))
  cat("\t")
  write.table(df, quote=FALSE,sep="\t", row.names = FALSE);sink(NULL)
  }
  invisible(p)
}

#' @title plot_perf_VS_rand
#' @description It output a histogram showing the performance of a real regression model and null models.
#' @param x The train data.
#' @param y The numeric labeling data.
#' @param nfolds The number of folds in the cross validation. If nfolds > length(y)
#' or nfolds==-1, uses leave-one-out cross-validation. If nfolds was a factor, it means customized folds (e.g., leave-one-group-out cv) were set for CV.
#' @param predicted_y The predicted values for y.
#' @param n_features The number of features in the training data.
#' @param prefix The prefix for the dataset in the training or testing.
#' @param target_field A string indicating the target field of the metadata for regression.
#' @param metric The regression performance metric applied, including MAE, RMSE, MSE, R_squared, Adj_R_squared, or Separman_rho.
#' @param permutation The permutation times for a random guess of regression performance.
#' @param outdir The output directory.
#' @examples
#' set.seed(123)
#' x <- data.frame(rbind(t(rmultinom(7, 75, c(.201,.5,.02,.18,.099))),
#'             t(rmultinom(8, 75, c(.201,.4,.12,.18,.099))),
#'             t(rmultinom(15, 75, c(.011,.3,.22,.18,.289))),
#'             t(rmultinom(15, 75, c(.091,.2,.32,.18,.209))),
#'             t(rmultinom(15, 75, c(.001,.1,.42,.18,.299)))))
#' y<- 1:60
#' rf_model<-rf.out.of.bag(x, y)
#' p<-plot_perf_VS_rand(y=y, predicted_y=rf_model$predicted, prefix="train",
#' permutation=100, metric="MAE", target_field="age", n_features=5)
#' p
#' @author Shi Huang
#' @export
plot_perf_VS_rand<-function(x, y, predicted_y, prefix="train", target_field, nfolds,
                            metric="MAE", permutation=100, n_features, outdir=NULL){
  # rand_perf <- function(y, permutation.=permutation, metric.=metric,
  #                       n_features.=n_features){
  #   set.seed(123)
  #   rand_y_mat <-replicate(permutation, sample(y, replace = FALSE))
  #   rand_perf_values <- apply(rand_y_mat, 2, function(x) get.reg.performance(x, y, n_features)[[metric]])
  #   rand_perf_values
  # }
  # rand_perf_values=rand_perf(y)
  shuffle_y_perf <- function(x, y, permutation.=permutation, metric.=metric, nfolds.=nfolds,
                        n_features.=n_features){
    set.seed(123)
    rand_y_mat <-replicate(permutation, sample(y, replace = FALSE))
    rand_perf_values <- apply(rand_y_mat, 2, function(rand_y){
      rand_rf <- rf.cross.validation(x, rand_y, nfolds=nfolds)
      perf <- get.reg.performance(rand_y, rand_rf$predicted, n_features)[[metric]]
    })
    rand_perf_values
  }
  rand_perf_values=shuffle_y_perf(x, y, permutation = 100)

  perf_value<-get.reg.performance(predicted_y, y, n_features)[[metric]]
  perf_values<-data.frame(perf_value, rand_perf_values)

  emp_p_value <- function(perf_value, rand_perf_values){
    k<-length(rand_perf_values)
    #p<-(sum(abs(perf_value >= rand_perf_values))+1)/(k+1)
    p<-(sum(abs(perf_value >= rand_perf_values)))/(k)
    p
  }
  emp_p_value<-emp_p_value(perf_value, rand_perf_values)
  label = paste(metric, ": ", as.character(round(perf_value, 2)), "\np-value = ", round(emp_p_value, 2), sep="")
  p<-ggplot(perf_values, aes(x=.data$rand_perf_values)) + geom_histogram(alpha=0.5) +
    xlab(metric)+
    ylab("count")+
    geom_vline(data=perf_values, aes(xintercept = .data$perf_value)) +
    annotate(geom="text", x=perf_value, y=Inf, label=label, color="red", vjust=2, hjust=0)+theme_bw()
  p
  if(!is.null(outdir)){
    ggsave(filename=paste(outdir, prefix, ".", target_field, ".", metric,
                          "_vs_rand.histogram.pdf",sep=""), plot=p, height=4, width=4)
  }
  p
}

#' @title plot_train_vs_test
#' @description Plot the residuals of observed and predicted values from a ranger model.
#' @param train_y The numeric labels for training data.
#' @param predicted_train_y The predicted values for training data.
#' @param train_SampleIDs The sample ids in the train data.
#' @param test_y The numeric labels for testing data.
#' @param predicted_test_y The predicted values for test data.
#' @param test_SampleIDs The sample ids in the test data.
#' @param train_prefix The prefix for the dataset in the training data.
#' @param test_prefix The prefix for the dataset in the testing data.
#' @param train_target_field A string indicating the target field of the training metadata for regression.
#' @param test_target_field A string indicating the target field of the testing metadata for regression.
#' @param outdir The output directory.
#' @examples
#' set.seed(123)
#' x <- data.frame(rbind(t(rmultinom(7, 75, c(.201,.5,.02,.18,.099))),
#'             t(rmultinom(8, 75, c(.201,.4,.12,.18,.099))),
#'             t(rmultinom(15, 75, c(.011,.3,.22,.18,.289))),
#'             t(rmultinom(15, 75, c(.091,.2,.32,.18,.209))),
#'             t(rmultinom(15, 75, c(.001,.1,.42,.18,.299)))))
#' y<- 1:60
#' rf_model<-rf.out.of.bag(x, y)
#' plot_perf_VS_rand(y=y, predicted_y=rf_model$predicted, prefix="train",
#' permutation=100, metric="MAE", target_field="age")
#' @author Shi Huang
#' @export
plot_train_vs_test<-function(train_y, predicted_train_y, test_y, predicted_test_y, train_SampleIDs=NULL, test_SampleIDs=NULL,
                             train_prefix="train", test_prefix="test", train_target_field, test_target_field, outdir=NULL){
  train.pred<-data.frame(value=train_y,predicted_value=predicted_train_y)
  if(nrow(train.pred)==length(train_SampleIDs)){
    train.pred<-data.frame(SampleIDs=train_SampleIDs, train.pred)
  }else if(length(train_SampleIDs)>0 & nrow(train.pred)!=length(train_SampleIDs)){
    stop("Please make sure that sample IDs match with y or predicted y in the train data.")
  }
  test.pred<-data.frame(value=test_y,predicted_value=predicted_test_y)
  if(nrow(test.pred)==length(test_SampleIDs)){
    test.pred<-data.frame(SampleIDs=test_SampleIDs, test.pred)
  }else if(length(test_SampleIDs)>0 & nrow(test.pred)!=length(test_SampleIDs)){
    stop("Please make sure that sample IDs match with y or predicted y in the test data.")
  }
  data_name<-c(rep(train_prefix,nrow(train.pred)),rep(test_prefix,nrow(test.pred)))
  pred<-data.frame(data=data_name,rbind(train.pred,test.pred))
  pred$data<-factor(pred$data,levels=c(train_prefix,test_prefix),ordered=TRUE)
  l<-levels(pred$data); l_sorted<-sort(levels(pred$data))
  #Mycolor <- c("#D55E00", "#0072B2")
  #if(identical(order(l), order(l_sorted))){Mycolor=Mycolor }else{Mycolor=rev(Mycolor)}
  p<-ggplot(pred,aes(x=.data$value,y=.data$predicted_value))+
    ylab(paste("Predicted ",train_target_field,sep=""))+
    xlab(paste("Observed ",train_target_field,sep=""))+
    geom_point(aes(color=.data$data), alpha=0.1)+
    geom_smooth(aes(color=.data$data), method="loess",span=1)+
    #scale_color_manual(values = Mycolor)+
    theme_bw() +
    facet_wrap(~.data$data)+
    theme(legend.position="none")
  #coord_flip()+
  if(!is.null(outdir)){
  ggsave(filename=paste(outdir, train_prefix,"-",test_prefix,".",test_target_field, ".train_test_ggplot.pdf",sep=""),plot=p, height=3, width=6)
  sink(paste(outdir, train_prefix,"-",test_prefix,".",test_target_field, ".train_test_results.xls",sep=""));
  cat("\t");write.table(pred, quote=FALSE,sep="\t", row.names = FALSE);sink(NULL)
  }
  invisible(p)
}

#' @title plot_reg_feature_selection
#' @description Plot the regression performance against the reduced number of features used in the modeling.
#' @param x The data frame or data matrix for model training.
#' @param y The numeric values for labeling data.
#' @param nfolds The number of folds in the cross-validation for each feature set.
#' @param unit The unit of numeric metadata variables that can be printed in the output figure.
#' @param rf_reg_model The rf regression model from \code{rf.out.of.bag}
#' @param metric The regression performance metric applied.
#' This must be one of "MAE", "RMSE", "MSE", "MAPE".
#' @param outdir The output directory.
#' @examples
#' set.seed(123)
#' x <- data.frame(rbind(t(rmultinom(7, 7*10, c(.201,.5,.02,.18,.099, .201,.5,.02,.18,.099))),
#'             t(rmultinom(8, 8*10, c(.201,.4,.12,.18,.099, .201,.5,.02,.18,.099))),
#'             t(rmultinom(15, 15*10, c(.011,.3,.22,.18,.289, .091,.2,.32,.18,.209))),
#'             t(rmultinom(15, 15*10, c(.091,.2,.32,.18,.209, .011,.3,.22,.18,.289))),
#'             t(rmultinom(15, 15*10, c(.001,.1,.42,.18,.299, .001,.1,.42,.18,.299)))))
#' y<- 1:60
#' rf_reg_model<-rf.out.of.bag(x, y)
#' rf_reg_model<-rf.cross.validation(x, y)
#' plot_reg_feature_selection(x, y, rf_reg_model, metric="MAE", outdir="./")
#' @author Shi Huang
#' @export
plot_reg_feature_selection <- function(x, y, rf_reg_model, nfolds=5, metric="MAE", unit=NA, outdir=NULL){
  if(class(rf_reg_model)=="rf.cross.validation"){
    rank_mat <- apply(rf_reg_model$importances, 2, function(x){rank(-x)})
    rf_imp_rank <- rank(apply(rank_mat, 1, median))
    }else if(class(rf_reg_model)=="rf.out.of.bag"){
      rf_imp_rank<-rank(-(rf_reg_model$importances))
    }else{
      stop("The class of input rf model should be rf.out.of.bag or rf.cross.validation.")
    }
  top_n_rf_list <-list()
  max_n<-max(rf_imp_rank, na.rm = TRUE)
  n_total_features<-ncol(x)
  n_features<-c(2, 4, 8, 16, 32, 64, 128, 256, 512, 1024)
  n_features<-n_features[1:which.min(abs(n_features-n_total_features))]
  top_n_perf<-matrix(NA, ncol=6, nrow=length(n_features)+1)
  colnames(top_n_perf)<-c("n_features", "MSE", "RMSE", "MAE", "MAPE", "Spearman_rho")
  rownames(top_n_perf)<-top_n_perf[,1]<-c(n_features, max_n)
  for(i in 1:length(n_features)){
    idx<-which(rf_imp_rank<=n_features[i])
    #top_n_features<-names(rf_imp_rank[idx])
    x_n<-x[, idx]
    y_n<-y
    top_n_rf_list[[i]] <- top_n_rf<-rf.cross.validation(x_n, y_n, nfolds=5, ntree=500)
    top_n_rf_perf<-get.reg.performance(top_n_rf$predicted, y_n)
    top_n_perf[i, 1]<-n_features[i]
    top_n_perf[i, 2]<-top_n_rf_perf$MSE
    top_n_perf[i, 3]<-top_n_rf_perf$RMSE
    top_n_perf[i, 4]<-top_n_rf_perf$MAE
    top_n_perf[i, 5]<-top_n_rf_perf$MAPE
    top_n_perf[i, 6]<-top_n_rf_perf$Spearman_rho
  }

  all_rf_perf<-get.reg.performance(rf_reg_model$predicted, y)
  top_n_perf[length(n_features)+1, ]<-as.numeric(c(max_n, all_rf_perf[c("MSE", "RMSE", "MAE", "MAPE", "Spearman_rho")]))
  top_n_perf<-data.frame(top_n_perf)
  breaks<-top_n_perf$n_features
  y_label <- ifelse(is.na(unit), metric, paste(metric, " (", unit ,")", sep=""))
  p<-ggplot(top_n_perf, aes(x=.data$n_features, y=get(metric))) +
    xlab("# of features used")+
    ylab(y_label)+
    scale_x_continuous(trans = "log", breaks=breaks)+
    geom_point() + geom_line()+
    theme_bw()+
    theme(axis.line = element_line(color="black"),
          axis.title = element_text(size=18),
          strip.background = element_rect(colour = "white"),
          panel.border = element_blank())
  if(!is.null(outdir)){
  ggsave(filename=paste(outdir,"rf__",metric,"__top_rankings.scatterplot.pdf",sep=""), plot=p, width=5, height=4)
  }
  res <- list()
  res$top_n_perf <- top_n_perf
  res$top_n_rf <-top_n_rf_list
  res
}

#' @title calc_rel_predicted
#' @description Calculate the relative predicted values to the spline fit.
#' @param train_y The numeric labels for training data.
#' @param predicted_train_y The predicted values for training data.
#' @param train_SampleIDs The sample ids in the train data.
#' @param test_y The numeric labels for testing data.
#' @param predicted_test_y The predicted values for test data.
#' @param test_SampleIDs The sample ids in the test data.
#' @param train_prefix The prefix for the dataset in the training data.
#' @param test_prefix The prefix for the dataset in the testing data.
#' @param train_target_field A string indicating the target field of the training metadata for regression.
#' @param test_target_field A string indicating the target field of the testing metadata for regression.
#' @param outdir The output directory.
#' @examples
#' set.seed(123)
#' train_x <- data.frame(rbind(t(rmultinom(7, 75, c(.201,.5,.02,.18,.099))),
#'             t(rmultinom(8, 75, c(.201,.4,.12,.18,.099))),
#'             t(rmultinom(15, 75, c(.011,.3,.22,.18,.289))),
#'             t(rmultinom(15, 75, c(.091,.2,.32,.18,.209))),
#'             t(rmultinom(15, 75, c(.001,.1,.42,.18,.299)))))
#' train_y<- 1:60
#' test_x <- data.frame(rbind(t(rmultinom(7, 75, c(.201,.5,.02,.18,.099))),
#'             t(rmultinom(8, 75, c(.201,.4,.12,.18,.099))),
#'             t(rmultinom(15, 75, c(.011,.3,.22,.18,.289))),
#'             t(rmultinom(15, 75, c(.091,.2,.32,.18,.209)))))
#' test_y<- 1:45
#' train_rf_model<-rf.out.of.bag(train_x, train_y)
#' predicted_test_y<-predict(train_rf_model$rf.model, test_x)$predictions
#' calc_rel_predicted(train_y, predicted_train_y=train_rf_model$predicted)
#' calc_rel_predicted(train_y=train_y, predicted_train_y=train_rf_model$predicted,
#'                    #train_SampleIDs=as.character(1:60),
#'                    test_y=test_y, predicted_test_y=predicted_test_y,
#'                    #test_SampleIDs=as.character(1:45),
#'                    train_target_field="y",  test_target_field="test_y", outdir=NULL)
#' calc_rel_predicted(train_y=train_y, predicted_train_y=train_rf_model$predicted,
#'                    train_SampleIDs=as.character(1:60),
#'                    test_y=test_y, predicted_test_y=predicted_test_y,
#'                    test_SampleIDs=as.character(1:45),
#'                    train_target_field="y",  test_target_field="test_y", outdir=NULL)
#' @author Shi Huang
#' @export
calc_rel_predicted<-function(train_y, predicted_train_y, train_SampleIDs=NULL,
                             test_y=NULL, predicted_test_y=NULL, test_SampleIDs=NULL,
                             train_prefix="train", test_prefix="test",
                             train_target_field="y", test_target_field="y", outdir=NULL){
  spl_train <- smooth.spline(train_y, predicted_train_y)
  train_relTrain <- residuals(spl_train)
  train_fittedTrain <- fitted(spl_train)
  relTrain_data<-train_relTrain_data <- data.frame(y=train_y, predicted_y=predicted_train_y,
                                                   rel_predicted_y=train_relTrain,
                                                   fitted_predicted_y=train_fittedTrain)
  if(!is.null(train_SampleIDs)){
    if(nrow(relTrain_data)!=length(train_SampleIDs)){
      stop("Please make sure that sample IDs match with y or predicted y in the train data.")
    }else{
      relTrain_data<-train_relTrain_data<-data.frame(SampleIDs=train_SampleIDs, relTrain_data)
    }
  }
  sink(paste(outdir, train_prefix,".Relative_",train_target_field,".results.xls",sep=""));cat("\t");
  write.table(relTrain_data,quote=FALSE,sep="\t", row.names = FALSE);sink(NULL)

  if(!is.null(test_y) & !is.null(predicted_test_y)){
    test_fitted<-predict(spl_train, test_y)$y
    test_relTrain <- predicted_test_y - test_fitted
    test_relTrain_data <- data.frame(y=test_y, predicted_y=predicted_test_y,
                                     rel_predicted_y=test_relTrain,
                                     fitted_predicted_y=test_fitted)
    if(!is.null(test_SampleIDs)){
      if(nrow(test_relTrain_data)!=length(test_SampleIDs)){
        stop("Please make sure that sample IDs match with y or predicted y in the test data.")
        }else{
          test_relTrain_data<-data.frame(SampleIDs=test_SampleIDs, test_relTrain_data)
        }
    }
    relTrain_data<-rbind(train_relTrain_data, test_relTrain_data)
    DataSet <- factor(c(rep(train_prefix,length(train_relTrain)), rep(test_prefix,length(test_relTrain)) ))
    relTrain_data<-data.frame(relTrain_data, DataSet)
    sink(paste(outdir, train_prefix,"-",test_prefix,".Relative_",train_target_field,".results.xls",sep=""));
    cat("\t");
    write.table(relTrain_data,quote=FALSE,sep="\t", row.names = FALSE);sink(NULL)
  }
    return(relTrain_data)
}

#' @title plot_rel_predicted
#' @description Calculate the relative predicted values to the spline fit in the training data.
#' @param relTrain_data The output dataframe of \code{calc_rel_predicted}.
#' @param prefix The prefix of a train dataset.
#' @param target_field A string indicating the target field in the metadata for regression.
#' @param outdir The output directory.
#' @examples
#' set.seed(123)
#' train_x <- data.frame(rbind(t(rmultinom(7, 75, c(.201,.5,.02,.18,.099))),
#'             t(rmultinom(8, 75, c(.201,.4,.12,.18,.099))),
#'             t(rmultinom(15, 75, c(.011,.3,.22,.18,.289))),
#'             t(rmultinom(15, 75, c(.091,.2,.32,.18,.209))),
#'             t(rmultinom(15, 75, c(.001,.1,.42,.18,.299)))))
#' train_y<- 1:60
#' test_x <- data.frame(rbind(t(rmultinom(7, 75, c(.201,.5,.02,.18,.099))),
#'             t(rmultinom(8, 75, c(.201,.4,.12,.18,.099))),
#'             t(rmultinom(15, 75, c(.011,.3,.22,.18,.289))),
#'             t(rmultinom(15, 75, c(.091,.2,.32,.18,.209)))))
#' test_y<- 1:45
#' train_rf_model<-rf.out.of.bag(train_x, train_y)
#' predicted_test_y<-predict(train_rf_model$rf.model, test_x)$predictions
#' relTrain_data<-calc_rel_predicted(train_y, train_rf_model$predicted, test_y, predicted_test_y)
#' plot_rel_predicted(relTrain_data)
#' @author Shi Huang
#' @export
plot_rel_predicted <- function(relTrain_data, prefix="train", target_field="value", outdir=NULL){
  if(length(prefix)==1){
    relTrain_data<-subset(relTrain_data, .data$DataSet==prefix)
    p<-ggplot(relTrain_data, aes(x=.data$y, y=.data$rel_predicted_y))+
      ylab(paste("Relative prediceted ",target_field,sep=""))+
      xlab(paste("Observed ",target_field,sep=""))+
      geom_point(alpha=0.1)+
      geom_hline(yintercept=0)+
      #geom_smooth(method="loess",span=span)+
      theme_bw()
    #coord_flip()+
  }else{
    p<-ggplot(relTrain_data, aes(x=.data$y, y=.data$rel_predicted_y, color=.data$DataSet))+
      ylab(paste("Relative prediceted ",target_field,sep=""))+
      xlab(paste("Observed ",target_field,sep=""))+
      geom_point(alpha=0.1)+
      geom_hline(yintercept=0)+
      #geom_smooth(method="loess",span=span)+
      theme_bw()
    prefix=paste(prefix, collapse = "-")
  }
  ggsave(filename=paste(outdir, prefix, ".", target_field, ".obs_vs_relative_pred.scatterplot.pdf",sep=""), plot=p, height=4, width=4)
  invisible(p)
}

#' @title boxplot_rel_predicted_train_vs_test
#' @description Make the boxplot the relative predicted values in both train and test datasets.
#' @param relTrain_data The output dataframe of \code{calc_rel_predicted}.
#' @param train_prefix The prefix for the dataset in the training data.
#' @param test_prefix The prefix for the dataset in the testing data.
#' @param train_target_field A string indicating the target field of the training metadata for regression.
#' @param outdir The output directory.
#' @examples
#' set.seed(123)
#' train_x <- data.frame(rbind(t(rmultinom(7, 75, c(.201,.5,.02,.18,.099))),
#'             t(rmultinom(8, 75, c(.201,.4,.12,.18,.099))),
#'             t(rmultinom(15, 75, c(.011,.3,.22,.18,.289))),
#'             t(rmultinom(15, 75, c(.091,.2,.32,.18,.209))),
#'             t(rmultinom(15, 75, c(.001,.1,.42,.18,.299)))))
#' train_y<- 1:60
#' test_x <- data.frame(rbind(t(rmultinom(7, 75, c(.201,.5,.02,.18,.099))),
#'             t(rmultinom(8, 75, c(.201,.4,.12,.18,.099))),
#'             t(rmultinom(15, 75, c(.011,.3,.22,.18,.289))),
#'             t(rmultinom(15, 75, c(.091,.2,.32,.18,.209)))))
#' test_y<- 1:45
#' train_rf_model<-rf.out.of.bag(train_x, train_y)
#' predicted_test_y<-predict(train_rf_model$rf.model, test_x)$predictions
#' relTrain_data<-calc_rel_predicted(train_y, predicted_train_y=train_rf_model$predicted,
#'                                   train_SampleIDs=NULL,
#'                                   test_y, predicted_test_y,
#'                                   test_SampleIDs=NULL)
#' boxplot_rel_predicted_train_vs_test(relTrain_data)
#' @author Shi Huang
#' @export
boxplot_rel_predicted_train_vs_test<-function(relTrain_data, train_target_field="value",
                                              train_prefix="train", test_prefix="test", outdir=NULL){
    NoTestDataset<-all(grepl("DataSet", colnames(relTrain_data))==FALSE)
    if(NoTestDataset) stop("Test dataset should be included for residuals comparison between train and test datasets!")
    p_w<-formatC(stats::wilcox.test(rel_predicted_y~DataSet, data = relTrain_data)$p.value,digits=4,format="g")
    # l<-levels(relTrain_data$DataSet); l_sorted<-sort(levels(relTrain_data$DataSet))
    # Mycolor <- rep(c("#D55E00", "#0072B2"), length.out=length(l))
    # if(identical(order(l), order(l_sorted))){
    #   Mycolor=Mycolor; l_ordered=l
    # }else{Mycolor=rev(Mycolor); l_ordered=l_sorted}
    p<-ggplot(relTrain_data, aes(x=.data$DataSet, y=.data$rel_predicted_y)) +
      geom_violin(aes(color=.data$DataSet))+
      geom_boxplot(outlier.shape = NA, width=0.4)+
      geom_jitter(position=position_jitter(width=0.2),alpha=0.1) + # aes(color=DataSet),
      geom_hline(yintercept=0)+
      theme_bw()+
      ggtitle(paste("Wilcoxon Rank Sum Test:\n P=", p_w, sep=""))+
      xlab("Data sets") +
      ylab(paste("Relative predicted", train_target_field))+
      theme(legend.position="none")
    # p<-p+scale_color_manual(values = Mycolor, labels=l_ordered)
    ggsave(filename= paste(outdir, train_prefix,"-",test_prefix, ".Relative_",train_target_field,".boxplot.pdf",sep=""), width=3, height=4)
    invisible(p)
}

