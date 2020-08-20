#' @importFrom doMC registerDoMC
#' @importFrom foreach foreach %dopar%
#' @importFrom parallel detectCores
#' @importFrom caret confusionMatrix


#' @title rf_clf.pairwise
#' @description Perform pairwise rf classfication for a data matrix between all pairs of group levels
#' @param df A data matrix or data.frame
#' @param f A factor with more than two levels.
#' @param nfolds The number of folds in the cross validation.
#' @param ntree The number of trees.
#' @param verbose The boolean value indicating if the computation status and estimated runtime shown.
#' @return A summary table containing the performance of Random Forest classification models
#' @seealso ranger
#' @examples
#' df <- t(rmultinom(16,160,c(.001,.6,.2,.3,.299))) + 0.65
#' f<-factor(c(rep("A", 4), rep("B", 4), rep("C", 4), rep("D", 4)))
#' f0<-factor(c(rep("A", 8), rep("B", 8)))
#' rf_clf.pairwise(df, f)
#' @author Shi Huang
#' @export
rf_clf.pairwise <- function (df, f, nfolds=3, ntree=5000, verbose=FALSE) {
  ## TODO: check the class of inputs
  rf_compare_levels <- function(df, f, i=1, j=2, nfolds=3, ntree=500, verbose=FALSE) {
    df_ij <- df[which(as.integer(f) == i | as.integer(f) == j), ]
    f_ij <- factor(f[which(as.integer(f) == i | as.integer(f) == j)])
    if(nfolds==3){
      oob<-rf.out.of.bag(df_ij, f_ij, verbose=verbose, ntree=ntree)
    }else{
      oob<-rf.cross.validation(df_ij, f_ij, nfolds=nfolds, verbose=verbose, ntree=ntree)
    }
    positive_class=levels(f)[i]
    cat("\nTraining dataset: ", levels(f)[i], "-", levels(f)[j] ,"\n\n")
    conf<-caret::confusionMatrix(data=oob$predicted, f_ij, positive=positive_class)
    acc<-conf$overall[1]
    kappa_oob<-conf$overall[2]
    cat("Accuracy in the cross-validation: ", acc ,"\n")
    if(nlevels(f_ij)==2){rf_AUC<-get.auroc(oob$probabilities, f_ij, positive_class)}else{rf_AUC<-NA}
    cat("AUC in the cross-validation: ", rf_AUC ,"\n") # wired value of "1" kept showing
    c("AUC"=rf_AUC, acc, kappa_oob, conf$byClass)
  }
  if(nlevels(f)==2){
    out_summ<-rf_compare_levels(df, f, i=1, j=2, nfolds, ntree, verbose)
  }else{
    level_names<-levels(f)
    ix <- stats::setNames(seq_along(level_names), level_names)
    out_list<-outer(ix[-1L], ix[-length(ix)],
                    function(ivec, jvec) sapply(seq_along(ivec),
                                                function(k) {
                                                  i <- ivec[k]
                                                  j <- jvec[k]
                                                  if (i > j) {
                                                    rf_compare_levels(df, f, i, j, nfolds, ntree, verbose)
                                                  }else{NA}
                                                }))
    dataset_name<-outer(dimnames(out_list)[[1]],dimnames(out_list)[[2]], paste, sep="__")
    out_rownames<-dataset_name[lower.tri(dataset_name, diag = T)]
    out_summ<-do.call(rbind, out_list[lower.tri(out_list, diag = T)]); rownames(out_summ)<-out_rownames
  }
  class(out_summ)<-"rf_clf.pairwise"
  out_summ
}

#' @title rf_clf.by_datasets
#' @description It runs standard random forests with oob estimation for classification of
#' c_category in each the sub-datasets splited by the s_category,
#' and apply the model to all the other datasets. The output includes
#' accuracy, auc and Kappa statistics.
#' @param df Training data: a data.frame.
#' @param metadata Sample metadata with at least two columns.
#' @param s_category A string indicates the category in the sample metadata: a ‘factor’ defines the sample grouping for data spliting.
#' @param c_category A indicates the category in the sample metadata as a responsive vector: if a 'factor', rf classification is performed in each of splited datasets.
#' @param positive_class A string indicates one class in the 'c_category' column of metadata.
#' @param nfolds The number of folds in the cross validation.
#' @param clr_transform A boolean value indicating if the clr-transformation applied.
#' @param rf_imp_pvalues A boolean value indicating if compute both importance score and pvalue for each feature.
#' @param verbose A boolean value indicating if show computation status and estimated runtime.
#' @param ntree The number of trees.
#' @param p.adj.method The p-value correction method, default is "bonferroni".
#' @param q_cutoff The cutoff of q values for features, the default value is 0.05.
#' @return ...
#' @seealso ranger
#' @examples
#' df <- data.frame(rbind(t(rmultinom(7, 75, c(.21,.6,.12,.38,.099))),
#'             t(rmultinom(8, 75, c(.001,.6,.42,.58,.299))),
#'             t(rmultinom(15, 75, c(.011,.6,.22,.28,.289))),
#'             t(rmultinom(15, 75, c(.091,.6,.32,.18,.209))),
#'             t(rmultinom(15, 75, c(.001,.6,.42,.58,.299)))))
#' df0 <- data.frame(t(rmultinom(60, 300,c(.001,.6,.2,.3,.299))))
#' metadata<-data.frame(f_s=factor(c(rep("A", 15), rep("B", 15), rep("C", 15), rep("D", 15))),
#'                      f_c=factor(c(rep("C", 7), rep("H", 8), rep("C", 7), rep("H", 8),
#'                                   rep("C", 7), rep("H", 8), rep("C", 7), rep("H", 8))),
#'                      f_d=factor(rep(c(rep("a", 5), rep("b", 5), rep("c", 5)), 4)))
#' system.time(rf_clf.by_datasets(df, metadata, s_category='f_s',
#'             c_category='f_c', positive_class="C"))
#' rf_clf.by_datasets(df, metadata, s_category='f_s', c_category='f_c',
#'                    positive_class="C", rf_imp_pvalues=TRUE)
#' rf_clf.by_datasets(df, metadata, s_category='f_s', c_category='f_d')
#' rf_clf.by_datasets(df, metadata, s_category='f_s', c_category='f_d', rf_imp_pvalues=TRUE)
#' @author Shi Huang
#' @export
rf_clf.by_datasets<-function(df, metadata, s_category, c_category, positive_class=NA,
                             rf_imp_pvalues=FALSE, clr_transform=TRUE, nfolds=3, verbose=FALSE, ntree=500,
                             p.adj.method = "BH", q_cutoff=0.05){
  ## TODO: check the class of inputs
  y_list<-split(metadata[, c_category], metadata[, s_category])
  x_list<-split(df, metadata[, s_category])
  datasets<-levels(factor(metadata[, s_category]))
  L<-length(y_list)
  positive_class<-ifelse(is.na(positive_class), levels(factor(y_list[[1]]))[1], positive_class)
  # 1. sample size of all datasets
  sample_size<-as.numeric(table(metadata[, s_category]))
  nCores <- parallel::detectCores()
  doMC::registerDoMC(nCores)
  # comb function for parallelization using foreach
  comb <- function(x, ...) {
    lapply(seq_along(x),
           function(i) c(x[[i]], lapply(list(...), function(y) y[[i]])))
  }
  oper<-foreach::foreach(i=1:L, .combine='comb', .multicombine=TRUE, .init=list(list(), list(), list())) %dopar% {
    x<-x_list[[i]]
    y<-factor(y_list[[i]])
    # 2. AUC of random forest model
    if(nfolds==3){
      oob<-rf.out.of.bag(x, y, verbose=verbose, ntree=ntree, imp_pvalues = rf_imp_pvalues)
      rf_imps=oob$importances
    }else{
      oob<-rf.cross.validation(x, y, nfolds=nfolds, verbose=verbose, ntree=ntree, imp_pvalues = rf_imp_pvalues)
      rf_imps=rowMeans(oob$importances)
    }
    if(nlevels(y)==2){rf_AUC<-get.auroc(oob$probabilities, y, positive_class) }else{ rf_AUC <- NA }
    # 3. # of significantly differential abundant features between health and disease
    out<-BetweenGroup.test(x, y, clr_transform=clr_transform, positive_class=positive_class, p.adj.method = p.adj.method, q_cutoff=q_cutoff)
    feature_imps<-data.frame(feature=rownames(out), dataset=rep(datasets[i], ncol(x)),
                             rf_imps=rf_imps, out)
    list(oob=oob, rf_AUC=rf_AUC, feature_imps=feature_imps)
  }

  result<-list()
  result$x_list<-x_list
  result$y_list<-y_list
  result$datasets<-datasets
  result$sample_size<-sample_size
  result$rf_model_list<-oper[[1]]
  result$rf_AUC<-unlist(oper[[2]])
  result$feature_imps_list<-oper[[3]]
  class(result)<-"rf_clf.by_datasets"
  return(result)
}

#' @title rf_reg.by_datasets
#' @description It runs standard random forests with oob estimation for regression of
#' \code{c_category} in each the sub-datasets splited by the \code{s_category},
#' and apply the model to all the other datasets. The output includes
#' accuracy, auc and Kappa statistics.
#' @param df Training data: a data.frame.
#' @param metadata Sample metadata with at least two columns.
#' @param s_category A string indicates the category in the sample metadata: a ‘factor’ defines the sample grouping for data spliting.
#' @param c_category A indicates the category in the sample metadata: a 'factor' used as sample label for rf classification in each of splited datasets.
#' @param rf_imp_pvalues A boolean value indicate if compute both importance score and pvalue for each feature.
#' @param ntree The number of trees.
#' @param nfolds The number of folds in the cross validation.
#' @param verbose Show computation status and estimated runtime.
#' @return ...
#' @seealso ranger
#' @examples
#' df <- data.frame(rbind(t(rmultinom(7, 75, c(.21,.6,.12,.38,.099))),
#'             t(rmultinom(8, 75, c(.001,.6,.42,.58,.299))),
#'             t(rmultinom(15, 75, c(.011,.6,.22,.28,.289))),
#'             t(rmultinom(15, 75, c(.091,.6,.32,.18,.209))),
#'             t(rmultinom(15, 75, c(.001,.6,.42,.58,.299)))))
#' df0 <- data.frame(t(rmultinom(60, 300,c(.001,.6,.2,.3,.299))))
#' metadata<-data.frame(f_s=factor(c(rep("A", 30), rep("B", 30))),
#'                      f_s1=factor(c(rep(TRUE, 30), rep(FALSE, 30))),
#'                      f_c=factor(c(rep("C", 15), rep("H", 15), rep("D", 15), rep("P", 15))),
#'                      age=c(1:30, 2:31)
#'                      )
#'
#' reg_res<-rf_reg.by_datasets(df, metadata, nfolds=5, s_category='f_s', c_category='age')
#' reg_res
#' @author Shi Huang
#' @export
rf_reg.by_datasets<-function(df, metadata, s_category, c_category, nfolds=3,
                             rf_imp_pvalues=FALSE, verbose=FALSE, ntree=500){
  ## TODO: check the class of inputs
  #as.numeric.factor <- function(y) {as.numeric(levels(y))[y]}
  y<-metadata[, c_category]
  if(is.factor(metadata[, c_category])) y<-as.numeric(as.character(y))
  y_list<-split(y, metadata[, s_category])
  x_list<-split(df, metadata[, s_category])
  # data check
  try(if(!all(unlist(lapply(y_list, length)) == unlist(lapply(x_list, nrow))))
    stop("# of samples in x and length of y should match to each other within all datasets!") )
  datasets<-levels(factor(metadata[, s_category]))
  L<-length(y_list)
  # sample size of all datasets
  sample_size<-as.numeric(table(metadata[, s_category]))
  nCores <- parallel::detectCores()
  doMC::registerDoMC(nCores)
  # comb function for parallelization using foreach
  comb <- function(x, ...) {
    lapply(seq_along(x),
           function(i) c(x[[i]], lapply(list(...), function(y) y[[i]])))
  }
  if(nfolds==3){
    oper_len=14
    out_list<-list(list())[rep(1, oper_len)]
    oper_names<-c("rf.model", "y", "predicted", "MSE", "RMSE", "nRMSE",
                  "MAE", "MAPE", "MASE", "R_squared", "Adj_R_squared",
                  "importances", "params", "error.type")
  }else if(nfolds!=3){
                  oper_len=15
                  out_list<-list(list())[rep(1, oper_len)]
                  oper_names<-c("rf.model", "y", "predicted", "MSE", "RMSE", "nRMSE",
                                "MAE", "MAPE", "MASE", "R_squared", "Adj_R_squared",
                                "importances", "params", "error.type", "nfolds")}
  oper<-foreach::foreach(i=1:L, .combine='comb', .multicombine=TRUE,
                .init=out_list) %dopar% {
                             if(nfolds==3){
                               oob<-rf.out.of.bag(x_list[[i]], y_list[[i]],
                                                  verbose=verbose, ntree=ntree,
                                                  imp_pvalues = rf_imp_pvalues)
                             }else{
                               oob<-rf.cross.validation(x_list[[i]], y_list[[i]], nfolds=nfolds,
                                                        verbose=verbose, ntree=ntree,
                                                        imp_pvalues = rf_imp_pvalues)
                             }
                  oob
                }
  names(oper)<-oper_names
  result<-list()
  result$x_list<-x_list
  result$y_list<-y_list
  result$datasets<-datasets
  result$sample_size<-sample_size
  result$rf_model_list<-oper$rf.model; names(result$rf_model_list)<-result$datasets
  result$rf_predicted<-oper$predicted; names(result$rf_predicted)<-result$datasets
  result$feature_imps_list<-oper$importances
  if(nfolds!=3) result$feature_imps_list<-lapply(result$feature_imps_list, rowMeans)
  names(result$feature_imps_list)<-result$datasets
  result$rf_MSE<-unlist(oper$MSE)
  result$rf_RMSE<-unlist(oper$RMSE)
  result$rf_nRMSE<-unlist(oper$nRMSE)
  result$rf_MAE<-unlist(oper$MAE)
  result$rf_MAPE<-unlist(oper$MAPE)
  result$rf_MASE<-unlist(oper$MASE)
  result$rf_R_squared<-unlist(oper$R_squared)
  result$rf_Adj_R_squared<-unlist(oper$Adj_R_squared)
  class(result)<-"rf_reg.by_datasets"
  return(result)
}

#' @title rf_clf.cross_appl
#' @description Based on pre-computed rf models classifying 'c_category' in each the sub-datasets splited by the 's_category',
#' perform cross-datasets application of the rf models. The inputs are precalculated
#' rf models, and the outputs include accuracy, auc and Kappa statistics.
#' @param rf_model_list A list of rf.model objects from \code{rf.out.of.bag}.
#' @param x_list A list of training datasets usually in the format of data.frame.
#' @param y_list A list of responsive vector for regression in the training datasets.
#' @param positive_class A string indicates one common class in each of elements in the y_list.
#' @return A object of class rf_clf.cross_appl including a list of performance summary and predicted values of all predictions
#' @seealso ranger
#' @examples
#' df <- data.frame(rbind(t(rmultinom(7, 75, c(.21,.6,.12,.38,.099))),
#'             t(rmultinom(8, 75, c(.001,.6,.42,.58,.299))),
#'             t(rmultinom(15, 75, c(.011,.6,.22,.28,.289))),
#'             t(rmultinom(15, 75, c(.091,.6,.32,.18,.209))),
#'             t(rmultinom(15, 75, c(.001,.6,.42,.58,.299)))))
#' df0 <- data.frame(t(rmultinom(60, 300,c(.001,.6,.2,.3,.299))))
#' metadata<-data.frame(f_s=factor(c(rep("A", 15), rep("B", 15), rep("C", 15), rep("D", 15))),
#'                      f_s0=factor(c(rep("A", 30), rep("B", 30))),
#'                      f_c=factor(c(rep("C", 7), rep("H", 8), rep("C", 7), rep("H", 8),
#'                                   rep("C", 7), rep("H", 8), rep("C", 7), rep("H", 8))),
#'                      age=c(1:30, 2:31)
#'                      )
#' res_list<-rf_clf.by_datasets(df, metadata, s_category='f_s', nfolds=5,
#'                              c_category='f_c', positive_class="C")
#' rf_model_list<-res_list$rf_model_list
#'
#' rf_clf.cross_appl(rf_model_list, res_list$x_list, res_list$y_list, positive_class="C")
#' #--------------------
#' comp_group="A"
#' comps_res<-rf_clf.comps(df, f=metadata[, 'f_s'], comp_group, verbose=FALSE,
#'                         ntree=500, p.adj.method = "bonferroni", q_cutoff=0.05)
#' comps_res
#' rf_clf.cross_appl(comps_res$rf_model_list,
#'                   x_list=comps_res$x_list,
#'                   y_list=comps_res$y_list,
#'                   positive_class=comp_group)
#' @author Shi Huang
#' @export
rf_clf.cross_appl<-function(rf_model_list, x_list, y_list, positive_class=NA){
  ## TODO: check the class of inputs
  L<-length(rf_model_list)
  y_list<-lapply(y_list,factor)
  positive_class<-ifelse(is.na(positive_class), levels(factor(y_list[[1]]))[1], positive_class)
  try(if(!identical(L, length(x_list), length(y_list))) stop("The length of x list, y list and rf model list should be identical."))
  perf_summ<-data.frame(matrix(NA, ncol=17, nrow=L*L))
  colnames(perf_summ)<-c("Train_data", "Test_data", "Validation_type", "Accuracy", "AUC", "Kappa",
                         "Sensitivity", "Specificity", "Pos_Pred_Value","Neg_Pred_Value", "Precision", "Recall",
                         "F1", "Prevalence", "Detection_Rate", "Detection_Prevalence", "Balanced_Accuracy")
  predicted<-matrix(list(), ncol=2, nrow=L*L)
  colnames(predicted)<-c("predicted", "probabilities")
  for(i in 1:L){
    y<-y_list[[i]]
    x<-x_list[[i]]
    try(if(nlevels(y)==1) stop("Less than one level in the subgroup for classification"))
    #rf_model<-randomForest(x, y, ntree=5000, importance=T)
    #oob<-rf.out.of.bag(x, y, nfolds=nfolds, verbose=verbose, ntree=ntree)
    oob<-rf_model_list[[i]]
    #---
    #  RF Training accuracy
    #---
    cat("\nTraining dataset: ", names(x_list)[i] ,"\n\n")
    conf<-confusionMatrix(data=oob$predicted, oob$y, positive=positive_class)
    acc<-conf$overall[1]
    kappa_oob<-conf$overall[2]
    cat("Accuracy in the self-validation: ", acc ,"\n")
    #---
    #  AUC computation using "pROC" package
    #---
    auc<-get.auroc(oob$probabilities, oob$y, positive_class)
    cat("AUC in the self-validation: ", auc ,"\n")

    D<-1:L
    T<-D[D!=i]
    a=1+(i-1)*L
    perf_summ[a, 1:3]<-c(names(x_list)[i], names(x_list)[i], "self_validation")
    perf_summ[a, 4:17]<-c(acc, auc, kappa_oob, conf$byClass)
    predicted[a, 1][[1]]<-data.frame(test_y=y, pred_y=oob$predicted)
    predicted[a, 2][[1]]<-oob$probabilities
    #rownames(predicted)[a]<-paste(names(x_list)[i], names(x_list)[i], sep="__VS__")
    loop_num<-1
    for(j in T){
      if(nrow(x_list[[j]])>0){
        newx<-x_list[[j]]
        newy<-y_list[[j]]
        if(class(oob)=="rf.out.of.bag"){
           oob_rf.model<-oob$rf.model
          }else{
            oob_rf.model<-oob$rf.model[[1]]} # pick one sout of K rf models in the "rf.cross.validation"
        #pred_prob<-predict(oob$rf.model, x_list[[j]], type="prob") # for regular rf.out.of.bag (randomForest) function
        pred_prob <- get.predict.probability.from.forest(oob_rf.model, newx) # ranger only
        pred_prob<-pred_prob[,order(colnames(pred_prob))] # to avoid unanticipated order of numeric levels of factor y
        pred_newy<-factor(predict(oob_rf.model, newx, type="response")$predictions) # ranger only
        if(identical(levels(newy), levels(oob$y))){
          colnames(pred_prob)<- levels(newy)
          levels(pred_newy)<- levels(newy)
          #---Accuracy
          cat("Test dataset: ", names(x_list)[j] ,"\n")
          test_conf<-confusionMatrix(data=pred_newy, newy, positive=positive_class)
          test_acc<-test_conf$overall[1]
          test_kappa<-test_conf$overall[2]
          cat("Accuracy in the cross-applications: ", test_acc ,"\n")
          #---AUC
          test_auc<-get.auroc(pred_prob, newy, positive_class)
          cat("AUC in the cross-applications: ", test_auc ,"\n")
          perf_summ[a+loop_num, 4:17]<-c(test_acc, test_auc, test_kappa, test_conf$byClass)
        }else{
          colnames(pred_prob)<- levels(oob$y)
          levels(pred_newy)<- levels(oob$y)
          perf_summ[a+loop_num, 4:17]<-rep(NA, 14)
        }
        perf_summ[a+loop_num, 1:3]<-c(names(x_list)[i], names(x_list)[j], "cross_application")
        predicted[a+loop_num, 1][[1]]<-data.frame(test_y=newy, pred_y=pred_newy)
        predicted[a+loop_num, 2][[1]]<-pred_prob
        #rownames(predicted)[a+loop_num]<-paste(names(x_list)[i], names(x_list)[j], sep="__VS__")
        loop_num<-loop_num+1
      }
    }
  }
  res<-list()
  res$perf_summ<-perf_summ
  res$predicted<-predicted
  class(res)<-"rf_clf.cross_appl"
  res
}

#' @title rf_reg.cross_appl
#' @description Based on pre-computed rf models regressing \code{c_category} in each the sub-datasets splited by the \code{s_category},
#' perform cross-datasets application of the rf models. The inputs are precalculated
#' rf regression models, x_list and y_list.
#' @param rf_list A list of rf.model objects from \code{rf.out.of.bag}.
#' @param x_list A list of training datasets usually in the format of data.frame.
#' @param y_list A list of responsive vector for regression in the training datasets.
#' @return ...
#'
#' @seealso ranger
#' @examples
#' df <- data.frame(rbind(t(rmultinom(7, 75, c(.21,.6,.12,.38,.099))),
#'             t(rmultinom(8, 75, c(.001,.6,.42,.58,.299))),
#'             t(rmultinom(15, 75, c(.011,.6,.22,.28,.289))),
#'             t(rmultinom(15, 75, c(.091,.6,.32,.18,.209))),
#'             t(rmultinom(15, 75, c(.001,.6,.42,.58,.299)))))
#' df0 <- data.frame(t(rmultinom(60, 300,c(.001,.6,.2,.3,.299))))
#' metadata<-data.frame(f_s=factor(c(rep("A", 30), rep("B", 30))),
#'                      f_c=factor(c(rep("C", 15), rep("D", 15), rep("E", 15), rep("F", 15))),
#'                      age=c(1:30, 2:31)
#'                      )
#' table(metadata[, 'f_c'])
#' reg_res<-rf_reg.by_datasets(df, metadata, s_category='f_c', c_category='age')
#' rf_reg.cross_appl(reg_res, x_list=reg_res$x_list, y_list=reg_res$y_list)
#' reg_res<-rf_reg.by_datasets(df, metadata, nfolds=5, s_category='f_c', c_category='age')
#' rf_reg.cross_appl(reg_res, x_list=reg_res$x_list, y_list=reg_res$y_list)
#' @author Shi Huang
#' @export
rf_reg.cross_appl<-function(rf_list, x_list, y_list){
  ## TODO: check the class of inputs
  L<-length(rf_list$rf_model_list)
  try(if(!identical(L, length(x_list), length(y_list))) stop("The length of x list, y list and rf model list should be identical."))
  try(if(all(unlist(lapply(y_list, mode))!="numeric")) stop("All elements in the y list should be numeric for regression."))
  perf_summ<-data.frame(matrix(NA, ncol=14, nrow=L*L))
  colnames(perf_summ)<-c("Train_data", "Test_data", "Validation_type", "Sample_size", "Min_acutal_value", "Max_acutal_value", "Min_predicted_value", "Max_predicted_value",
                         "MSE", "RMSE", "MAE", "MAPE", "R_squared", "Adj_R_squared")
  predicted<-list()
  for(i in 1:L){
    y<-y_list[[i]]
    x<-x_list[[i]]
    try(if(nlevels(y)==1) stop("Less than one level in the subgroup for classification"))
    oob<-rf_list$rf_model_list[[i]]
    #---  RF Training performance: MSE, MAE and R_squared
    cat("\nTraining dataset: ", names(x_list)[i] ,"\n\n")
    train_sample_size<-length(y)
    train_y_min<-range(y)[1] #paste0(range(y), collapse = "-")
    train_y_max<-range(y)[2]
    train_pred_y_min<-range(oob$predictions)[1]
    train_pred_y_max<-range(oob$predictions)[2]#paste0(range(oob$predictions), collapse = "-")
    train_MSE<-rf_list$rf_MSE[[i]]
    train_RMSE<-rf_list$rf_RMSE[[i]]
    train_MAE<-rf_list$rf_MAE[[i]]
    train_MAPE<-rf_list$rf_MAPE[[i]]
    train_R_squared<-rf_list$rf_R_squared[[i]]
    train_Adj_R_squared<-rf_list$rf_Adj_R_squared[[i]]
    cat("MSE in the self-validation: ", train_MSE ,"\n")
    cat("RMSE in the self-validation: ", train_RMSE ,"\n")
    cat("MAE in the self-validation: ", train_MAE ,"\n")
    cat("MAE percentage in the self-validation: ", train_MAPE ,"\n")
    cat("R squared in the self-validation: ", train_R_squared ,"\n")
    cat("Adjusted R squared in the self-validation: ", train_Adj_R_squared ,"\n")
    D<-1:L
    T<-D[D!=i]
    a=1+(i-1)*L
    perf_summ[a, 1:3]<-c(names(x_list)[i], names(x_list)[i], "self_validation")
    perf_summ[a, 4:14]<-c(train_sample_size, train_y_min, train_y_max, train_pred_y_min, train_pred_y_max,
                          train_MSE, train_RMSE, train_MAE, train_MAPE, train_R_squared, train_Adj_R_squared)
    predicted[[a]]<-data.frame(test_y=y, pred_y=oob$predictions)
    names(predicted)[a]<-paste(names(x_list)[i], names(x_list)[i], sep="__VS__")
    loop_num<-1
    for(j in T){
      if(nrow(x_list[[j]])>0){
        newx<-x_list[[j]]
        newy<-y_list[[j]]
        if(class(oob)=="ranger"){oob_rf.model<-oob
          }else{
           oob_rf.model<-oob[[1]]} # pick one out of K rf models in the ranger list
        pred_newy<-predict(oob_rf.model, newx, type="response")$predictions # ranger only
        #---  RF test performance: MSE, MAE and R_squared
        cat("Test dataset: ", names(x_list)[j] ,"\n")
        test_sample_size<-length(newy)
        test_y_min<-range(newy)[1] #paste0(range(newy), collapse = "-")
        test_y_max<-range(newy)[2]
        test_pred_y_min<-range(pred_newy)[1]
        test_pred_y_max<-range(pred_newy)[2]#paste0(range(pred_newy), collapse = "-")
        MSE<-function(y, pred_y){ mean((y-pred_y)^2)}
        RMSE<-function(y, pred_y){ sqrt(mean((y-pred_y)^2))}
        MAE<-function(y, pred_y){ mean(sqrt((y-pred_y)^2))}
        R2<-function(y, pred_y){ 1-(sum((y-pred_y)^2) / sum((y-mean(y))^2)) }
        MAPE<-function(y, pred_y){ mean(sqrt((y-pred_y)^2)/y)}
        adj.R2<-function(y, pred_y, k){
          n=length(y);
          1-(1-R2(y, pred_y)^2)*(n-1)/(n-k-1) # k is # of predictors
        }
        test_MSE<-MSE(newy, pred_newy)
        test_RMSE<-RMSE(newy, pred_newy)
        test_MAE<-MAE(newy, pred_newy)
        test_MAPE<-MAPE(newy, pred_newy)
        test_R_squared<-R2(newy, pred_newy)
        test_Adj_R_squared<-adj.R2(newy, pred_newy, k=ncol(newx))
        cat("MSE in the cross-applications: ", test_MSE ,"\n")
        cat("RMSE in the cross-applications: ", test_RMSE ,"\n")
        cat("MAE in the cross-applications: ", test_MAE ,"\n")
        cat("MAE percentage in the cross-applications: ", test_MAPE ,"\n")
        cat("R squared in the cross-application: ", test_R_squared ,"\n")
        cat("Adjusted R squared in the cross-application: ", test_Adj_R_squared ,"\n")
        perf_summ[a+loop_num, 1:3]<-c(names(x_list)[i], names(x_list)[j], "cross_application")
        perf_summ[a+loop_num, 4:14]<-c(test_sample_size, test_y_min, test_y_max, test_pred_y_min, test_pred_y_max,
                                       test_MSE, test_RMSE, test_MAE, test_MAPE, test_R_squared, test_Adj_R_squared)
        predicted[[a+loop_num]]<-data.frame(test_y=newy, pred_y=pred_newy)
        names(predicted)[a+loop_num]<-paste(names(x_list)[i], names(x_list)[j], sep="__VS__")
        loop_num<-loop_num+1
      }
    }
  }
  res<-list()
  res$perf_summ<-perf_summ
  res$predicted<-predicted
  class(res)<-"rf_reg.cross_appl"
  res
}

#' @title generate.comps_datalist
#' @description To generate sub-datasets and sub-metadata by pairing
#' one level and each of other levels of the specified category in the metadata.
#' @param df Training data: a data.frame.
#' @param f A factor in the metadata with at least two levels (groups).
#' @param comp_group A string indicates the group in the f
#' @return A object of list with elements: a list of data.frame and a list of factor
#' @author Shi Huang
#'
generate.comps_datalist<-function(df, f, comp_group){
  all_other_groups<-levels(f)[which(levels(f)!=comp_group)]
  L<-length(all_other_groups)
  f_list<-list()
  df_list<-list()
  for(i in 1:L){
    f_list[[i]]<-factor(f[which(f==comp_group | f==all_other_groups[i])])
    df_list[[i]]<-df[which(f==comp_group | f==all_other_groups[i]), ]
  }
  names(f_list)<-names(df_list)<-paste(comp_group, all_other_groups, sep="_VS_")
  res<-list()
  res$f_list<-f_list
  res$df_list<-df_list
  res
}

#' @title rf_clf.comps
#' @description It runs standard random forests with oob estimation for classification of
#' one level VS all other levels of one category in the datasets.
#' The output includes a list of rf models for the sub datasets
#' and all important statistics for each of features.
#'
#' @param df Training data: a data.frame.
#' @param f A factor in the metadata with at least two levels (groups).
#' @param comp_group A string indicates the group in the f
#' @param clr_transform A boolean value indicating if the clr-transformation applied.
#' @param rf_imp_values A boolean value indicating if compute both importance score and pvalue for each feature.
#' @param verbose A boolean value indicating if show computation status and estimated runtime.
#' @param ntree The number of trees.
#' @param p.adj.method The p-value correction method, default is "bonferroni".
#' @param q_cutoff The cutoff of q values for features, the default value is 0.05.
#' @return ...
#' @seealso ranger
#' @examples
#' df0 <- data.frame(t(rmultinom(60, 300,c(.001,.6,.2,.3,.299))))
#' df <- data.frame(rbind(t(rmultinom(7, 75, c(.21,.6,.12,.38,.099))),
#'             t(rmultinom(8, 75, c(.001,.6,.42,.58,.299))),
#'             t(rmultinom(15, 75, c(.011,.6,.22,.28,.289))),
#'             t(rmultinom(15, 75, c(.091,.6,.32,.18,.209))),
#'             t(rmultinom(15, 75, c(.001,.6,.42,.58,.299)))))
#' f=factor(c(rep("A", 15), rep("B", 15), rep("C", 15), rep("D", 15)))
#' comp_group="A"
#' comps_res<-rf_clf.comps(df, f, comp_group, verbose=FALSE, ntree=500,
#'                         p.adj.method = "bonferroni", q_cutoff=0.05)
#' comps_res
#' @author Shi Huang
#' @export
rf_clf.comps<-function(df, f, comp_group, verbose=FALSE, clr_transform=TRUE,
                       rf_imp_values=FALSE,ntree=500, p.adj.method = "bonferroni",
                       q_cutoff=0.05){
  f<-factor(f)
  all_other_groups<-levels(f)[which(levels(f)!=comp_group)]
  L<-length(all_other_groups)
  nCores <- parallel::detectCores()
  doMC::registerDoMC(nCores)
  # comb function for parallelization using foreach
  comb <- function(x, ...) {
    lapply(seq_along(x),
           function(i) c(x[[i]], lapply(list(...), function(y) y[[i]])))
  }
  require('foreach')
  oper<-foreach::foreach(i=1:L, .combine='comb', .multicombine=TRUE,
                         .init=list(list(), list(), list(), list(), list(), list(), list())) %dopar% {
    sub_f<-factor(f[which(f==comp_group | f==all_other_groups[i])])
    sub_df<-df[which(f==comp_group | f==all_other_groups[i]), ]
    print(levels(factor(sub_f)))
    # 1. sample size of all datasets
    sample_size<-length(factor(sub_f))
    dataset<-paste(comp_group, all_other_groups[i], sep="_VS_")
    # 2. AUC of random forest model
    oob <- rf.out.of.bag(sub_df, factor(sub_f), verbose=verbose, ntree=ntree, imp_pvalues = rf_imp_values)
    if(nlevels(factor(sub_f))==2){rf_AUC <- get.auroc(oob$probabilities, factor(sub_f), comp_group)}else{rf_AUC <- NA}
    # 3. # of significantly differential abundant features between health and disease
    out<-BetweenGroup.test(sub_df, factor(sub_f), clr_transform=clr_transform, q_cutoff=q_cutoff,
                           positive_class=comp_group, p.adj.method = p.adj.method)
    wilcox<-data.frame(feature=rownames(out),
                       dataset=rep(dataset, ncol(sub_df)), rf_imps=oob$importances,
                       out)
    list(x=sub_df, y=sub_f, sample_size=sample_size, datasets=dataset, oob=oob, rf_AUC=rf_AUC, wilcox=wilcox)
  }
  names(oper[[1]])<-names(oper[[2]])<-names(oper[[5]])<-unlist(oper[[4]])
  result<-list()
  result$x_list<-oper[[1]]
  result$y_list<-oper[[2]]
  result$sample_size<-unlist(oper[[3]])
  result$datasets<-unlist(oper[[4]])
  result$rf_model_list<-oper[[5]]
  result$rf_AUC<-unlist(oper[[6]])
  result$feature_imps_list<-oper[[7]]
  class(result)<-"rf_clf.comps"
  return(result)
}


