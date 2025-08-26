#' @importFrom stats median
#' @importFrom pROC plot.roc ci.se
#' @importFrom grDevices dev.off hcl pdf
#' @importFrom graphics plot text
#' @importFrom rlang .data
#' @importFrom utils data write.table
#' @importFrom ggplot2 ggplot aes xlab ylab theme_bw coord_flip ggsave geom_boxplot geom_histogram geom_violin geom_jitter geom_point geom_line geom_hline geom_vline xlim ylim scale_x_continuous scale_color_manual theme element_rect element_line element_text element_blank
#' @importFrom PRROC roc.curve pr.curve
#' @importFrom gridExtra arrangeGrob
#' @importFrom reshape2 melt
#' @importFrom gtools rdirichlet
#' @import viridis

#' @title plot_clf_pROC
#' @param y A factor of classes to be used as the true results
#' @param rf_clf_model A list object of a random forest model.
# @param predictor A numeric or ordered vector of the same length than response, containing the predicted value of each observation.
#' @param positive_class an optional character string for the factor level that corresponds to a "positive" result (if that makes sense for your data).
#' If there are only two factor levels, the first level will be used as the "positive" result.
#' @param prefix The prefix of data set.
#' @param outdir The output directory.
#' @examples
#' set.seed(123)
#' x <- data.frame(rbind(t(rmultinom(7, 75, c(.201,.5,.02,.18,.099))),
#'             t(rmultinom(8, 75, c(.201,.4,.12,.18,.099))),
#'             t(rmultinom(15, 75, c(.011,.3,.22,.18,.289))),
#'             t(rmultinom(15, 75, c(.091,.2,.32,.18,.209))),
#'             t(rmultinom(15, 75, c(.001,.1,.42,.18,.299)))))
#' y<-factor(c(rep("A", 30), rep("C", 30)))
#' rf_clf_model<-rf.out.of.bag(x, y)
#' positive_class="A"
#' plot_clf_pROC(y, rf_clf_model, outdir='./')
#' @author Shi Huang
#' @export
plot_clf_pROC<-function(y, rf_clf_model, positive_class=NA, prefix="train", outdir=NULL){
  if(class(rf_clf_model)!="rf.out.of.bag" & class(rf_clf_model)!="rf.cross.validation")stop("The rf_clf_model is not rf.out.of.bag/rf.cross.validation.")
  if(nlevels(y)!=2) stop("pROC only support for ROC analysis in the binary classification!")
  positive_class<-ifelse(is.na(positive_class), levels(y)[1], positive_class)
  predictor<-rf_clf_model$probabilities[, positive_class]
  rocobj <- pROC::plot.roc(y, predictor, main="", percent=TRUE,ci=TRUE) # print the AUC (will contain the CI)
  ciobj <- pROC::ci.se(rocobj, specificities=seq(0, 100, 5)) # over a select set of specificities
  if(!is.null(outdir)){
  pdf(paste(outdir, prefix, ".rf_clf_pROC.ci.pdf",sep=""), width=4, height=4)
  plot(ciobj, type="shape", col="#1c61b6AA", print.auc=TRUE) # plot as a blue shape
  text(70,25, paste0(levels(y), collapse = " VS "), pos=4)
  text(70,15, paste("AUC = ",formatC(rocobj$auc,digits=2,format="f"),sep=""),pos=4)
  ci.lower<-formatC(rocobj$ci[1],digits=2,format="f")
  ci.upper<-formatC(rocobj$ci[3],digits=2,format="f")
  text(70,5, paste("95% CI: ",ci.lower,"-",ci.upper,sep=""),pos=4)
  dev.off()
  }
  result<-list()
  result$rocobj<-rocobj
  result$ciobj<-ciobj
  invisible(result)
}
#' @title plot_PRC
#' @param y A factor of classes to be used as the true results
#' @param rf_clf_model A list object of a random forest model.
# @param predictor A numeric or ordered vector of the same length than response, containing the predicted value of each observation.
#' @param positive_class an optional character string for the factor level that corresponds to a "positive" result (if that makes sense for your data).
#' If there are only two factor levels, the first level will be used as the "positive" result.
#' @param prefix The prefix of data set.
#' @param outdir The output directory.
#' @examples
#' set.seed(123)
#' x <- data.frame(rbind(t(rmultinom(7, 75, c(.201,.5,.02,.18,.099))),
#'             t(rmultinom(8, 75, c(.201,.4,.12,.18,.099))),
#'             t(rmultinom(15, 75, c(.011,.3,.22,.18,.289))),
#'             t(rmultinom(15, 75, c(.091,.2,.32,.18,.209))),
#'             t(rmultinom(15, 75, c(.001,.1,.42,.18,.299)))))
#' x0 <- data.frame(rbind(t(rmultinom(7, 75, c(.011,.3,.22,.18,.289))),
#'             t(rmultinom(8, 75, c(.011,.3,.22,.18,.289))),
#'             t(rmultinom(15, 75, c(.011,.3,.22,.18,.289))),
#'             t(rmultinom(15, 75, c(.011,.3,.22,.18,.289))),
#'             t(rmultinom(15, 75, c(.011,.3,.22,.18,.289)))))
#' y<-factor(c(rep("A", 30), rep("B", 30)))
#' y0<-factor(c(rep("A", 5), rep("B", 55)))
#' rf_clf_model<-rf.out.of.bag(x, y)
#' plot_clf_PRC(y, rf_clf_model, positive_class="A", outdir='./A')
#' plot_clf_PRC(y, rf_clf_model, positive_class="B", outdir='./B')
#' rf_clf_model0<-rf.out.of.bag(x0, y0)
#' plot_clf_PRC(y0, rf_clf_model0, positive_class="A", outdir='./A')
#' plot_clf_PRC(y0, rf_clf_model0, positive_class="B", outdir='./B')
#' @author Shi Huang
#' @export
plot_clf_PRC<-function(y, rf_clf_model, positive_class=NA, prefix="train", outdir=NULL){
  #require('PRROC')
  if(class(rf_clf_model)!="rf.out.of.bag" & class(rf_clf_model)!="rf.cross.validation")stop("The rf_clf_model is not rf.out.of.bag/rf.cross.validation.")
  # if(nlevels(y)!=2) stop("PPROC only support for PRC analysis in the binary classification!")
  positive_class<-ifelse(is.na(positive_class), levels(y)[1], positive_class)
  predictor<-rf_clf_model$probabilities[, positive_class]
  df<-data.frame(y, predictor)
  prob_pos<-df[df$y==positive_class, "predictor"]
  prob_neg<-df[df$y!=positive_class,"predictor"]
  pr<-pr.curve(prob_pos, prob_neg, curve=TRUE, max.compute = TRUE, min.compute = TRUE, rand.compute = TRUE)
  auprc<-round(pr$auc.integral, 3)
  rel_auprc<-round((pr$auc.integral-pr$min$auc.integral)/(pr$max$auc.integral-pr$min$auc.integral), 3)
  if(!is.null(outdir)){
    pdf(paste(outdir, prefix, ".rf_clf_PRC.pdf",sep=""), width=5, height=5)
    plot(pr, max.plot = TRUE, min.plot = TRUE, rand.plot = TRUE, fill.area = TRUE,
         auc.main = FALSE, main = paste("AUPRC=",auprc,"\n", "Relative AUPRC=", rel_auprc, sep=""))
    dev.off()
  }
  invisible(pr)
}

#' @title plot_ROC
#' @param y A factor of classes to be used as the true results
#' @param rf_clf_model A list object of a random forest model.
# @param predictor A numeric or ordered vector of the same length than response, containing the predicted value of each observation.
#' @param positive_class an optional character string for the factor level that corresponds to a "positive" result (if that makes sense for your data).
#' If there are only two factor levels, the first level will be used as the "positive" result.
#' @param prefix The prefix of data set.
#' @param outdir The output directory.
#' @examples
#' set.seed(123)
#' x <- data.frame(rbind(t(rmultinom(7, 75, c(.201,.5,.02,.18,.099))),
#'             t(rmultinom(8, 75, c(.201,.4,.12,.18,.099))),
#'             t(rmultinom(15, 75, c(.011,.3,.22,.18,.289))),
#'             t(rmultinom(15, 75, c(.091,.2,.32,.18,.209))),
#'             t(rmultinom(15, 75, c(.001,.1,.42,.18,.299)))))
#' x0 <- data.frame(rbind(t(rmultinom(7, 75, c(.011,.3,.22,.18,.289))),
#'             t(rmultinom(8, 75, c(.011,.3,.22,.18,.289))),
#'             t(rmultinom(15, 75, c(.011,.3,.22,.18,.289))),
#'             t(rmultinom(15, 75, c(.011,.3,.22,.18,.289))),
#'             t(rmultinom(15, 75, c(.011,.3,.22,.18,.289)))))
#' y<-factor(c(rep("A", 20), rep("B", 20), rep("C", 20)))
#' y0<-factor(c(rep("A", 5), rep("B", 55)))
#' rf_clf_model<-rf.out.of.bag(x, y)
#' plot_clf_ROC(y, rf_clf_model, positive_class="A", outdir='./A')
#' plot_clf_ROC(y, rf_clf_model, positive_class="B", outdir='./B')
#' plot_clf_PRC(y, rf_clf_model, positive_class="A", outdir='./A')
#' plot_clf_PRC(y, rf_clf_model, positive_class="B", outdir='./B')
#' rf_clf_model0<-rf.out.of.bag(x0, y0)
#' plot_clf_ROC(y0, rf_clf_model0, positive_class="A", outdir='./A')
#' plot_clf_ROC(y0, rf_clf_model0, positive_class="B", outdir='./B')
#' plot_clf_PRC(y0, rf_clf_model0, positive_class="A", outdir='./A')
#' plot_clf_PRC(y0, rf_clf_model0, positive_class="B", outdir='./B')
#' pred_df<-data.frame(y=y0, prediction=rf_clf_model0$probabilities[, positive_class])
#' ggplot(pred_df, aes(prediction, fill = y)) +
#' geom_histogram(alpha = 0.5, position = 'identity')
#' @author Shi Huang
#' @export
plot_clf_ROC<-function(y, rf_clf_model, positive_class=NA, prefix="train", outdir=NULL){
  if(class(rf_clf_model)!="rf.out.of.bag" & class(rf_clf_model)!="rf.cross.validation")stop("The rf_clf_model is not rf.out.of.bag/rf.cross.validation.")
  #require('PRROC')
  #if(nlevels(y)!=2) stop("PPROC only support for ROC analysis in the binary classification!")
  positive_class<-ifelse(is.na(positive_class), levels(y)[1], positive_class)
  predictor<-rf_clf_model$probabilities[, positive_class]
  df<-data.frame(y, predictor)
  prob_pos<-df[df$y==positive_class, "predictor"]
  prob_neg<-df[df$y!=positive_class,"predictor"]
  roc<-roc.curve(prob_pos, prob_neg, curve=TRUE, max.compute = TRUE, min.compute = TRUE, rand.compute = TRUE)
  if(!is.null(outdir)){
    pdf(paste(outdir, prefix, ".rf_clf_ROC.pdf",sep=""), width=5, height=5)
    plot(roc, rand.plot = TRUE, fill.area = TRUE)
    dev.off()
  }
  invisible(roc)
}



#' @title plot_clf_probabilities
#' @param y A factor of classes to be used as the true results
#' @param rf_clf_model The rf classification model from \code{rf.out.of.bag}
#' @param positive_class an optional character string for the factor level that corresponds to a "positive" result (if that makes sense for your data).
#' If there are only two factor levels, the first level will be used as the "positive" result.
#' @param prefix The prefix of data set.
#' @param outdir The output directory.
#' @examples
#' set.seed(123)
#' x <- data.frame(rbind(t(rmultinom(7, 75, c(.201,.5,.02,.18,.099))),
#'             t(rmultinom(8, 75, c(.201,.4,.12,.18,.099))),
#'             t(rmultinom(15, 75, c(.011,.3,.22,.18,.289))),
#'             t(rmultinom(15, 75, c(.091,.2,.32,.18,.209))),
#'             t(rmultinom(15, 75, c(.001,.1,.42,.18,.299)))))
#' y<-factor(c(rep("1", 20), rep("A", 20), rep("C", 20)))
#' rf_clf_model<-rf.out.of.bag(x, y)
#' positive_class="1"
#' plot_clf_probabilities(y, rf_clf_model, positive_class, outdir='./')
#' @author Shi Huang
#' @export
plot_clf_probabilities<-function(y, rf_clf_model, positive_class=NA, prefix="train", outdir=NULL){
  if(class(rf_clf_model)!="rf.out.of.bag" & class(rf_clf_model)!="rf.cross.validation")stop("The rf_clf_model is not rf.out.of.bag/rf.cross.validation.")
  positive_class<-ifelse(is.na(positive_class), levels(y)[1], positive_class)
  l<-levels(y); l_sorted<-sort(levels(y))
  Mycolor <- rep(c("#D55E00", "#0072B2"), length.out=length(l))
  if(identical(order(l), order(l_sorted))){
    Mycolor=Mycolor; l_ordered=l
  }else{Mycolor=rev(Mycolor); l_ordered=l_sorted}
  y_prob<-data.frame(y, predictor=rf_clf_model$probabilities[,positive_class])
  p<-ggplot(y_prob, aes(x=.data$y, y=.data$predictor)) +
    geom_violin()+
    geom_jitter(position=position_jitter(width=0.2), alpha=0.1) +
    geom_boxplot(outlier.shape = NA, width=0.4, alpha=0.01)+
    geom_hline(yintercept=0.5, linetype="dashed")+
    ylim(0, 1)+
    theme_bw()+
    #ggtitle(paste("Wilcoxon Rank Sum Test:\n P=",p_mf,sep=""))+
    xlab("") +
    ylab(paste("Probability of ", positive_class))+
    theme(legend.position="none")+
    theme(axis.line = element_line(color="black"),
          strip.background = element_rect(colour = "white"),
          panel.border = element_blank())+
    scale_color_manual(values = Mycolor, labels=l_ordered)
  if(!is.null(outdir)){
  ggsave(filename=paste(outdir,prefix,".probability_",positive_class,".boxplot.pdf",sep=""), plot=p, width=3, height=4)
  }
}

#' @title plot_topN_imp_scores
#' @param rf_model The rf model from \code{rf.out.of.bag} or \code{rf.cross.validation}
#' @param topN A number indicating how many top important feature need to visualize in the barplot.
#' @param feature_md an optional data.frame including feature IDs and feature annotations.
#' @param feature_id_col The Feature_ID column in the feature metadata.
#' @param bar_color The column name in the feature metadata for coloring the bar plot.
#' @param plot_width plotting parameter.
#' @param plot_height plotting parameter.
#' @param outdir The output directory.
#' @examples
#' set.seed(123)
#' x <- data.frame(rbind(t(rmultinom(7, 75, c(.201,.5,.02,.18,.099))),
#'             t(rmultinom(8, 75, c(.201,.4,.12,.18,.099))),
#'             t(rmultinom(15, 75, c(.011,.3,.22,.18,.289))),
#'             t(rmultinom(15, 75, c(.091,.2,.32,.18,.209))),
#'             t(rmultinom(15, 75, c(.001,.1,.42,.18,.299)))))
#' y<-factor(c(rep("1", 20), rep("A", 20), rep("C", 20)))
#' feature_md <- data.frame(Feature_ID=c("X1", "X2", "X3", "X4", "X5"),
#'                          Taxon=c("AA", "Bd", "CC", "CC", "AA"))
#' rf_model<-rf.out.of.bag(x, y)
#' rf_model<-rf.cross.validation(x, y, nfolds=5)
#' plot_topN_imp_scores(rf_model, topN=4, feature_md, outdir=NULL, bar_color="Taxon")
#' @author Shi Huang
#' @export
plot_topN_imp_scores<-function(rf_model, topN=4, feature_md=NULL, feature_id_col="Feature_ID", bar_color=NA,
                               outdir=NULL, plot_height=8, plot_width=5){
  if(class(rf_model)!="rf.out.of.bag" & class(rf_model)!="rf.cross.validation")
    stop("The rf_model is not rf.out.of.bag/rf.cross.validation.")
  if(class(rf_model)=="rf.cross.validation"){
    rank_mat <- apply(rf_model$importances, 2, function(x){rank(-x, na.last = "keep")})
    rf_imp_rank <- rank(apply(rank_mat, 1, median), na.last = "keep")
    imps_df<-data.frame(Feature_ID=rownames(rf_model$importances),
                     rf_model$importances,
                     Imps=apply(rf_model$importances, 1, median),
                     ImpsSD=apply(rf_model$importances, 1, sd),
                     ImpsSE=apply(rf_model$importances, 1, sd)/sqrt(ncol(rf_model$importances)),
                     Rank=rf_imp_rank)
  }else{
    imps_df<-data.frame(Feature_ID=rownames(rf_model$importances),
                    Imps=rf_model$importances,
                    ImpsSD=0,
                    ImpsSE=0,
                    Rank=rank(-rf_model$importances, na.last = "keep"))
  }
  # add feature metadata
  merge_feature_md<-function(df, feature_md, df_id_col=1, fmd_id_col=1){
    matched_idx<-which(feature_md[, fmd_id_col] %in% df[, df_id_col])
    uniq_features_len<-length(unique(df[, df_id_col]))
    if(uniq_features_len  > length(matched_idx)){
      warning("# of features has no matching IDs in the feature metadata file: ", uniq_features_len-length(matched_idx), "\n")
    }
    feature_md_matched<-feature_md[matched_idx,]
    out<-merge(df, feature_md_matched, by.x=df_id_col, by.y=fmd_id_col)
    out
  }
  if(!is.null(feature_md)){
    imps_df<- merge_feature_md(imps_df, feature_md, df_id_col = "Feature_ID", fmd_id_col = feature_id_col)
  }
  imps_df <- imps_df[order(imps_df$Rank), ]
  top_n_imps_df<-imps_df[1:topN, ]
  top_n_imps_df<-top_n_imps_df[order(top_n_imps_df$Imps, decreasing = T), ]
  # bar plot
  p <- ggplot(top_n_imps_df, aes(x=reorder(top_n_imps_df$Feature_ID, top_n_imps_df$Imps), y=top_n_imps_df$Imps)) +
    geom_bar(stat="identity", alpha=0.5) +
    scale_fill_viridis(discrete=TRUE) +
    ylab("RF importance score") +
    xlab("Feature ID") +
    ylim(0, max(top_n_imps_df$Imps, na.rm=T)) +
    theme_minimal() +
    coord_flip()
  if(!is.na(bar_color)){
    p <-p + geom_bar(data=top_n_imps_df, aes(fill=get(bar_color)), stat="identity", alpha=0.5) + labs(fill=bar_color)
  }
  if(!is.null(outdir)){
    ggsave(filename=paste(outdir, "Top_", topN,".imps.barplot.pdf",sep=""), plot=p, width=plot_width, height=plot_height)
  }

  res <- list()
  res$imps_df <- imps_df
  res$plot <- p
  res

}

#' @title plot_clf_feature_selection
#' @description Plot the classification performance against the gradually reduced number of features used in the modeling.
#' @param x The data frame or data matrix for model training.
#' @param y A factor related to the responsive vector for training data.
#' @param nfolds The number of folds in the cross-validation for each feature set.
#' @param rf_clf_model The rf classification model from \code{rf.out.of.bag}
#' @param positive_class A class of the y.
#' @param metric The classification performance metric applied.
#' If binary classification, this must be one of "AUROC", "Accuracy", "Kappa", "F1".
#' If multi-class classification, this must be one of "Accuracy", "Kappa".
#' @param outdir The output directory.
#' @examples
#' set.seed(123)
#' require("gtools")
#' n_features <- 100
#' prob_vec <- rdirichlet(5, sample(n_features))
#' x <- data.frame(rbind(t(rmultinom(7, 7*n_features, prob_vec[1, ])),
#'             t(rmultinom(8, 8*n_features, prob_vec[2, ])),
#'             t(rmultinom(15, 15*n_features, prob_vec[3, ])),
#'             t(rmultinom(15, 15*n_features, prob_vec[4, ])),
#'             t(rmultinom(15, 15*n_features, prob_vec[5, ]))))
#' y<-factor(c(rep("A", 30), rep("C", 30)))
#' s<-factor(rep(c("B1", "B2", "B3", "B4"), 15))
#' rf_model<-rf.cross.validation(x, y, nfolds=5)
#' summ <- plot_clf_feature_selection(x, y, nfolds=5, rf_model, metric="AUROC", outdir=NULL)
#' summ <- plot_clf_feature_selection(x, y, nfolds=s, rf_model, metric="AUROC", outdir=NULL)
#' res<-replicate(10, plot_clf_feature_selection(x, y,
#' nfolds=5, rf_model, metric="AUROC", outdir=NULL))
#' do.call(rbind, res["top_n_perf", ])
#' @author Shi Huang
#' @export
plot_clf_feature_selection <- function(x, y, nfolds=5, rf_clf_model, #log10_n_features=FALSE,
                                       metric="AUROC", positive_class=NA, outdir=NULL){
  if(class(rf_clf_model)=="rf.cross.validation"){
    rank_mat <- apply(rf_clf_model$importances, 2, function(x){rank(-x, na.last = "keep")})
    rf_imp_rank <- rank(apply(rank_mat, 1, median), na.last = "keep")
  }else if(class(rf_clf_model)=="rf.out.of.bag"){
    rf_imp_rank<-rank(-(rf_clf_model$importances), na.last = "keep")
  }else{
    stop("The class of input rf model should be rf.out.of.bag or rf.cross.validation.")
  }
  positive_class<-ifelse(is.na(positive_class), levels(y)[1], positive_class)
  max_n<-max(rf_imp_rank, na.rm = TRUE)
  n_total_features<-ncol(x)
  n_features<-c(2, 4, 8, 16, 32, 64, 128, 256, 512, 1024)
  min_dist_idx <- which.min(abs(n_features-n_total_features))
  maginal_idx <- ifelse(n_features[min_dist_idx]>n_total_features, min_dist_idx-1, min_dist_idx)
  n_features<-n_features[1:maginal_idx]
  top_n_perf<-matrix(NA, ncol=6, nrow=length(n_features)+1)
  colnames(top_n_perf)<-c("n_features", "AUROC", "AUPRC", "Accuracy", "Kappa", "F1")
  rownames(top_n_perf)<-top_n_perf[,1]<-c(n_features, max_n)
  top_n_rf_list <-vector(mode="list", length = length(n_features))
  names(top_n_rf_list) <- n_features
  for(i in 1:length(n_features)){
    idx<-which(rf_imp_rank<=n_features[i])
    #top_n_features<-names(rf_imp_rank[idx])
    x_n<-x[, idx]
    y_n<-y
    top_n_rf_list[[i]] <- top_n_rf <- rf.cross.validation(x_n, y_n, nfolds=nfolds, ntree=500)
    top_n_AUROC<-plot_clf_ROC(y_n, top_n_rf, positive_class=positive_class)$auc
    top_n_AUPRC<-plot_clf_PRC(y_n, top_n_rf, positive_class=positive_class)$auc.integral
    top_n_conf<-caret::confusionMatrix(data=top_n_rf$predicted, top_n_rf$y, positive=positive_class)
    top_n_perf[i, 1]<-n_features[i]
    top_n_perf[i, 2]<-top_n_AUROC
    top_n_perf[i, 3]<-top_n_AUPRC
    top_n_perf[i, 4]<-top_n_conf$overall[1] # Accuracy
    top_n_perf[i, 5]<-top_n_conf$overall[2]# kappa conf$byClass["F1"]
    top_n_perf[i, 6]<-top_n_conf$byClass["F1"]
  }

  all_AUROC<-plot_clf_ROC(y, rf_clf_model, positive_class=positive_class)$auc
  all_AUPRC<-plot_clf_PRC(y, rf_clf_model, positive_class=positive_class)$auc.integral
  all_conf<-caret::confusionMatrix(data=rf_clf_model$predicted, rf_clf_model$y, positive=positive_class)
  top_n_perf[length(n_features)+1, ]<-c(max_n, all_AUROC, all_AUPRC, all_conf$overall[1],
                                        all_conf$overall[2], all_conf$byClass["F1"])
  top_n_perf<-data.frame(top_n_perf)
  best_idx <- which.max(top_n_perf[, metric])
  best_n_features <- top_n_perf$n_features[best_idx]
  breaks<-top_n_perf$n_features
  p<-ggplot(top_n_perf, aes(x=.data$n_features, y=get(metric))) +
    xlab("# of features used")+
    ylab(metric)+
    scale_x_continuous(trans = "log10", breaks=breaks)+
    geom_point() + geom_line()+ ylim(0.5, 1)+
    geom_vline(xintercept = best_n_features, linetype="dashed", color="blue")+
    geom_text(aes(x=best_n_features, y=top_n_perf[best_idx, metric], label=best_n_features), color="blue", vjust=-1)+
    theme_bw()+
    theme(axis.line = element_line(color="black"),
          axis.title = element_text(size=18),
          strip.background = element_rect(colour = "white"),
          panel.border = element_blank())
  if(!is.null(outdir)){
  ggsave(filename=paste(outdir,"train.",metric,"_VS_top_ranking_features.scatterplot.pdf",sep=""), plot=p, width=5, height=4)
  }
  res <- list()
  res$top_n_perf <- top_n_perf
  res$best_n_features <- best_n_features
  res$top_n_rf <-top_n_rf_list
  res$plot <- p
  res
}




