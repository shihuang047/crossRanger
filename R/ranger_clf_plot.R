#' @importFrom pROC plot.roc ci.se
#' @importFrom grDevices dev.off hcl pdf
#' @importFrom graphics plot text
#' @importFrom utils data write.table
#' @import ggplot2
#' @import PRROC
#' @importFrom gridExtra arrangeGrob
#' @importFrom reshape2 melt

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
#' y<-factor(c(rep("A", 30), rep("C", 30)))
#' rf_clf_model<-rf.out.of.bag(x, y)
#' positive_class="A"
#' plot_clf_PRC(y, rf_clf_model, outdir='./')
#' @author Shi Huang
#' @export
plot_clf_PRC<-function(y, rf_clf_model, positive_class=NA, prefix="train", outdir=NULL){
  require('PRROC')
  if(class(rf_clf_model)!="rf.out.of.bag" & class(rf_clf_model)!="rf.cross.validation")stop("The rf_clf_model is not rf.out.of.bag/rf.cross.validation.")
  if(nlevels(y)!=2) stop("PPROC only support for PRC analysis in the binary classification!")
  positive_class<-ifelse(is.na(positive_class), levels(y)[1], positive_class)
  predictor<-rf_clf_model$probabilities[, positive_class]
  df<-data.frame(y, predictor)
  prob_pos<-df[df$y==positive_class, "predictor"]
  prob_neg<-df[df$y!=positive_class,"predictor"]
  pr<-pr.curve(prob_pos, prob_neg, curve=TRUE)
  if(!is.null(outdir)){
    pdf(paste(outdir, prefix, ".rf_clf_PRC.pdf",sep=""), width=5, height=5)
    plot(pr)
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
#' y<-factor(c(rep("A", 30), rep("C", 30)))
#' rf_clf_model<-rf.out.of.bag(x, y)
#' positive_class="A"
#' plot_clf_ROC(y, rf_clf_model, outdir='./')
#' @author Shi Huang
#' @export
plot_clf_ROC<-function(y, rf_clf_model, positive_class=NA, prefix="train", outdir=NULL){
  if(class(rf_clf_model)!="rf.out.of.bag" & class(rf_clf_model)!="rf.cross.validation")stop("The rf_clf_model is not rf.out.of.bag/rf.cross.validation.")
  require('PRROC')
  if(nlevels(y)!=2) stop("PPROC only support for ROC analysis in the binary classification!")
  positive_class<-ifelse(is.na(positive_class), levels(y)[1], positive_class)
  predictor<-rf_clf_model$probabilities[, positive_class]
  df<-data.frame(y, predictor)
  prob_pos<-df[df$y==positive_class, "predictor"]
  prob_neg<-df[df$y!=positive_class,"predictor"]
  roc<-roc.curve(prob_pos, prob_neg, curve=TRUE)
  if(!is.null(outdir)){
    pdf(paste(outdir, prefix, ".rf_clf_ROC.pdf",sep=""), width=5, height=5)
    plot(roc)
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
#' y<-factor(c(rep("1", 30), rep("C", 30)))
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
  p<-ggplot(y_prob, aes(x=y, y=predictor)) +
    geom_violin()+
    geom_jitter(position=position_jitter(width=0.2),alpha=0.1) +
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

#' @title plot_clf_feature_selection
#' @description Plot the classification performance against the gradually reduced number of features used in the modeling.
#' @param x The data frame or data matrix for model training.
#' @param y A factor related to the responsive vector for training data.
#' @param rf_clf_model The rf classification model from \code{rf.out.of.bag}
#' @param positive_class A class of the y.
#' @param metric The classification performance metric applied.
#' If binary classification, this must be one of "AUROC", "Accuracy", "Kappa", "F1".
#' If multi-class classification, this must be one of "Accuracy", "Kappa".
#' @param outdir The output directory.
#' @examples
#' set.seed(123)
#' x <- data.frame(rbind(t(rmultinom(7, 75, c(.201,.5,.02,.18,.099))),
#'             t(rmultinom(8, 75, c(.201,.4,.12,.18,.099))),
#'             t(rmultinom(15, 75, c(.011,.3,.22,.18,.289))),
#'             t(rmultinom(15, 75, c(.091,.2,.32,.18,.209))),
#'             t(rmultinom(15, 75, c(.001,.1,.42,.18,.299)))))
#' y<-factor(c(rep("A", 30), rep("C", 30)))
#' rf_model<-rf.out.of.bag(x, y)
#' plot_clf_feature_selection(x, y, rf_model, metric="AUROC", outdir=NULL)
#' @author Shi Huang
#' @export
plot_clf_feature_selection <- function(x, y, rf_clf_model, metric="AUROC", positive_class=NA, outdir=NULL){
  if(class(rf_clf_model)!="rf.out.of.bag" & class(rf_clf_model)!="rf.cross.validation")stop("The rf_clf_model is not rf.out.of.bag/rf.cross.validation.")
  positive_class<-ifelse(is.na(positive_class), levels(y)[1], positive_class)
  if(rf_clf_model$error.type=="oob"){
    rf_imp_rank<-rank(-(rf_clf_model$importances)) #rf_all$importances
  }else{
    rf_imp_rank<-rank(-(rowMeans(rf_clf_model$importances))) #rf_all$importances
  }
  max_n<-max(rf_imp_rank, na.rm = TRUE)
  n_total_features<-ncol(x)
  n_features<-c(2, 4, 8, 16, 32, 64, 128, 256, 512, 1024)
  n_features<-n_features[1:which.min(abs(n_features-n_total_features))]
  if(nlevels(y)==2){
    top_n_perf<-matrix(NA, ncol=5, nrow=length(n_features)+1)
    colnames(top_n_perf)<-c("n_features", "AUROC", "Accuracy", "Kappa", "F1")
    rownames(top_n_perf)<-top_n_perf[,1]<-c(n_features, max_n)
    for(i in 1:length(n_features)){
      idx<-which(rf_imp_rank<=n_features[i])
      top_n_features<-names(rf_imp_rank[idx])
      x_n<-x[, top_n_features]
      y_n<-y
      top_n_rf<-rf.out.of.bag(x_n, y_n, ntree=500) # it depends on the rf.out.of.bag defined before
      top_n_AUC<-get.auroc(top_n_rf$probabilities, y_n, positive_class=positive_class)
      top_n_conf<-caret::confusionMatrix(data=top_n_rf$predicted, top_n_rf$y, positive=positive_class)
      top_n_perf[i, 1]<-n_features[i]
      top_n_perf[i, 2]<-top_n_AUC
      top_n_perf[i, 3]<-top_n_conf$overall[1] # Accuracy
      top_n_perf[i, 4]<-top_n_conf$overall[2]# kappa conf$byClass["F1"]
      top_n_perf[i, 5]<-top_n_conf$byClass["F1"]
    }
    all_AUC<-get.auroc(rf_clf_model$probabilities, y_n, positive_class=positive_class)
    all_conf<-caret::confusionMatrix(data=rf_clf_model$predicted, top_n_rf$y, positive=positive_class)
    top_n_perf[length(n_features)+1, ]<-c(max_n, all_AUC, all_conf$overall[1],
                                          all_conf$overall[2], all_conf$byClass["F1"])
  }else{
    top_n_perf<-matrix(NA, ncol=3, nrow=length(n_features)+1)
    colnames(top_n_perf)<-c("n_features", "Accuracy", "Kappa")
    rownames(top_n_perf)<-top_n_perf[,1]<-c(n_features, max_n)
    for(i in 1:length(n_features)){
      idx<-which(rf_imp_rank<=n_features[i])
      top_n_features<-names(rf_imp_rank[idx])
      x_n<-x[, idx]
      y_n<-y
      top_n_rf<-rf.out.of.bag(x_n, y_n, ntree=500) # it depends on the rf.out.of.bag defined before
      top_n_conf<-caret::confusionMatrix(data=top_n_rf$predicted, top_n_rf$y, positive=positive_class)
      top_n_perf[i, 1]<-n_features[i]
      top_n_perf[i, 2]<-top_n_conf$overall[1] # Accuracy
      top_n_perf[i, 3]<-top_n_conf$overall[2]# kappa
    }
    all_conf<-caret::confusionMatrix(data=rf_clf_model$predicted, top_n_rf$y, positive=positive_class)
    top_n_perf[length(n_features)+1, ]<-c(max_n, all_conf$overall[1], all_conf$overall[2])
  }
  top_n_perf<-data.frame(top_n_perf)
  breaks<-top_n_perf$n_features
  p<-ggplot(top_n_perf, aes(x=n_features, y=get(metric))) +
    xlab("# of features used")+
    ylab(metric)+
    scale_x_continuous(trans = "log",breaks=breaks)+
    geom_point() + geom_line()+ ylim(0.5, 1)+
    theme_bw()+
    theme(axis.line = element_line(color="black"),
          axis.title = element_text(size=18),
          strip.background = element_rect(colour = "white"),
          panel.border = element_blank())
  if(!is.null(outdir)){
  ggsave(filename=paste(outdir,"train.",metric,"_VS_top_ranking_features.scatterplot.pdf",sep=""), plot=p, width=5, height=4)
  }
}


