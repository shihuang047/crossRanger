#' @importFrom corrplot corrplot
#' @importFrom grDevices colorRampPalette
#' @importFrom stats cor reorder
#' @importFrom utils write.table
#' @importFrom plyr ldply
#' @importFrom rlang .data
#' @import ggplot2
#' @import viridis
#' @importFrom gridExtra arrangeGrob
#' @importFrom reshape2 melt

#' @title corr_datasets_by_imps
#' @description The correlation of importance scores between datasets
#' @param feature_imps_list a list of feature importance scores from the output of \code{rf_reg.by_datasets} or \code{rf_clf.by_datasets}
#' @param ranked if transform importance scores into rank for correlation analysis
#' @param plot if plot the correlation matrix
#' @seealso ranger rf_clf.by_datasets rf_reg.by_datasets
#' @examples
#'
#' df <- data.frame(rbind(t(rmultinom(14, 14*5, c(.21,.6,.12,.38,.099))),
#'             t(rmultinom(16, 16*5, c(.001,.6,.42,.58,.299))),
#'             t(rmultinom(30, 30*5, c(.011,.6,.22,.28,.289))),
#'             t(rmultinom(30, 30*5, c(.091,.6,.32,.18,.209))),
#'             t(rmultinom(30, 30*5, c(.001,.6,.42,.58,.299)))))
#' metadata<-data.frame(f_s=factor(c(rep("A", 15), rep("B", 15), rep("C", 15), rep("D", 15),
#'                                   rep("A", 15), rep("B", 15), rep("C", 15), rep("D", 15))),
#'                      f_c=factor(c(rep("C", 60), rep("D", 60))),
#'                      f_d=factor(c(rep("A", 30), rep("B", 30), rep("C", 30), rep("D", 30))),
#'                      age=c(1:60, 2:61)
#'                      )
#' reg_res<-rf_reg.by_datasets(df, metadata, s_category='f_d', c_category='age')
#' corr_datasets_by_imps(reg_res$feature_imps_list, plot=TRUE)
#' clf_res<-rf_clf.by_datasets(df, metadata, s_category='f_s', c_category='f_c')
#' feature_imps_list <- lapply(clf_res$feature_imps_list, function(x) x[,"rf_imps"])
#' names(feature_imps_list) <- levels(metadata[, 'f_s'])
#' corr_datasets_by_imps(feature_imps_list, plot=TRUE)
#' @author Shi Huang
#' @export
corr_datasets_by_imps<-function(feature_imps_list, ranked=TRUE, plot=FALSE){
  cor.mtest <- function(mat, ...) {
    mat <- as.matrix(mat)
    n <- ncol(mat)
    p.mat<- matrix(NA, n, n)
    diag(p.mat) <- 0
    for (i in 1:(n - 1)) {
      for (j in (i + 1):n) {
        tmp <- cor.test(mat[, i], mat[, j], ...)
        p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
      }
    }
    colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
    p.mat
  }
  if(ranked){
    ranked_list<-lapply(feature_imps_list, function(x) rank(-x))
    imp_rank_mat <- do.call(cbind, ranked_list)
    corr_mat<-stats::cor(imp_rank_mat)
    # matrix of the p-value of the correlation
    p_mat <- cor.mtest(imp_rank_mat)
  }else{
    imp_mat <- do.call(cbind, feature_imps_list)
    corr_mat<-stats::cor(imp_mat)
    p_mat <- cor.mtest(imp_mat)
  }

  if(plot){
    #require(corrplot)
    col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))
    corrplot(corr = corr_mat, method="color", col=col(200),
             type="upper", addCoef.col = "black",
             tl.col="black", tl.srt=45,
             p.mat = p_mat, sig.level = 0.05, insig = "pch",
             diag = FALSE
             )
  }
  res <-list()
  res$corr_mat <- corr_mat
  res$p_mat <- p_mat
  res
}

#' @title plot_clf_res_list
#' @description Plot a summary of the performance of a list of classification models
#' output from \code{rf_clf.by_datasets}.
#' @param clf_res_list A list, the output of the function \code{rf_clf.by_datasets}.
#' @param q_cutoff A number indicating the cutoff of q values after fdr correction.
#' @param p.adj.method A string indicating the p-value correction method.
#' @param p_cutoff A number indicating the cutoff of p values.
#' @param outdir The output directory.
#' @param plot_height The plot height.
#' @seealso ranger rf_clf.by_datasets
#' @examples
#' df <- data.frame(rbind(t(rmultinom(7, 75, c(.21,.6,.12,.38,.099))),
#'             t(rmultinom(8, 75, c(.001,.6,.42,.58,.299))),
#'             t(rmultinom(15, 75, c(.011,.6,.22,.28,.289))),
#'             t(rmultinom(15, 75, c(.091,.6,.32,.18,.209))),
#'             t(rmultinom(15, 75, c(.001,.6,.42,.58,.299)))))
#' df0 <- data.frame(t(rmultinom(60, 300,c(.001,.6,.2,.3,.299))))
#' metadata<-data.frame(f_s=factor(c(rep("A", 15), rep("B", 15), rep("A", 15), rep("B", 15))),
#'                      f_c=factor(c(rep("C", 30), rep("D", 30))),
#'                      f_d=factor(c(rep("A", 15), rep("B", 15), rep("C", 15), rep("D", 15))),
#'                      age=c(1:30, 2:31)
#'                      )
#' clf_res_list<-rf_clf.by_datasets(df, metadata, s_category="f_c", c_category="f_s")
#' plot_clf_res_list(clf_res_list)
#' @author Shi Huang
#' @export
plot_clf_res_list<-function(clf_res_list, p_cutoff=0.05, p.adj.method = "bonferroni", q_cutoff=0.05, outdir=NULL, plot_height=NULL){
  datasets<-clf_res_list$datasets
  sample_size<-clf_res_list$sample_size
  rf_AUROC<-clf_res_list$rf_AUROC
  rf_AUPRC<-clf_res_list$rf_AUPRC
  feature_imps_list<-clf_res_list$feature_imps_list
  # Statistics summary of biomarkers discovery in multiple datasets
  num_sig_p.adj<-sapply(feature_imps_list, function(x) sum(x[, "non.param.test_p.adj"] < q_cutoff, na.rm=TRUE))
  num_sig_p<-sapply(feature_imps_list, function(x) sum(x[, "non.param.test_p"]< p_cutoff, na.rm=TRUE))
  num_enriched<-sapply(feature_imps_list, function(x) length(grep("enriched", x[,"Enr"])))
  num_depleted<-sapply(feature_imps_list, function(x) length(grep("depleted", x[,"Enr"])))
  # ggplot
  # require('ggplot2')
  summ<-data.frame(datasets=datasets, sample_size=sample_size, AUROC=rf_AUROC, AUPRC=rf_AUPRC, num_sig_p, num_sig_p.adj, num_enriched, num_depleted)
  names(summ)[1:2]<-c("Data_sets", "Sample_size")
  #summ_m<-reshape2::melt(summ)
  p_a<- ggplot(summ, aes(x=.data$Data_sets, y=.data$Sample_size)) + xlab("Data sets") + ylab("Sample size")+
    geom_bar(stat="identity", alpha=0.5, width=0.5)+
    coord_flip()+ # if want to filp coordinate
    theme_bw()
  p_b<- ggplot(summ, aes(x=.data$Data_sets, y=.data$AUROC)) + xlab("") + ylab("AUROC")+
    geom_point(shape="diamond", size=4)+ geom_bar(stat = "identity", alpha=0.5, width=0.01) +
    geom_hline(yintercept=0.5, linetype="dashed")+
    coord_flip()+ # if want to filp coordinate
    theme_bw()+
    scale_y_continuous(limits = c(0.5,1)) +
    theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank())
  p_c<- ggplot(summ, aes(x=.data$Data_sets, y=.data$num_sig_p.adj)) + xlab("") + ylab("# of sig. features")+
    geom_point(size=2)+ geom_bar(stat = "identity", alpha=0.5, width=0.01) +
    coord_flip()+ # if want to filp coordinate
    theme_bw()+
    theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank())
  summ_Enr<-reshape2::melt(summ[, c("Data_sets","num_enriched","num_depleted")])
  summ_Enr[summ_Enr$variable=="num_depleted",]$value<--summ_Enr[summ_Enr$variable=="num_depleted",]$value
  p_d<- ggplot(summ_Enr, aes(x=.data$Data_sets, y=.data$value, colour=.data$variable)) + xlab("") + ylab("Enrichment")+
    ylim(-max(abs(summ_Enr$value)), max(abs(summ_Enr$value)))+
    geom_hline(yintercept=0)+
    geom_point(show.legend=F, size=2)+
    geom_bar(show.legend=F, stat = "identity", alpha=0.5, width=0.01) +
    coord_flip()+ # if want to filp coordinate
    theme_bw()+
    theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank())
  #require('gridExtra')
  p1<-gridExtra::arrangeGrob(p_a, p_b, p_c, p_d, ncol = 4, nrow = 1, widths = c(4, 2, 2, 2))
  ggsave(filename=paste(outdir,"Datasets_AUROC.ggplot.pdf",sep=""),p1,
         width=9, height=ifelse(is.null(plot_height), 3+nrow(summ)*0.2, plot_height))
  # boxplot indicating p and p.adj values of sig. features
  feature_res<-plyr::ldply(feature_imps_list)
  feature_res_m<-reshape2::melt(feature_res[, c("feature","dataset", "Enr", "AUROC", "AUPRC","rf_imps", "non.param.test_p", "non.param.test_p.adj")])
  # customised colors for Enr's 3 factors
  gg_color_hue <- function(n) {
    hues = seq(15, 375, length = n + 1)
    hcl(h = hues, l = 65, c = 100)[1:n]
  }
  my3cols<-c("grey60", rev(gg_color_hue(2)))
  p2<-ggplot(feature_res_m, aes(x=.data$dataset, y=.data$value)) +
    geom_boxplot(outlier.shape = NA)+
    facet_grid(.~.data$variable, scales="free") + scale_color_manual(values = my3cols) +
    coord_flip()+ # if want to filp coordinate
    theme_bw()+
    geom_jitter(aes(color=.data$Enr), position=position_jitter() ,size=1, alpha=0.4)+ #jitter
    theme(axis.line = element_line(color="black"),
          strip.background = element_rect(colour = "white"),
          panel.border = element_blank())
  ggsave(filename=paste(outdir,"Datasets_feature_res.ggplot.pdf",sep=""),plot=p2, width=10, height=4)

  result<-list()
  result$summ<-summ
  result$summ_plot<-p1
  result$feature_res<-feature_res
  result$feature_res_plot<-p2
  return(result)
}

#' @title plot_reg_res_list
#' @description Plot a summary of the performance of a list of classification models
#' output from \code{rf_clf.by_datasets}.
#' @param reg_res_list A list object, the output of the function \code{rf_reg.by_datasets}.
#' @param plot_height The plot height.
#' @param outdir The output directory.
#' @seealso ranger rf_reg.by_datasets
#' @examples
#' df <- data.frame(rbind(t(rmultinom(7, 75, c(.21,.6,.12,.38,.099))),
#'             t(rmultinom(8, 75, c(.001,.6,.42,.58,.299))),
#'             t(rmultinom(15, 75, c(.011,.6,.22,.28,.289))),
#'             t(rmultinom(15, 75, c(.091,.6,.32,.18,.209))),
#'             t(rmultinom(15, 75, c(.001,.6,.42,.58,.299)))))
#' df0 <- data.frame(t(rmultinom(60, 300,c(.001,.6,.2,.3,.299))))
#' metadata<-data.frame(f_s=factor(c(rep("A", 15), rep("B", 15), rep("A", 15), rep("B", 15))),
#'                      f_c=factor(c(rep("C", 30), rep("D", 30))),
#'                      f_d=factor(c(rep("A", 15), rep("B", 15), rep("C", 15), rep("D", 15))),
#'                      age=c(1:30, 2:31)
#'                      )
#' reg_res_list<-rf_reg.by_datasets(df, metadata, s_category="f_c", c_category="age")
#' plot_reg_res_list(reg_res_list)
#' @author Shi Huang
#' @export
plot_reg_res_list<-function(reg_res_list, outdir=NULL, plot_height=NULL){
  datasets<-reg_res_list$datasets
  sample_size<-reg_res_list$sample_size
  rf_MSE<-reg_res_list$rf_MSE
  rf_RMSE<-reg_res_list$rf_RMSE
  rf_MAE<-reg_res_list$rf_MAE
  rf_MAPE<-reg_res_list$rf_MAPE
  rf_Spearman_rho<-reg_res_list$rf_Spearman_rho
  rf_R_squared<-reg_res_list$rf_R_squared
  rf_Adj_R_squared<-reg_res_list$rf_Adj_R_squared
  feature_imps_list<-reg_res_list$feature_imps_list
  # prevelance of features in subdatasets
  prev_df<-do.call(cbind, lapply(reg_res_list$x_list, function(x) apply(x, 2, function(a) sum(a==0)/length(a))))
  # ggplot
  summ<-data.frame(datasets=datasets, sample_size=sample_size,
                   MSE=rf_MSE, RMSE=rf_RMSE, MAE=rf_MAE, MAPE=rf_MAPE, Spearman_rho=rf_Spearman_rho,
                   R_squared=rf_R_squared, Adj_R_squared=rf_Adj_R_squared)
  names(summ)[1:2]<-c("Data_sets", "Sample_size")
  #summ_m<-reshape2::melt(summ)
  p_a<- ggplot(summ, aes(x=.data$Data_sets, y=.data$Sample_size)) + xlab("Data sets") + ylab("Sample size")+
    geom_bar(stat="identity", alpha=0.5, width=0.5)+
    coord_flip()+ # if want to filp coordinate
    theme_bw()
  p_b<- ggplot(summ, aes(x=.data$Data_sets, y=.data$RMSE)) + xlab("") + ylab("RMSE")+
    geom_point(size=2)+ geom_bar(stat = "identity", alpha=0.5, width=0.01) +
    coord_flip()+ # if want to filp coordinate
    theme_bw()+
    theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank())
  p_c<- ggplot(summ, aes(x=.data$Data_sets, y=.data$MAE)) + xlab("") + ylab("MAE")+
    geom_point(size=2)+ geom_bar(stat = "identity", alpha=0.5, width=0.01) +
    coord_flip()+ # if want to filp coordinate
    theme_bw()+
    theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank())
  p_d<- ggplot(summ, aes(x=.data$Data_sets, y=.data$R_squared)) + xlab("") + ylab("R_squared")+
    geom_point(size=2)+ geom_bar(stat = "identity", alpha=0.5, width=0.01) +
    coord_flip()+ # if want to filp coordinate
    theme_bw()+
    theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank())
  #require('gridExtra')
  summ_plot<-arrangeGrob(p_a, p_b, p_c, p_d, ncol = 4, nrow = 1, widths = c(3, 2, 2, 2))
  ggsave(filename=paste(outdir,"Datasets_perfs.ggplot.pdf",sep=""), summ_plot,
         width=9, height=ifelse(is.null(plot_height), 3+nrow(summ)*0.2, plot_height))
  # boxplot indicating p and p.adj values of sig. features
  #feature_res<-plyr::ldply(feature_imps_list)
  feature_res<-do.call(cbind, feature_imps_list)
  result<-list()
  result$summ<-summ
  result$summ_plot<-summ_plot
  result$feature_res<-feature_res
  return(result)
}


#' @title rf_clf.by_datasets.summ
#' @description A integrated pipeline for \code{rf_clf.by_datasets}. It runs standard random forests with oob estimation for classification of
#' c_category in each the sub-datasets splited by the s_category. The output includes a summary of rf models in the sub datasets
#' and all important statistics for each of features.
#' @param df Training data: a data.frame.
#' @param metadata A metadata with at least two categorical variables.
#' @param s_category A string indicates the category in the sample metadata: a ‘factor’ defines the sample grouping for data spliting.
#' @param c_category A indicates the category in the sample metadata: a 'factor' used as sample label for rf classification in each of splited datasets.
#' @param positive_class A string indicates one class in the 'c_category' column of metadata.
#' @param ntree The number of trees.
#' @param nfolds The number of folds in the cross validation.
#' @param p.adj.method The p-value correction method, default is "bonferroni".
#' @param verbose Show computation status and estimated runtime.
#' @param rf_imp_pvalues If compute both importance score and pvalue for each feature.
#' @param p_cutoff The cutoff of p values for features, the default value is 0.05.
#' @param q_cutoff The cutoff of q values for features, the default value is 0.05.
#' @param outdir The output directory, the default is "./".
#' @return A list includes a summary of rf models in the sub datasets,
#' all important statistics for each of features, and plots.
#' @seealso ranger rf_clf.by_datasets
#' @export
rf_clf.by_datasets.summ<-function(df, metadata, s_category, c_category, positive_class=NA,
                                  rf_imp_pvalues=FALSE, nfolds=3, verbose=FALSE, ntree=500,
                                  p_cutoff=0.05, p.adj.method = "bonferroni", q_cutoff=0.05,
                                  outdir=NULL){
  res_list<-rf_clf.by_datasets(df, metadata, s_category, c_category, positive_class,
                               rf_imp_pvalues, nfolds, verbose, ntree)
  expected_list_names<-c("x_list", "y_list","datasets","rf_model_list",
                         "sample_size","rf_AUROC","rf_AUPRC", "feature_imps_list")
  stopifnot(all(names(res_list) %in% expected_list_names==TRUE))
  plot_res_list<-plot_clf_res_list(res_list, p_cutoff=0.05,
                                   p.adj.method = "bonferroni", q_cutoff=0.05,
                                   outdir=outdir)
  result<-list()
  result$rf_models<-res_list$rf_model_list
  result$summ<-plot_res_list$summ
  result$summ_plot<-plot_res_list$summ_plot
  result$feature_res<-plot_res_list$feature_res
  result$feature_res_plot<-plot_res_list$feature_res_plot
  return(result)
}


#' @title rf_clf.comps.summ
#' @description Runs standard random forests with oob estimation for classification of
#' one level VS all other levels of one category in the datasets.
#' The output includes a summary of rf models in the sub datasets,
#' all important statistics for each of features, and plots.
#' @param df Training data: a data.frame.
#' @param f A factor in the metadata with at least two levels (groups).
#' @param comp_group A string indicates the group in the f
#' @param clr_transform A string indicates one class in the 'c_category' column of metadata.
#' @param ntree The number of trees.
#' @param nfolds The number of folds in the cross validation.
#' @param p.adj.method The p-value correction method, default is "bonferroni".
#' @param q_cutoff The cutoff of q values for features, the default value is 0.05.
#' @param p_cutoff The cutoff of p values for features, the default value is 0.05.
#' @param verbose A boolean value indicates if showing computation status and estimated runtime.
#' @param outdir The outputh directory, default is "./".
#' @return A list includes a summary of rf models in the sub datasets,
#' all important statistics for each of features, and plots.
#' @seealso ranger
#' @examples
#' df <- data.frame(t(rmultinom(60, 300,c(.001,.6,.2,.3,.299))))
#' f=factor(c(rep("A", 15), rep("B", 15), rep("C", 15), rep("D", 15)))
#' comp_group="A"
#' rf_clf.comps.summ(df, f, comp_group, verbose=FALSE, ntree=500, p_cutoff=0.05,
#'                   p.adj.method = "bonferroni", q_cutoff=0.05, outdir=NULL)
#'@export
rf_clf.comps.summ<-function(df, f, comp_group, clr_transform=TRUE, nfolds=3, verbose=FALSE, ntree=5000,
                           p_cutoff=0.05, p.adj.method = "bonferroni", q_cutoff=0.05, outdir=NULL){
  clf_comps_list<-rf_clf.comps(df, f, comp_group, verbose=verbose, ntree=ntree, clr_transform=clr_transform)
  expected_list_names<-c("x_list", "y_list", "datasets","rf_model_list","sample_size","rf_AUROC","rf_AUPRC", "feature_imps_list")
  stopifnot(all(names(clf_comps_list) %in% expected_list_names==TRUE))
  plot_res_list<-plot_clf_res_list(clf_res_list=clf_comps_list, p_cutoff=0.05,
                                   p.adj.method = "bonferroni", q_cutoff=0.05, outdir=outdir)
  result<-list()
  result$rf_models<-clf_comps_list$rf_model_list
  result$summ<-plot_res_list$summ
  result$summ_plot<-plot_res_list$summ_plot
  result$feature_res<-plot_res_list$feature_res
  result$feature_res_plot<-plot_res_list$feature_res_plot
  return(result)
}


#' @title plot_cross_appl
#' @description The heatmap indicating ML performance (e.g., MAE) in the self-validation and cross-applications.
#' @param cross_appl_res The output object from \code{rf_clf.cross_appl} or \code{rf_reg.cross_appl}.
#' @param metric The classification performance metric applied.
#' @param plot_height The height (inches) of heatmap.
#' @param plot_width The width (inches) of heatmap
#' @param outdir The outputh directory, default is NULL.
#' @return A heat map showing ML performance in self-validation and cross-applications.
#' @examples
#' df <- data.frame(rbind(t(rmultinom(14, 14*5, c(.21,.6,.12,.38,.099))),
#'             t(rmultinom(16, 16*5, c(.001,.6,.42,.58,.299))),
#'             t(rmultinom(30, 30*5, c(.011,.6,.22,.28,.289))),
#'             t(rmultinom(30, 30*5, c(.091,.6,.32,.18,.209))),
#'             t(rmultinom(30, 30*5, c(.001,.6,.42,.58,.299)))))
#' df0 <- data.frame(t(rmultinom(120, 600,c(.001,.6,.2,.3,.299))))
#' metadata<-data.frame(f_s=factor(rep(c("A", "B"), 60)),
#'                      f_s1=factor(c(rep(TRUE, 60), rep(FALSE, 60))),
#'                      f_c=factor(c(rep("C", 30), rep("H", 30), rep("D", 30), rep("P", 30))),
#'                      age=c(1:60, 2:61)
#'                      )
#'
#' table(metadata[, c('f_s', 'f_c')])
#' clf_res<-rf_clf.by_datasets(df, metadata, nfolds=5, s_category='f_c', c_category='f_s')

#' clf_cross_appl_res <- rf_clf.cross_appl(clf_res$rf_model_list,
#'                                         x_list=clf_res$x_list,
#'                                         y_list=clf_res$y_list)
#' plot_cross_appl(clf_cross_appl_res, metric="AUROC")
#'
#' reg_res<-rf_reg.by_datasets(df, metadata, nfolds=5, s_category='f_c', c_category='age')
#' reg_cross_appl_res <- rf_reg.cross_appl(reg_res,
#'                                         x_list=reg_res$x_list,
#'                                         y_list=reg_res$y_list)
#' plot_cross_appl(reg_cross_appl_res, metric="MAE")
#'@export
plot_cross_appl<-function(cross_appl_res, metric="MAE", outdir=NULL, plot_width=8, plot_height=7){
  if(!class(cross_appl_res)%in%c('rf_reg.cross_appl', 'rf_clf.cross_appl'))
    stop("The class of cross_appl_res should be 'rf_reg.cross_appl' or 'rf_clf.cross_appl'")
  perf_summ<-cross_appl_res$perf_summ
  # The heatmap indicating ML performance (e.g., MAE) in the self-validation and cross-applications
  # require(viridis)
  p<-ggplot(perf_summ, aes(x=as.factor(perf_summ$Test_data), y=as.factor(perf_summ$Train_data), z=get(metric))) +
    xlab("Test data")+ylab("Train data")+
    geom_tile(aes(fill = get(metric), color = perf_summ$Validation_type, width=0.9, height=0.9), size=1) + #
    scale_color_manual(values=c("white","grey80"))+
    geom_text(aes(label = round(get(metric), 2)), color = "white") +
    scale_fill_viridis()+
    labs(fill=metric) +
    theme_bw() + theme_classic() +
    theme(axis.line = element_blank(), axis.text.x = element_text(angle = 90),
          axis.ticks = element_blank())
  if(!is.null(outdir)){
    sink(paste(outdir,"crossRF_reg_perf_summ.xls",sep=""))
    write.table(perf_summ,quote=FALSE,sep="\t", row.names = F)
    sink()
    ggsave(filename=paste(outdir, metric, "_cross_appl_matrix.heatmap.pdf",sep=""),
           plot=p, width=plot_width, height=plot_height)
  }
  p
}

#' @title id_non_spcf_markers
#' @description Non-specific features across datasets (at least present in two of datasets): Last update: 20190130
#' @param feature_res the inheritant output from the function of plot.res_list, rf_clf.comps.summ or rf_clf.by_dataset.summ
#' @param positive_class A string indicates one class in the 'c_category' column of metadata.
#' @param other_class A string indicates the other class in the factor, such as 'health'.
#' @param p.adj.method The p-value correction method, default is "bonferroni".
#' @param outdir The outputh directory, default is "./".
#' @return ...
#' @seealso ranger
#' @examples
#' df <- data.frame(rbind(t(rmultinom(15, 75, c(.21,.6,.12,.38,.099))),
#'             t(rmultinom(15, 75, c(.011,.6,.22,.28,.289))),
#'             t(rmultinom(15, 75, c(.091,.6,.32,.18,.209))),
#'             t(rmultinom(15, 75, c(.001,.6,.42,.58,.299)))))
#' # A dataset with no feature significantly changed in all sub-datasets!
#' df0 <- data.frame(t(rmultinom(60, 300,c(.001,.6,.2,.3,.299))))
#' f=factor(c(rep("A", 15), rep("B", 15), rep("C", 15), rep("D", 15)))
#' comp_group="A"
#' res <-rf_clf.comps.summ(df, f, comp_group, verbose=FALSE, ntree=500,
#'                         p_cutoff=0.05, p.adj.method = "bonferroni",
#'                         q_cutoff=0.05, outdir=NULL)
#' feature_res<-res$feature_res
#' id_non_spcf_markers(feature_res, positive_class="disease",
#'                     other_class="health", p.adj.method= "BH", outdir=NULL)
#' @author Shi Huang
#' @export
id_non_spcf_markers <- function(feature_res, positive_class="disease",
                                other_class="health", p.adj.method = "BH", outdir=NULL){
  feature_res<-feature_res[order(feature_res$dataset, feature_res$feature), ]
  if(all(feature_res$Enr=='Neutral')) stop("No feature significantly changed in any sub-datasets!")
  Enr_dataset_df<-do.call(rbind, tapply(feature_res$Enr, feature_res$feature, summary))
  #Enr_dataset_df<-dcast(feature_res_sig, feature~Enr, fun=length, value.var="Enr")
  shared_global<-as.factor(apply(Enr_dataset_df, 1, function(x) ifelse(any(x[2]==sum(x) | x[3]==sum(x)), "all_shared", "unshared")))
  non_spcf_global<-as.factor(apply(Enr_dataset_df, 1, function(x) ifelse(any(x[2]==sum(x) | x[3]==sum(x)), "all_shared",
                                                                         ifelse(any(x[2]>=2 | x[3]>=2), "non_spcf",
                                                                                ifelse(any(x[3]==0 && x[2]==0), "non_marker", "spcf")))
  ))
  # customised colors for Enr's 3 factors
  gg_color_hue <- function(n) {
    hues = seq(15, 375, length = n + 1)
    hcl(h = hues, l = 65, c = 100)[1:n]
  }
  #---- statistics summary of non_spcf markers aross all datasets
  non_spcf_disease_global<-as.factor(apply(Enr_dataset_df, 1, function(x)
    ifelse(any(x[3]>=2 && x[2]==0), paste("non_spcf_", positive_class, sep=""),
           ifelse(any(x[3]==0 && x[2]>=2), paste("non_spcf_", other_class, sep=""),
                  ifelse(any(x[3]>=1 && x[2]>=1), "non_spcf_mixed",
                         ifelse(any(x[3]==0 && x[2]==0), "non_marker", "spcf"))))
  ))
  #---- calculate the overlap fraction of non-specific markers in each of datasets
  tmp<-data.frame(feature_res, shared_global=as.character(shared_global),
                  non_spcf_global=as.character(non_spcf_global),
                  non_spcf_disease_global=as.character(non_spcf_disease_global),
                  non_spcf_local=as.character(non_spcf_global),
                  non_spcf_disease_local=as.character(non_spcf_disease_global),
                  stringsAsFactors =F)
  tmp[which(tmp[, "Enr"]=="Neutral"),  c("non_spcf_local", "non_spcf_disease_local")]<-"non_marker"
  count_spcf<-do.call(rbind, tapply(tmp$non_spcf_local, tmp$dataset, function(x)
    c(sum(x=="all_shared"),
      sum(x=="non_spcf"),
      sum(x=="spcf"),
      sum(x!="all_shared" & x!="non_spcf" & x!="spcf")
    )
  )
  )
  zero_count <- tapply(tmp$mean_all, tmp$dataset, function(x) sum(x==0))
  count_spcf[, 4]<-count_spcf[, 4]-zero_count
  fac_spcf<-sweep(count_spcf, 1, rowSums(count_spcf), "/")
  colnames(fac_spcf)<-colnames(count_spcf)<-c("All_shared_markers","Non_specific_markers",
                                              "Specific_markers", "Others")
  # remove the columns with total zero values
  fac_spcf<-fac_spcf[, which(!apply(fac_spcf, 2, function(x) all(x==0)))]
  count_spcf<-count_spcf[, which(!apply(count_spcf, 2, function(x) all(x==0)))]
  fac_spcf<-data.frame(dataset=rownames(fac_spcf), fac_spcf)
  # output the fac_spcf table
  sink(paste(outdir,"Markers_fraction_overlap_with_non_specific_",p.adj.method,".txt",sep=""));
  cat("\t"); write.table(fac_spcf,sep='\t',quote=F)
  sink(NULL)
  # output the count_spcf table
  sink(paste(outdir,"Markers_number_overlap_with_non_specific_",p.adj.method,".txt",sep=""));
  cat("\t"); write.table(count_spcf,sep='\t',quote=F)
  sink(NULL)
  # ggplot the fraction overlap of non-specific markers in all datasets
  fac_spcf_m<-reshape2::melt(fac_spcf)
  fac_spcf_m$variable<-factor(fac_spcf_m$variable, levels = rev(levels(fac_spcf_m$variable)))
  mycolors<-c("grey80", gg_color_hue(nlevels(fac_spcf_m$variable)-1))
  p_summ<-ggplot(fac_spcf_m, aes(x=.data$dataset, y=.data$value, fill=.data$variable)) +
    geom_bar(stat="identity")+
    scale_fill_manual(values=mycolors)+
    ylab("Fraction overlap with non-specific markers")+
    coord_flip()+
    #ylim(0, 1)+
    theme_bw()
  ggsave(filename=paste(outdir,"Markers_fraction_overlap_with_non_specific_",p.adj.method,".pdf",sep=""),
         plot=p_summ, width=6, height=4)
  #---- summary of the abundance and occurence rate of
  #---- non-specific disease, non-specific health, non-specific mixed and non markers
  #---- across all patients in all datasets
  feature_res_spcf<-tmp
  p_spcf_abd<-ggplot(feature_res_spcf, aes(x=.data$non_spcf_disease_global, y=log(.data$mean_all)))  +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(position=position_jitter(width = 0.2) ,size=1, alpha=0.4) + #jitter
    xlab("")+ ylab("log10(mean abundance)")+
    coord_flip()+
    theme_bw()
  ggsave(filename=paste(outdir,"Markers_specific_VS_mean_",p.adj.method,".pdf",sep=""),
         plot=p_spcf_abd, width=5, height=3)
  p_spcf_OccRate<-ggplot(feature_res_spcf, aes(x=.data$non_spcf_disease_global, y=.data$OccRate_all))  +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(position=position_jitter(width = 0.2) ,size=1, alpha=0.4) + #jitter
    xlab("")+ ylab("Ubiquity")+
    coord_flip()+
    theme_bw()
  ggsave(filename=paste(outdir,"Markers_specific_VS_ubiquity_",p.adj.method, ".pdf",sep=""),
         plot=p_spcf_OccRate, width=5, height=3)

  result<-list()
  result$fac_spcf<-fac_spcf
  result$count_spcf<-count_spcf
  result$feature_res_spcf<-feature_res_spcf
  result$plot_spcf_summ<-p_summ
  result$plot_spcf_abd<-p_spcf_abd
  result$plot_spcf_OccRate<-p_spcf_OccRate
  return(result)
}

