#' @import foreach
#' @importFrom parallel detectCores
#' @importFrom doMC registerDoMC
#' @importFrom rlang .data
#' @importFrom plyr mapvalues
#' @importFrom stats sd var t.test wilcox.test bartlett.test oneway.test kruskal.test p.adjust


#' @title BetweenGroup.test
#' @description It runs standard univarate statistical tests (such as t.test, wilcox.test,
#' oneway.test, kruskal.test)  for all variates (columns) in the data.frame/data.matrix.
#' @param x A data.martrix or data.frame including multiple numeric vectors.
#' @param y A factor with two or more levels.
#' @param clr_transform A logical value indicates if clr transformation before statistical analysis of compositional microbiome data.
#' @param p.adj.method A string indicating the p-value correction method.
#' @param positive_class A string indicating the specified class in the factor y.
#' @param q_cutoff A number indicating the cutoff of q values after fdr correction.
#' @param paired A logical indicating if paired between-group comparison is desired.
#' @return ...
#' @seealso clr, t.test, wilcox.test, oneway.test, kruskal.test
#' @examples
#' x0 <- data.frame(t(rmultinom(16,160,c(.001,.5,.3,.3,.299))) + 0.65)
#' x <- data.frame(rbind(t(rmultinom(7, 75, c(.201,.5,.02,.18,.099))),
#'             t(rmultinom(8, 75, c(.201,.4,.12,.18,.099))),
#'             t(rmultinom(15, 75, c(.011,.3,.22,.18,.289))),
#'             t(rmultinom(15, 75, c(.091,.2,.32,.18,.209))),
#'             t(rmultinom(15, 75, c(.001,.1,.42,.18,.299)))))
#' # clr_x<-compositions::clr(x)
#' y<-factor(c(rep("A", 30), rep("B", 30)))
#' y1<-factor(c(rep("A", 15), rep("B", 15), rep("C", 15), rep("D", 15)))
#' system.time(BetweenGroup.test(x, y, clr_transform=FALSE))
#' system.time(BetweenGroup.test(x, y, clr_transform=TRUE))
#' system.time(BetweenGroup.test(x, y1, clr_transform=TRUE))
#' x_ <- data.frame(rbind(t(rmultinom(7, 7500, rep(c(.201,.5,.02,.18,.099), 100))),
#'             t(rmultinom(8, 7500, rep(c(.201,.4,.12,.18,.099), 100))),
#'             t(rmultinom(15, 7500, rep(c(.011,.3,.22,.18,.289), 100))),
#'             t(rmultinom(15, 7500, rep(c(.091,.2,.32,.18,.209), 100))),
#'             t(rmultinom(15, 7500, rep(c(.001,.1,.42,.18,.299), 100)))))
#' y_<-factor(c(rep("A", 30), rep("B", 30)))
#' y_1<-factor(c(rep("A", 15), rep("B", 15), rep("C", 15), rep("D", 15)))
#' system.time(BetweenGroup.test(x_, y_))
#' system.time(BetweenGroup.test(x_, y_1))
#' @author Shi Huang
#' @export
BetweenGroup.test <-function(x, y, clr_transform=FALSE, p.adj.method="bonferroni", positive_class=NA, q_cutoff=0.2, paired=FALSE){
  # p.adjust.methods
  # c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none")
  if(clr_transform){ xx<-compositions::clr(x) }else{xx<-x}
  n_group<-nlevels(y)
  positive_class<-ifelse(is.na(positive_class), levels(y)[1], positive_class)
  if(!is.numeric(n_group) | n_group==1)
    stop("group must be a numeric and up to two levels\n")
    test_output<-mttest(xx, y, p.adj.method=p.adj.method, paired=paired)
    #  descriptive statistics of each feature (column) for all samples
    desc_stats_all_df<-desc_stats_all(x, y, clr_transform=clr_transform)
    #  descriptive statistics of each feature (column) for samples (rows) grouped by y
    desc_stats_by_group_df<-desc_stats_by_group(x, y, clr_transform=clr_transform, positive_class=positive_class)
    #-------------------------------Enrichment
    Enr_all<-Enr_by_q_cutoff(test_output, desc_stats_by_group_df, positive_class, q_cutoff=q_cutoff)
    #-------------------------------
    output1<-data.frame(desc_stats_all_df, desc_stats_by_group_df, test_output,  Enr_all)
    output1
}


#' @title log.mat
#' @description Log transformation of a data matrix with pseudo count.
#' @param mat A data.martrix or data.frame including multiple numeric vectors.
#' @param base The base in the log-transformation.
#' @examples
#' x <- data.frame(rbind(t(rmultinom(7, 75, c(.201,.5,.02,.18,.099))),
#'             t(rmultinom(8, 75, c(.201,.4,.12,.18,.099))),
#'             t(rmultinom(15, 75, c(.011,.3,.22,.18,.289))),
#'             t(rmultinom(15, 75, c(.091,.2,.32,.18,.209))),
#'             t(rmultinom(15, 75, c(.001,.1,.42,.18,.299)))))
#' log.mat(x, base=10)
#' @rdname log.mat
#' @author Shi Huang
#' @export
log.mat<-function(mat, base=2){
  if(any(mat == 0)) {
    if(sum(mat == 0) == length(unlist(mat))) {
      mat <- mat + 1e-5
    }
    else{
      v <- as.vector(mat)
      minval <- min(v[v > 0])/2
      mat <- mat + minval
    }
  }
  out<-log(mat, base)
  return(out)
}


#' @title mttest
#' @description Perform the univariate test for all features in the compositional microbiome data.
#' @param x A data.martrix or data.frame including multiple numeric vectors.
#' @param y A factor with two or more levels.
#' @param p.adj.method A string indicating the p-value correction method.
#' @param paired A logical indicating if paired between-group comparison is desired.
#' @examples
#' y <-factor(c(rep("A", 30), rep("B", 30)))
#' y <-factor(c(rep("A", 15), rep("B", 15), rep("C", 15), rep("D", 15)))
#' x <- data.frame(rbind(t(rmultinom(7, 75, c(.201,.5,.02,.18,.099))),
#'             t(rmultinom(8, 75, c(.201,.4,.12,.18,.099))),
#'             t(rmultinom(15, 75, c(.011,.3,.22,.18,.289))),
#'             t(rmultinom(15, 75, c(.091,.2,.32,.18,.209))),
#'             t(rmultinom(15, 75, c(.001,.1,.42,.18,.299)))))
#' mttest(x, y, clr_transform=FALSE)
#' mttest(x, y, clr_transform=TRUE)
#' @author Shi Huang
#' @export
mttest<-function(x, y, p.adj.method="bonferroni", paired=FALSE){
  test_output<-matrix(NA, ncol=4, nrow=ncol(x))
  rownames(test_output)<-colnames(x)
  colnames(test_output)<- c("param.test_p","non.param.test_p","param.test_p.adj","non.param.test_p.adj")
  # comb function for parallelization using foreach
  comb <- function(x, ...) {
    lapply(seq_along(x),
           function(i) c(x[[i]], lapply(list(...), function(y) y[[i]])))
  }
  oper<-foreach::foreach(i=1:ncol(x), .combine='comb', .multicombine=TRUE, .init=list(c(), c())) %dopar% {
    if(stats::var(x[, i])==0){test_out1<-test_out2<-1
    }else{
      if(nlevels(y)==2){
        test_out1<-t.test(x[,i]~y, paired=paired)$p.value
        test_out2<-wilcox.test(x[,i]~y, paired=paired, conf.int=TRUE, exact=FALSE, correct=FALSE)$p.value
      }else{
        test_out1<-oneway.test(x[,i]~y, var.equal=FALSE)$p.value
        test_out2<-kruskal.test(x[,i]~y)$p.value
      }
    }
    out<-c(test_out1, test_out2)
  }
  test_output[,1]<-unlist(oper[[1]])
  test_output[,2]<-unlist(oper[[2]])
  test_output[,3]<-p.adjust(test_output[,1], method = p.adj.method, n = ncol(x))
  test_output[,4]<-p.adjust(test_output[,2], method = p.adj.method, n = ncol(x))
  test_output
}



#' @title desc_stats_all
#' @description The descriptive statistics summary of all features of compositional microbiome data.
#' @param x A data.martrix or data.frame including multiple numeric vectors.
#' @param y A factor with two or more levels.
#' @param clr_transform A logical value indicates if clr transformation before statistical analysis of compositional microbiome data.
#' @param positive_class A string indicating the specified class in the factor y. The positive class should be specified in the case-control design.
#' @examples
#' y <-factor(c(rep("A", 30), rep("B", 30)))
#' x <- data.frame(rbind(t(rmultinom(7, 75, c(.201,.5,.02,.18,.099))),
#'             t(rmultinom(8, 75, c(.201,.4,.12,.18,.099))),
#'             t(rmultinom(15, 75, c(.011,.3,.22,.18,.289))),
#'             t(rmultinom(15, 75, c(.091,.2,.32,.18,.209))),
#'             t(rmultinom(15, 75, c(.001,.1,.42,.18,.299)))))
#' desc_stats_all(x, y, clr_transform=FALSE)
#' desc_stats_all(x, y, clr_transform=TRUE)
#' @author Shi Huang
#' @export
desc_stats_all<-function(x, y, positive_class=NA, clr_transform=FALSE){
  OccRate<-function(x) sum(x!=0)/length(x)
  func_all_list<-c("mean_all", "var_all", "sd_all", "OccRate_all", "AUROC", "AUPRC")
  positive_class<-ifelse(is.na(positive_class), levels(y)[1], positive_class)
  desc_stats_all_df<-data.frame(t(apply(x,2,function(a)
    c(mean(a), stats::var(a), stats::sd(a), OccRate(a), get.auroc(a,y, positive_class), get.auprc(a,y, positive_class)))))
  colnames(desc_stats_all_df)<-func_all_list
  if(clr_transform){
    clr_x<-compositions::clr(x)
    func_all_list<-c("clr_mean_all", "clr_var_all", "clr_sd_all", "clr_AUROC", "clr_AUPRC")
    clr_desc_stats_all_df<-data.frame(t(apply(clr_x, 2, function(a)
      c(mean(a), stats::var(a), stats::sd(a), get.auroc(a,y, positive_class), get.auprc(a,y, positive_class)))))
    colnames(clr_desc_stats_all_df)<-func_all_list
    desc_stats_all_df<-data.frame(desc_stats_all_df, clr_desc_stats_all_df)
  }
  desc_stats_all_df
}


#' @title desc_stats_by_group_using_mean_logfc
#' @description The descriptive statistics summary of all features by grouping of compositional microbiome data.
#' @param x A data.martrix or data.frame including multiple numeric vectors.
#' @param y A factor with two or more levels.
#' @param clr_transform A logical value indicates if clr transformation before statistical analysis of compositional microbiome data.
#' @param positive_class A string indicating the specified class in the factor y. The positive class should be specified in the case-control design.
#' @examples
#' y <-factor(c(rep("A", 30), rep("B", 30)))
#' x <- data.frame(rbind(t(rmultinom(7, 75, c(.201,.5,.02,.18,.099))),
#'             t(rmultinom(8, 75, c(.201,.4,.12,.18,.099))),
#'             t(rmultinom(15, 75, c(.011,.3,.22,.18,.289))),
#'             t(rmultinom(15, 75, c(.091,.2,.32,.18,.209))),
#'             t(rmultinom(15, 75, c(.001,.1,.42,.18,.299)))))
#' desc_stats_by_group(x, y, clr_transform=FALSE)
#' desc_stats_by_group(x, y, clr_transform=TRUE)
#' @author Shi Huang
#' @export
desc_stats_by_group<-function(x, y, clr_transform=FALSE, positive_class=NA){
  positive_class<-ifelse(is.na(positive_class), levels(factor(y))[1], positive_class)
  OccRate<-function(x) sum(x!=0)/length(x)
  log10_median<-function(x, base=10) log.mat(stats::median(x), base=10)
  log10_mean<-function(x, base=10) log.mat(stats::median(x), base=10)

  mean_logfc<-function(x, y, base=2, positive_class.=positive_class){
    if(nlevels(y)>2) levels(y)[levels(y)!=positive_class] <- "Others"
    logMeanAbd<-log.mat(t(apply(x,2,function(x) tapply(x, y, mean))), base=base)
    out<-logMeanAbd[, positive_class]-logMeanAbd[, colnames(logMeanAbd)!=positive_class]
    out
  }

  median_logfc <- function(x, y, base=2, positive_class.=positive_class){
    if(nlevels(y)>2) levels(y)[levels(y)!=positive_class] <- "Others"
    logMedianAbd<-log.mat(t(apply(x,2,function(x) tapply(x, y, median))), base=base)
    out<-logMedianAbd[, positive_class]-logMedianAbd[, colnames(logMedianAbd)!=positive_class]
    out
  }

  new_quantile <- function(x) {
    quantile(x, probs = seq(0, 1, 0.1))
  }
  generalized_logfc <- function(x, y, base=2, positive_class.=positive_class){
    if(nlevels(y)>2) levels(y)[levels(y)!=positive_class] <- "Others"
    QuantileAbd <- apply(x,2,function(x) tapply(x, y, new_quantile))
    results <- rep(0, ncol(x))
    names(results) <- colnames(x)
    for(i in 1:length(QuantileAbd)) {
      QuantileAbd[[i]] <- do.call(data.frame, QuantileAbd[[i]])
      logQuantileAbd <- log.mat(QuantileAbd[[i]], base = base)
      out<-logQuantileAbd[, positive_class]-logQuantileAbd[, colnames(logQuantileAbd)!=positive_class]
      out<-mean(out)
      results[i] <- out
    }

    results
  }

  func_by_group_list<-c("mean", "sd", "median", "log10_median", "OccRate")
  tmp<-apply(x,2,function(x) tapply(x, y, function(x) c(mean(x), stats::sd(x), stats::median(x), log10_median(x), OccRate(x))));
  desc_stats_by_group_df<-t(sapply(tmp, unlist));
  colnames(desc_stats_by_group_df)<-unlist(lapply(levels(y), function(x) paste(func_by_group_list, x, sep="__")))
  desc_stats_by_group_df<-data.frame(desc_stats_by_group_df,
                                     mean_logfc=mean_logfc(x, y),
                                     median_logfc = median_logfc(x, y),
                                     generalized_logfc = generalized_logfc(x, y))
  if(clr_transform){
    clr_x<-compositions::clr(x)
    func_by_group_list<-c("clr_mean", "clr_sd", "clr_median")
    tmp<-apply(clr_x, 2, function(x) tapply(x, y, function(x) c(mean(x), stats::sd(x), stats::median(x))));
    clr_desc_stats_by_group_df<-t(sapply(tmp, unlist));
    colnames(clr_desc_stats_by_group_df)<-unlist(lapply(levels(y), function(x) paste(func_by_group_list, x, sep="__")))
    desc_stats_by_group_df<-data.frame(desc_stats_by_group_df, clr_desc_stats_by_group_df)
  }
  desc_stats_by_group_df
}

#' @title Enr_by_q_cutoff
#' @description The significance of all features by q values.
#' @param test_output The output from function \code{mttest}.
#' @param desc_stats_by_group_df The output from function \code{desc_stats_by_group}.
#' @param positive_class A string indicating the specified class in the factor y.
#' @param q_cutoff A number indicating the cutoff of q values after fdr correction.
#' @examples
#' y <-factor(c(rep("A", 30), rep("B", 30)))
#' y1<-factor(c(rep("A", 15), rep("B", 15), rep("C", 15), rep("D", 15)))
#' x <- data.frame(rbind(t(rmultinom(7, 75, c(.201,.5,.02,.18,.099))),
#'             t(rmultinom(8, 75, c(.201,.4,.12,.18,.099))),
#'             t(rmultinom(15, 75, c(.011,.3,.22,.18,.289))),
#'             t(rmultinom(15, 75, c(.091,.2,.32,.18,.209))),
#'             t(rmultinom(15, 75, c(.001,.1,.42,.18,.299)))))
#' test_output<-mttest(x, y, p.adj.method="bonferroni", paired=FALSE)
#' test_output<-mttest(x, y1, p.adj.method="bonferroni", paired=FALSE)
#' desc_stats_by_group_df<-desc_stats_by_group(x, y, clr_transform=FALSE)
#' desc_stats_by_group_df<-desc_stats_by_group(x, y1, clr_transform=FALSE)
#' Enr_by_q_cutoff(test_output, desc_stats_by_group_df, q_cutoff=0.05, "A")
#' @export
Enr_by_q_cutoff<-function(test_output, desc_stats_by_group_df, positive_class=NA, q_cutoff=0.05){
  # positive_class<-ifelse(is.na(positive_class), levels(factor(y))[1], positive_class)
  # require("plyr")
  IfSig<-as.factor(ifelse(test_output[, "non.param.test_p.adj"]< q_cutoff, "Sig", "NotSig"))
  Enr0<-factor(ifelse(desc_stats_by_group_df$mean_logfc>0, paste(positive_class, "enriched", sep="_"), paste(positive_class, "depleted", sep="_")))
  IfSigEnr<-interaction(IfSig, Enr0)
  Enr<-plyr::mapvalues(IfSigEnr,
                       c(paste("NotSig.",positive_class,"_depleted",sep=""),
                         paste("NotSig.",positive_class,"_enriched",sep=""),
                         paste("Sig.",positive_class,"_depleted",sep=""),
                         paste("Sig.",positive_class,"_enriched",sep="")),
                       c("Neutral", "Neutral",
                         paste(positive_class,"_depleted",sep=""),
                         paste(positive_class,"_enriched",sep="")))
  data.frame(IfSig,IfSigEnr,Enr)
}


