# Loading libraries
# p <- c("ade4","caret","reshape2","MASS","ranger","plyr",  "foreach", "doMC", "ALEPlot", "vegan", "compositions",
#        "ggplot2","squash","RColorBrewer","pheatmap","pROC","ROCR","caTools")
# usePackage <- function(p) {
#     if (!is.element(p, installed.packages()[,1]))
#         install.packages(p, dep = TRUE, repos = "http://cran.us.r-project.org")
#     suppressWarnings(suppressMessages(invisible(require(p, character.only = TRUE))))
# }
# invisible(lapply(p, usePackage))

#' @title rf.out.of.bag
#' @description Runs standard random forests with out-of-bag error estimation for both classification and regression using \code{ranger}.
#' This is merely a wrapper that extracts relevant info from \code{ranger} output.
#' @importFrom ranger ranger
#' @importFrom Matrix Matrix
#' @importFrom stats model.matrix predict quantile sd
#' @param x Training data: data.matrix or data.frame.
#' @param y Data label for each row.
#' @param ntree The number of trees.
#' @param verbose Show computation status and estimated runtime.
#' @param imp_pvalues If compute both importance score and pvalue for each feature.
#' @return Object of class \code{rf.out.of.bag} with elements including the ranger object and critical metrics for model evaluation.
#'
#' @details
#' Ranger is a fast implementation of random forests (Breiman 2001) or recursive partitioning,
#' particularly suited for high dimensional data.
#' Classification, regression, and survival forests are supported.
#' @references
#' Wright, M. N. & Ziegler, A. (2017). ranger: A fast implementation of random forests for
#' high dimensional data in C++ and R. Journal of Statistical Software 77:1-17.
#' @seealso ranger
#' @examples
#' set.seed(123)
#' x <- data.frame(rbind(t(rmultinom(7, 75, c(.201,.5,.02,.18,.099))),
#'             t(rmultinom(8, 75, c(.201,.4,.12,.18,.099))),
#'             t(rmultinom(15, 75, c(.011,.3,.22,.18,.289))),
#'             t(rmultinom(15, 75, c(.091,.2,.32,.18,.209))),
#'             t(rmultinom(15, 75, c(.001,.1,.42,.18,.299)))))
#' y<-factor(c(rep("A", 15), rep("B", 15), rep("C", 15), rep("D", 15)))
#' rf.out.of.bag(x, y, imp_pvalues=FALSE)
#' rf.out.of.bag(x, y, imp_pvalues=TRUE)
#' x_ <- data.frame(rbind(t(rmultinom(7, 7500, rep(c(.201,.5,.02,.18,.099), 10000))),
#'             t(rmultinom(8, 7500, rep(c(.201,.4,.12,.18,.099), 10000))),
#'             t(rmultinom(15, 7500, rep(c(.011,.3,.22,.18,.289), 10000))),
#'             t(rmultinom(15, 7500, rep(c(.091,.2,.32,.18,.209), 10000))),
#'             t(rmultinom(15, 7500, rep(c(.001,.1,.42,.18,.299), 10000)))))
#' y_<-factor(c(rep("A", 15), rep("B", 15), rep("C", 15), rep("D", 15)))
#' system.time(rf.out.of.bag(x_, y_, imp_pvalues=FALSE))
#' #   user  system elapsed
#' #100.824   0.263  42.947
#' system.time(rf.out.of.bag(x, y, imp_pvalues=TRUE))
#' x_ <- data.frame(rbind(t(rmultinom(7000, 75000, c(.201,.5,.02,.18,.099))),
#'             t(rmultinom(8000, 75000, c(.201,.4,.12,.18,.099))),
#'             t(rmultinom(15000, 75000, c(.011,.3,.22,.18,.289))),
#'             t(rmultinom(15000, 75000, c(.091,.2,.32,.18,.209))),
#'             t(rmultinom(15000, 75000, c(.001,.1,.42,.18,.299)))))
#' y_ <- factor(c(rep("A", 15000), rep("B", 15000), rep("C", 15000), rep("D", 15000)))
#' system.time(rf.out.of.bag(x_, y_, imp_pvalues=FALSE))
#' system.time(rf.out.of.bag(x_, y_, imp_pvalues=TRUE))
#' y0<-factor(c(rep("old", 30), rep("young", 30)))
#' rf.out.of.bag(x, y0, imp_pvalues=FALSE)
#' rf.out.of.bag(x, y0, imp_pvalues=TRUE)
#' y<- 1:60
#' rf.out.of.bag(x, y)
#' rf.out.of.bag(x, y, imp_pvalues=TRUE)
#' @author Shi Huang
#' @export
"rf.out.of.bag" <-function(x, y, ntree=500, verbose=FALSE, sparse = FALSE, imp_pvalues=FALSE, ...){
  set.seed(123)
  data<-data.frame(y, x)
  require(Matrix)
  if(sparse){
    sparse_data <- Matrix::Matrix(data.matrix(data), sparse = sparse)
    rf.model <- ranger::ranger(dependent.variable.name="y", data=sparse_data, keep.inbag=TRUE, importance='permutation',
                               classification=ifelse(is.factor(y), TRUE, FALSE), num.trees=ntree, verbose=verbose, probability = FALSE)
  }else{
    rf.model<-ranger::ranger(y~., data=data, keep.inbag=TRUE, importance='permutation',
                             classification=ifelse(is.factor(y), TRUE, FALSE), num.trees=ntree, verbose=verbose, probability = FALSE)
  }
  result <- list()
  result$rf.model <- rf.model
  result$y <- y
  if(is.factor(y)){
    if(sparse){
      y_numeric<-as.factor(sparse_data[,'y'])
      result$predicted <- factor(rf.model$predictions,levels=levels(y_numeric)); levels(result$predicted)<-levels(y)
    }else{
      result$predicted <- factor(rf.model$predictions,levels=levels(y))
    }
    result$probabilities <- get.oob.probability.from.forest(rf.model, x); colnames(result$probabilities)<-levels(result$y)
    result$confusion.matrix <- t(sapply(levels(y), function(level) table(result$predicted[y==level])))
    result$errs <- mean(result$predicted != result$y)
    #if(nlevels(y)==2) result$auroc <- get.rf.auroc(result$probabilities, y, positive_class = levels(y)[1])
  }else{
    result$predicted <- rf.model$predictions
    reg_perf<-get.reg.oob.performance(rf.model, y)
    result$MSE <- reg_perf$MSE
    result$RMSE <- reg_perf$RMSE
    result$MAE <- reg_perf$MAE
    result$MAE_perc <- reg_perf$MAE_perc
    result$R_squared <- reg_perf$R_squared
    result$Adj_R_squared <- reg_perf$Adj_R_squared
    cat("Mean squared residuals: ", reg_perf$MSE, "\n")
    cat("Mean absolute error: ", reg_perf$MAE, "\n")
    cat("pseudo R-squared (%explained variance): ", reg_perf$R_squared, "\n")
  }
  result$params <- list(ntree=ntree)
  if(imp_pvalues==FALSE){
    result$importances <- rf.model$variable.importance
    result$importances[colSums(x)==0]<-NA
  }else{
    result$importances <- ranger::importance_pvalues(rf.model, method = "altmann", formula = y ~., data=data)
  }
  result$error.type <- "oob"
  class(result) <- "rf.out.of.bag"
  return(result)
}

#' @title balanced.folds
#' @description Get balanced folds where each fold has close to overall class ratio
#' @importFrom stats quantile
#' @param y Data label for each row.
#' @param nfolds The number of folds in the cross validation.
#' @return folds
#' @examples
#' y<-factor(c(rep("A", 15), rep("B", 15), rep("C", 15), rep("D", 15)))
#' print(balanced.folds(y, nfolds=3))
#' y<- 1:60
#' print(balanced.folds(y, nfolds=3))
#' @author Shi Huang
#' @export
balanced.folds <- function(y, nfolds=3){
    folds = rep(0, length(y))
    if(is.factor(y)){
      classes<-levels(y)
    }else{
      y<-factor(findInterval(y, quantile(y, seq(0, 1, by=1/nfolds), type=5), rightmost.closed=TRUE))
      classes<-levels(y)
    }
    # size of each class
    Nk = table(y)
    # -1 or nfolds = len(y) means leave-one-out
    if (nfolds == -1 || nfolds == length(y)){
        invisible(1:length(y))
    }
    else{
    # Can't have more folds than there are items per class
    nfolds = min(nfolds, max(Nk))
    # Assign folds evenly within each class, then shuffle within each class
        for (k in 1:length(classes)){
            ixs <- which(y==classes[k])
            folds_k <- rep(1:nfolds, ceiling(length(ixs) / nfolds))
            folds_k <- folds_k[1:length(ixs)]
            folds_k <- sample(folds_k)
            folds[ixs] = folds_k
        }
        invisible(folds)
    }
}

#' @title rf.cross.validation
#' @description Runs standard random forests with n-folds cross-validation error estimation for both classification and regression using rf.out.of.bag.
#' @param x Training data: data.matrix or data.frame.
#' @param y Data label for each row.
#' @param ntree The number of trees.
#' @param nfolds The number of folds in the cross validation. If nfolds > length(y) or nfolds==-1, uses leave-one-out cross-validation.
#' @param verbose Show computation status and estimated runtime.
#' @param imp_pvalues If compute both importance score and pvalue for each feature.
#' @return Object of class \code{rf.cross.validation} with elements including a \code{ranger} object and mutiple metrics for model evaluations.
#' @seealso ranger
#' @examples
#'
#' x0 <- data.frame(t(rmultinom(16,160,c(.001,.5,.3,.3,.299))) + 0.65)
#' x <- data.frame(rbind(t(rmultinom(7, 75, c(.201,.5,.02,.18,.099))),
#'             t(rmultinom(8, 75, c(.201,.4,.12,.18,.099))),
#'             t(rmultinom(15, 75, c(.011,.3,.22,.18,.289))),
#'             t(rmultinom(15, 75, c(.091,.2,.32,.18,.209))),
#'             t(rmultinom(15, 75, c(.001,.1,.42,.18,.299)))))
#' y<-factor(c(rep("A", 15), rep("B", 15), rep("C", 15), rep("D", 15)))
#' y0<-factor(c(rep("A", 10), rep("B", 30), rep("C", 5), rep("D", 15)))
#' system.time(rf.cross.validation(x, y, imp_pvalues=FALSE))
#' system.time(rf.cross.validation(x, y, imp_pvalues=TRUE))
#' rf.cross.validation(x, y0, imp_pvalues=FALSE)
#' y<- 1:60
#' rf.cross.validation(x, y, nfolds=5, imp_pvalues=FALSE)
#' rf.cross.validation(x, y, nfolds=5, imp_pvalues=TRUE)
#' @author Shi Huang
#' @export
"rf.cross.validation" <- function(x, y, nfolds=3, ntree=500, verbose=FALSE, sparse = FALSE, imp_pvalues=FALSE, ...){
    if(nfolds==-1) nfolds <- length(y)
    folds <- balanced.folds(y, nfolds=nfolds)
    result <- list()
    if(is.factor(y)){
      result$y <- as.factor(y)
      result$errs <- numeric(length(unique(folds)))
      result$probabilities <- matrix(0, nrow=length(result$y), ncol=length(levels(result$y)))
      rownames(result$probabilities) <- rownames(x)
      colnames(result$probabilities) <- levels(result$y)
    }else{
      result$y <- y
      result$MSE <- numeric(length(unique(folds)))
      result$MAE <- numeric(length(unique(folds)))
      result$R_squared <- numeric(length(unique(folds)))
    }
    result$models<-list()
    result$predicted <- result$y
    result$importances <- matrix(0, nrow=ncol(x), ncol=nfolds)
    if(imp_pvalues==TRUE){ result$importance_pvalues <- matrix(0, nrow=ncol(x), ncol=nfolds) }
    # K-fold cross-validation
    for(fold in sort(unique(folds))){
        cat("The # of folds: ", fold, "\n") # if(verbose)
        foldix <- which(folds==fold)
        if(is.factor(y)) y_tr<-factor(result$y[-foldix]) else y_tr<-result$y[-foldix]
        data<-data.frame(y=y_tr, x[-foldix,])
        require(Matrix)
        if(sparse){
          sparse_data <- Matrix::Matrix(data.matrix(data), sparse = TRUE)
          result$models[[fold]]<- model <- ranger::ranger(dependent.variable.name="y", data=sparse_data, classification=ifelse(is.factor(y_tr), TRUE, FALSE),
                                                          keep.inbag=TRUE, importance='permutation', verbose=verbose, num.trees=ntree)
        }else{
          result$models[[fold]]<- model <- ranger::ranger(y~., data=data, keep.inbag=TRUE, importance='permutation',
                                                         classification=ifelse(is.factor(y_tr), TRUE, FALSE), num.trees=ntree, verbose=verbose)
        }

        newx <- x[foldix,]
        if(length(foldix)==1) newx <- matrix(newx, nrow=1)
        predicted_foldix<-predict(model, newx)$predictions
        if(is.factor(y)){
          if(sparse){
            y_numeric<-as.factor(sparse_data[,'y'])
            predicted_foldix <- factor(predicted_foldix, levels=levels(y_numeric))
            levels(predicted_foldix)=levels(y)
          }else{
            predicted_foldix <- factor(predicted_foldix, levels=levels(y))
          }
          result$predicted[foldix] <- predicted_foldix
          probs <- get.predict.probability.from.forest(model, newx); colnames(probs)<-levels(result$y)
          result$probabilities[foldix, colnames(probs)] <- probs
          result$errs[fold] <- mean(result$predicted[foldix] != result$y[foldix])
          result$confusion.matrix <- t(sapply(levels(y), function(level) table(result$predicted[y==level])))
          cat("Error rate: ", result$errs[fold], "\n")
        }else{
          result$predicted[foldix] <- predicted_foldix
          reg_perf<-get.reg.predict.performance(model, newx, newy=result$y[foldix])
          result$MSE[fold] <- reg_perf$MSE
          result$RMSE[fold] <- reg_perf$RMSE
          result$MAE[fold] <- reg_perf$MAE
          result$MAE_perc[fold] <- reg_perf$MAE_perc
          result$R_squared[fold] <- reg_perf$R_squared
          result$Adj_R_squared[fold] <- reg_perf$Adj_R_squared
          cat("Mean squared residuals: ", reg_perf$MSE, "\n")
          cat("Mean absolute error: ", reg_perf$MAE, "\n")
          cat("pseudo R-squared (%explained variance): ", reg_perf$R_squared, "\n")
        }
        if(imp_pvalues==FALSE){
          result$importances[,fold] <- model$variable.importance
          }else{
          imp<-ranger::importance_pvalues(model, method = "altmann", formula = y ~., data=data)
          result$importances[,fold] <- imp[, 1]
          result$importance_pvalues[,fold] <- imp[, 2]
        }

    }
    result$nfolds <- nfolds
    result$params <- list(...)
    result$error.type <- "cv"
    class(result) <- "rf.cross.validation"
    return(result)
}



#' @title get.oob.probability.from.forest
#' @description Get probability of each class using only out-of-bag predictions from RF
#' @param model A object of random forest out-of-bag model.
#' @param x The training data.
#' @return A table of probabilities \code{probs}
#' @author Shi Huang
#' @export
"get.oob.probability.from.forest" <- function(model, x){
    # get aggregated class votes for each sample using only OOB trees
    votes <- get.oob.votes.from.forest(model,x)
    # convert to probs
    probs <- sweep(votes, 1, apply(votes, 1, sum), '/')
    rownames(probs) <- rownames(x)
    colnames(probs) <- model$forest$class.values

    return(invisible(probs))
}

#' @title get.oob.votes.from.forest
#' @description get votes for each class using only out-of-bag predictions from RF
#' @param model A object of class \code{rf.out.of.bag} with elements.
#' @param x The training data.
#' @return votes
#' @author Shi Huang
#' @export
"get.oob.votes.from.forest" <- function(model, x){
    # get aggregated class votes for each sample using only OOB trees
    n_classes <- length(model$forest$class.values)
    votes <- matrix(0, nrow=nrow(x), ncol=n_classes)

    rf.pred <- predict(model, data.frame(x), type="response", predict.all=TRUE)
    model_inbag<-t(do.call(rbind, model$inbag)) # convert the inbag list (length=ntree) into a inbag matrix
    for(i in 1:nrow(x)){
        # find which trees are not inbag for this sample
        outofbag <- model_inbag[i,]==0
        # get oob predictions for this sample
        votes[i,] <- table(factor(rf.pred$predictions[i,][outofbag],levels=seq_along(model$forest$class.values)))
    }
    rownames(votes) <- rownames(x)
    colnames(votes) <- model$forest$levels

    return(invisible(votes))
}

#' @title get.predict.probability.from.forest
#' @description to get votes from the forest in the external test. Last update: Apr. 7, 2019
#' @param model A object of class \code{rf.out.of.bag} with elements.
#' @param newx The external test data.
#' @return \code{probs}
#' @author Shi Huang
#' @export
"get.predict.probability.from.forest" <- function(model, newx){
  # get aggregated class votes for each sample using only OOB trees
  n_classes <- length(model$forest$class.values) #model$forest$levels
  votes <- matrix(0, nrow=nrow(newx), ncol=n_classes)
  rf.pred <- predict(model, data.frame(newx), type="response", predict.all=TRUE)
  for(i in 1:nrow(newx)){
    votes[i,] <- table(factor(rf.pred$predictions[i,], levels=seq_along(model$forest$class.values)))
  }
  rownames(votes) <- rownames(newx)
  colnames(votes) <- seq_along(model$forest$class.values) # model$forest$levels; to correct unanticipated order of class values
  probs <- sweep(votes, 1, apply(votes, 1, sum), '/')
  return(invisible(probs))
}

#' @title get.reg.oob.performance
#' @description Get oob performance from forests for training data. Last update: Apr. 7, 2019
#' @param model A object of class \code{rf.out.of.bag} with elements.
#' @param y The dependent variable in training data.
#' @return \code{perf}
#' @author Shi Huang
#' @export

"get.reg.oob.performance" <- function(model, y){
  pred_y<-as.numeric(model$predictions)
  y<-as.numeric(y)
  MSE <- mean((y-pred_y)^2)
  RMSE <- sqrt(mean((y-pred_y)^2))
  MAE <- mean(sqrt((y-pred_y)^2))
  MAE_perc <- mean(sqrt((y-pred_y)^2)/y)
  R_squared <- 1 - (sum((y-pred_y)^2) / sum((y-mean(y))^2))
  Adj_R_squared <- Adj_R_squared(y, pred_y, k = model$num.independent.variables)
  perf<-list()
  perf$MSE<-MSE
  perf$RMSE<-RMSE
  perf$MAE<-MAE
  perf$MAE_perc<-MAE_perc
  perf$R_squared<-R_squared
  perf$Adj_R_squared<-Adj_R_squared
  return(invisible(perf))
}


#' @title get.reg.predict.performance
#' @description Get prediction performance from forests for the external test data. Last update: Apr. 7, 2019
#' @param model A object of random forest out-of-bag regression model.
#' @param newx The external test data.
#' @param newy The dependent variable in the external test data.
#' @return \code{perf}
#' @author Shi Huang
#' @export
"get.reg.predict.performance" <- function(model, newx, newy){
  rf.pred <- predict(model, data.frame(newx))
  pred_y <- as.numeric(rf.pred$predictions)
  y <- as.numeric(newy)
  MSE <- mean((y-pred_y)^2)
  RMSE <- sqrt(mean((y-pred_y)^2))
  MAE <- mean(sqrt((y-pred_y)^2))
  MAE_perc <- mean(sqrt((y-pred_y)^2)/y)
  R_squared <- 1 - (sum((y-pred_y)^2) / sum((y-mean(y))^2))
  Adj_R_squared <- Adj_R_squared(y, pred_y, k = ncol(newx))
  perf<-list()
  perf$MSE<-MSE
  perf$RMSE<-RMSE
  perf$MAE<-MAE
  perf$MAE_perc<-MAE_perc
  perf$R_squared<-R_squared
  perf$Adj_R_squared<-Adj_R_squared
  return(invisible(perf))
}

#' @title get.auroc
#' @description calculates the area under the ROC curve
#' @param obs A binary factor vector indicates observed values.
#' @param prob A two-column numeric matrix of the same # of rows as the length of observed values,
#'             containing the predicted value of each observation.
#' @importFrom stats model.matrix predict quantile sd
#' @return \code{auroc}
#' @examples
#' obs<-factor(c(rep("A", 31), rep("B", 29)))
#' pred <- c(runif(30, 0.5, 0.9), runif(30, 0, 0.6))
#' prob <-data.frame(A=pred, B=1-pred)
#' get.auroc(obs, prob, positive_class="A")
#' @author Shi Huang
#' @export
get.auroc <- function(prob, obs, positive_class) {
  pos_class<-function(f, positive_class){
    if(nlevels(f)>2) stop("Only two-level factor allowed.")
    idx<-which(levels(f)==positive_class)
    levels(f)[-idx]<-0
    levels(f)[idx]<-1
    factor(f)
  }
  require(ROCR)
  obs<-pos_class(obs, positive_class=positive_class)
  pred <- ROCR::prediction(prob[, positive_class], obs)
  auroc  <- ROCR::performance(pred, "auc")@y.values[[1]]
  auroc
  return(auroc)
}

#' @title get.auprc
#' @description calculates the area under Precision-recall curve
#' @param obs A binary factor vector indicates observed values.
#' @param prob A two-column numeric matrix of the same # of rows as the length of observed values,
#'             containing the predicted value of each observation.
#' @param positive_class A class in the binary factor.
#' @return auprc
#' @examples
#' obs<-factor(c(rep("A", 10), rep("B", 50)))
#' pred <- c(runif(10, 0.4, 0.9), runif(50, 0, 0.6))
#' prob <-data.frame(A=pred, B=1-pred)
#' get.auprc(obs, prob, "A")
#' get.auprc(obs, prob, "B")
#' @author Shi Huang
#' @export
get.auprc <- function(obs, prob, positive_class) {
  pos_class<-function(f, positive_class){
    if(nlevels(f)>2) stop("Only two-level factor allowed.")
    idx<-which(levels(f)==positive_class)
    levels(f)[-idx]<-0
    levels(f)[idx]<-1
    factor(f)
  }
  obs<-pos_class(obs, positive_class=positive_class)
  #require(ROCR)
  #require(caTools)
  xx.df <- ROCR::prediction(prob[, positive_class], obs)
  perf  <- ROCR::performance(xx.df, "prec", "rec")
  xy    <- data.frame(recall=perf@x.values[[1]], precision=perf@y.values[[1]])

  # take out division by 0 for lowest threshold
  xy <- subset(xy, !is.nan(xy$precision))

  # Designate recall = 0 as precision = x...arbitrary
  xy <- rbind(c(0, 0), xy)
  #xy <- xy[!(rowSums(xy)==0), ]

  auprc   <- caTools::trapz(xy$recall, xy$precision)
  auprc
}


#' @title get.mislabel.scores
#' @description Get mislabelled scores of samples
#' @param y Data label for each row.
#' @param y.prob The probability output from the RF model.
#' @return mislabelled scores
#' @examples
#' x <- data.frame(rbind(t(rmultinom(7, 75, c(.201,.5,.02,.18,.099))),
#'             t(rmultinom(8, 75, c(.201,.4,.12,.18,.099))),
#'             t(rmultinom(15, 75, c(.011,.3,.22,.18,.289))),
#'             t(rmultinom(15, 75, c(.091,.2,.32,.18,.209))),
#'             t(rmultinom(15, 75, c(.001,.1,.42,.18,.299)))))
#' y<-factor(c(rep("A", 15), rep("B", 15), rep("C", 15), rep("D", 15)))
#' y<- 1:60
#' rf_model<-rf.out.of.bag(x, y)
#' get.mislabel.scores(rf_model$y, rf_model$probabilities)
#' @author Shi Huang
#' @export
"get.mislabel.scores" <- function(y, y.prob){
  result <- matrix(0,nrow=length(y),ncol=3)
  # get matrices containing only p(other classes), and containing only p(class)
  mm <- stats::model.matrix(~0 + y)
  y.prob.other.max <- apply(y.prob * (1-mm),1,max)
  y.prob.alleged <- apply(y.prob * mm, 1, max)
  result <- cbind(y.prob.alleged, y.prob.other.max, y.prob.alleged - y.prob.other.max)
  rownames(result) <- rownames(y.prob)
  colnames(result) <- c('P(alleged label)','P(second best)','P(alleged label)-P(second best)')
  return(result)
}


#' @title save.rf_clf.results
#' @description Prints random forests results file
#' @importFrom  utils write.table
"save.rf_clf.results" <- function(result, feature.ids, rf.opts){
    save.rf_clf.results.summary(result, rf.opts, outdir=rf.opts$outdir)
    save.rf_clf.results.probabilities(result, outdir=rf.opts$outdir)
    save.rf_clf.results.mislabeling(result, outdir=rf.opts$outdir)
    save.rf_clf.results.importances(result, feature.ids, outdir=rf.opts$outdir)
    save.rf_clf.results.confusion.matrix(result, outdir=rf.opts$outdir)
}

#' @title save.rf_clf.results.summary
#' @description Print "summary" file
#' @importFrom  utils write.table
"save.rf_clf.results.summary" <- function(result, rf.opts, filename='train.summary.xls', outdir='.'){
    filepath <- sprintf('%s/%s',outdir,filename)
    err <- mean(result$errs)
    err.sd <- sd(result$errs)
    y<-result$y
    baseline.err <- 1-max(table(y))/length(y)
    filepath <- sprintf('%s/%s',outdir,filename)
    sink(filepath)
    cat(sprintf('Model\tRandom Forest\n'))
    cat(sprintf('Error type\t%s\n',result$error.type))
    if(result$error.type == 'oob'){
        cat(sprintf('Estimated error\t%.5f\n',err))
    } else {
        cat(sprintf('Estimated error (mean +/- s.d.)\t%.5f +/- %.5f\n',err,err.sd))
    }
    cat(sprintf('Baseline error (for random guessing)\t%.5f\n',baseline.err))
    cat(sprintf('Ratio baseline error to observed error\t%.5f\n',baseline.err / err))
    cat(sprintf('Number of trees\t%d\n',result$params$ntree))
    sink(NULL)
}

#' @title save.rf_clf.results.probabilities
#' @description Print "probabilities" file
#' @importFrom  utils write.table
"save.rf_clf.results.probabilities" <- function(result, filename='train.cv_probabilities.xls', outdir='.'){
    filepath <- sprintf('%s/%s',outdir,filename)
    sink(filepath)
    cat('SampleID\t')
    write.table(result$probabilities,sep='\t',quote=F)
    sink(NULL)
}

#' @title save.rf_clf.results.mislabeling
#' @description Print "mislabeling" file
#' @importFrom  utils write.table
"save.rf_clf.results.mislabeling" <- function(result, filename='train.mislabeling.xls', outdir='.'){
    filepath <- sprintf('%s/%s',outdir,filename)
    sink(filepath)
    cat('SampleID\t')
    write.table(get.mislabel.scores(result$y,result$probabilities),sep='\t',quote=F)
    sink(NULL)
}

#' @title save.rf_clf.results.importances
#' @description Print "feature_importance_scores" file
#' @importFrom  utils write.table
"save.rf_clf.results.importances" <- function(result, feature.ids, filename='train.feature_importance_scores.xls', outdir='.'){
    filepath <- sprintf('%s/%s',outdir,filename)
    if(is.null(dim(result$importances))){
        imp <- result$importances
        imp.sd <- rep(NA,length(imp))
    } else {
        imp <- rowMeans(result$importances)
        imp.sd <- apply(result$importances, 1, sd)
    }
    output.table <- cbind(imp, imp.sd)
    rownames(output.table) <- feature.ids
    output.table <- output.table[sort(imp,dec=T,index=T)$ix,]
    colnames(output.table) <- c('Mean_decrease_in_accuracy','Standard_deviation')

    sink(filepath)
    cat('Feature_id\t')
    write.table(output.table,sep='\t',quote=F)
    sink(NULL)
}

#' @title save.rf_clf.results.confusion.matrix
#' @description Print "confusion matrix" file
#' @importFrom  utils write.table
"save.rf_clf.results.confusion.matrix" <- function(result, filename='train.confusion_matrix.xls', outdir='.'){
    filepath <- sprintf('%s/%s',outdir,filename)
    # add class error column to each row
    x <- result$confusion.matrix
    class.errors <- rowSums(x * (1-diag(nrow(x)))) / rowSums(x)
    output <- cbind(result$confusion.matrix, class.errors)
    colnames(output)[ncol(output)] <- "Class error"
    sink(filepath)
    cat('True\\Predicted\t')
    write.table(output,quote=F,sep='\t')
    sink(NULL)
}
