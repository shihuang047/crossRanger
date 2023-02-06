#' @importFrom stats model.matrix predict quantile sd cor.test
#' @importFrom Matrix Matrix
#' @importFrom ranger ranger importance_pvalues
#' @importFrom PRROC pr.curve roc.curve
#' @importFrom ROCR prediction performance

#' @title rf.out.of.bag
#' @description It runs standard random forests with out-of-bag error estimation for both classification and regression using \code{ranger}.
#' This is merely a wrapper that extracts relevant info from \code{ranger} output.
#' @param x Training data: data.matrix or data.frame.
#' @param y A response vector. If a factor, classification is assumed, otherwise regression is assumed.
#' @param ntree The number of trees.
#' @param sparse A boolean value indicates if the input matrix transformed into sparse matrix for rf modeling.
#' @param verbose A boolean value indicates if showing computation status and estimated runtime.
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
#' y<-factor(c(rep("A", 20), rep("B", 20), rep("C", 20)))
#' rf.out.of.bag(x, y, imp_pvalues=FALSE)
#' rf.out.of.bag(x, y, imp_pvalues=TRUE)
#' x_ <- data.frame(rbind(t(rmultinom(7, 7500, rep(c(.201,.5,.02,.18,.099), 1000))),
#'             t(rmultinom(8, 750, rep(c(.201,.4,.12,.18,.099), 1000))),
#'             t(rmultinom(15, 750, rep(c(.011,.3,.22,.18,.289), 1000))),
#'             t(rmultinom(15, 750, rep(c(.091,.2,.32,.18,.209), 1000))),
#'             t(rmultinom(15, 750, rep(c(.001,.1,.42,.18,.299), 1000)))))
#' y_<-factor(c(rep("A", 15), rep("B", 15), rep("C", 15), rep("D", 15)))
#' rf.out.of.bag(x_, y_, imp_pvalues=FALSE)
#' rf.out.of.bag(x, y_, imp_pvalues=TRUE)
#' y0<-factor(c(rep("old", 30), rep("young", 30)))
#' rf.out.of.bag(x, y0, imp_pvalues=FALSE)
#' rf.out.of.bag(x, y0, imp_pvalues=TRUE)
#' y<- 1:60
#' rf.out.of.bag(x, y, imp_pvalues=FALSE)
#' @author Shi Huang
#' @export
"rf.out.of.bag" <-function(x, y, ntree=500, verbose=FALSE, sparse = FALSE, imp_pvalues=FALSE){
  set.seed(123)
  data<-data.frame(y, x)
  #require(Matrix)
  if(sparse){
    sparse_data <- Matrix::Matrix(data.matrix(data), sparse = sparse)
    rf.model <- ranger(dependent.variable.name="y", data=sparse_data, keep.inbag=TRUE, importance='permutation',
                               classification=ifelse(is.factor(y), TRUE, FALSE), num.trees=ntree, verbose=verbose, probability = FALSE)
  }else{
    rf.model<-ranger(dependent.variable.name="y", data=data, keep.inbag=TRUE, importance='permutation',
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
    result<-append(result, reg_perf)
    cat("Mean squared residuals: ", reg_perf$MSE, "\n")
    cat("Mean absolute error: ", reg_perf$MAE, "\n")
    cat("pseudo R-squared (%explained variance): ", reg_perf$R_squared, "\n")
  }
  if(imp_pvalues==FALSE){
    result$importances <- rf.model$variable.importance
    result$importances[colSums(x)==0]<-NA
  }else{
      result$importances <- importance_pvalues(rf.model, method = "altmann", formula = y ~., data=data) # if sparse_data, there will NOT pvalue output
  }
  result$params <- list(ntree=ntree)
  result$error.type <- "oob"
  class(result) <- "rf.out.of.bag"
  return(result)
}

#' @title balanced.folds
#' @description Get balanced folds where each fold has close to overall class ratio
#' @param y A response vector. If a factor, classification is assumed, otherwise regression is assumed.
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
    }else{
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
#' @description It runs standard random forests with n-folds cross-validation error estimation for both classification and regression using rf.out.of.bag.
#' @param x Training data: data.matrix or data.frame.
#' @param y A response vector. If a factor, classification is assumed, otherwise regression is assumed.
#' @param ntree The number of trees.
#' @param nfolds The number of folds in the cross validation. If nfolds > length(y)
#' or nfolds==-1, uses leave-one-out cross-validation. If nfolds was a factor, it means customized folds (e.g., leave-one-group-out cv) were set for CV.
#' @param sparse A boolean value indicates if the input matrix transformed into sparse matrix for rf modeling.
#' @param verbose A boolean value indicates if showing computation status and estimated runtime.
#' @param imp_pvalues If compute both importance score and pvalue for each feature.
#' @param ... Other parameters applicable to `ranger`.
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
#' s<-factor(c(rep("A", 15), rep("B", 15), rep("C", 15), rep("D", 15)))
#' y<-factor(rep(c("Y", "N"), 30))
#' y0<-factor(c(rep("A", 10), rep("B", 30), rep("C", 5), rep("D", 15)))
#' system.time(rf.cross.validation(x, y, imp_pvalues=FALSE))
#' system.time(rf.cross.validation(x, y, imp_pvalues=TRUE))
#' rf.cross.validation(x, y0, imp_pvalues=FALSE)
#' y_n<- 1:60
#' rf.cross.validation(x, y_n, nfolds=5, imp_pvalues=FALSE)
#' rf.cross.validation(x, y_n, nfolds=5, imp_pvalues=TRUE)
#' # when nfolds is a factor, it actually run a leave-one-group-out cv
#' rf.cross.validation(x, y, nfolds=s, imp_pvalues=TRUE)
#' rf.cross.validation(x, y_n, nfolds=s, imp_pvalues=FALSE)
#' @author Shi Huang
#' @export
"rf.cross.validation" <- function(x, y, nfolds=3, ntree=500, verbose=FALSE, sparse = FALSE, imp_pvalues=FALSE, ...){
  if(length(y)!=dim(x)[1]) stop("The target varible doesn't match the shape of train data matrix, please check!")
  if(is.factor(nfolds)){
    folds = factor(nfolds) # it means the customized folds were set
    nfolds = nlevels(nfolds) # how many folds in this customized vector
  }else{
    if(nfolds==-1) nfolds <- length(y)
    folds <- balanced.folds(y, nfolds=nfolds)
  }
  if(any(table(folds)<3)) stop("Less than 3 samples in at least one fold!\n")
  result <- list()
  result$rf.model<-list()
  if(is.factor(y)){
    result$y <- y <- factor(y)
    result$predicted <- y
    result$probabilities <- matrix(0, nrow=length(y), ncol=length(levels(y)))
    rownames(result$probabilities) <- rownames(x)
    colnames(result$probabilities) <- levels(y)
    result$confusion.matrix<-t(sapply(levels(y), function(level) table(y[y==level])))
    result$errs <- numeric(length(unique(folds))); names(result$errs) = sort(unique(folds))
  }else{
    result$y <- y
    result$predicted <- y
  }

  if(imp_pvalues==TRUE){ result$importance_pvalues <- matrix(0, nrow=ncol(x), ncol=nfolds) }
  # K-fold cross-validation
  for(fold in sort(unique(folds))){
    cat("The # of folds: ", fold, "\n") # if(verbose)
    foldix <- which(folds==fold)
    if(is.factor(y)) y_tr<-factor(result$y[-foldix]) else y_tr<-result$y[-foldix]
    data<-data.frame(y=y_tr, x[-foldix,], check.names = FALSE)
    if(sparse){
      sparse_data <- Matrix::Matrix(data.matrix(data), sparse = TRUE)
      result$rf.model[[fold]]<- model <- ranger(dependent.variable.name="y", data=sparse_data, classification=ifelse(is.factor(y_tr), TRUE, FALSE),
                                                keep.inbag=TRUE, importance='permutation', verbose=verbose, num.trees=ntree, ...)
    }else{
      result$rf.model[[fold]]<- model <- ranger(dependent.variable.name="y", data=data, keep.inbag=TRUE, importance='permutation',
                                                classification=ifelse(is.factor(y_tr), TRUE, FALSE), num.trees=ntree, verbose=verbose, ...)
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
      result$nRMSE[fold] <- reg_perf$nRMSE
      result$MAE[fold] <- reg_perf$MAE
      result$MAPE[fold] <- reg_perf$MAPE
      result$MASE[fold] <- reg_perf$MASE
      result$Spearman_rho[fold] <- reg_perf$Spearman_rho
      result$R_squared[fold] <- reg_perf$R_squared
      result$Adj_R_squared[fold] <- reg_perf$Adj_R_squared

      cat("Mean squared residuals: ", reg_perf$MSE, "\n")
      cat("Mean absolute error: ", reg_perf$MAE, "\n")
      cat("pseudo R-squared (%explained variance): ", reg_perf$R_squared, "\n")
    }
  }
  result$importances <- matrix(0, nrow=ncol(x), ncol=nfolds)
  colnames(result$importances) <- sort(unique(folds))
  rownames(result$importances) <- colnames(x)
  if(imp_pvalues) result$importance_pvalues <- result$importances
  for(fold in sort(unique(folds))){
    if(imp_pvalues==FALSE){
      result$importances[,fold] <- result$rf.model[[fold]]$variable.importance
    }else{
      imp<-importance_pvalues(result$rf.model[[fold]], method = "altmann", formula = y ~., data=data)
      result$importances[,fold] <- imp[, 1]
      result$importance_pvalues[,fold] <- imp[, 2]
    }
  }
  result$params <- list(ntree=ntree, nfolds=nfolds)
  result$error.type <- "cv"
  result$nfolds <- nfolds
  class(result) <- "rf.cross.validation"
  return(result)
}


#' @title get.oob.probability.from.forest
#' @description Get probability of each class using only out-of-bag predictions from RF
#' @param model A model object in the class of `rf.out.of.bag`.
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
  rf.pred <- predict(model, data.frame(newx, check.names = FALSE), type="response", predict.all=TRUE)
  for(i in 1:nrow(newx)){
    votes[i,] <- table(factor(rf.pred$predictions[i,], levels=seq_along(model$forest$class.values)))
  }
  rownames(votes) <- rownames(newx)
  colnames(votes) <- seq_along(model$forest$class.values) # model$forest$levels; to correct unanticipated order of class values
  probs <- sweep(votes, 1, apply(votes, 1, sum), '/')
  return(invisible(probs))
}

#' @title get.reg.performance
#' @description Get regression performance
#' @param pred_y The predicted y in the numeric format.
#' @param y A response vector for regression.
#' @param n_features The number of features.
#' @return A list of regression performance metrics.
#' @author Shi Huang
#' @export
"get.reg.performance" <- function(pred_y, y, n_features=NA){
  pred_y<-as.numeric(pred_y)
  y<-as.numeric(y)
  diff_y<-diff(range(y))
  mean_y<-mean(y)
  MSE <- mean((y-pred_y)^2)
  RMSE <- sqrt(mean((y-pred_y)^2))
  nRMSE <- RMSE/mean_y # <https://en.wikipedia.org/wiki/Root-mean-square_deviation#Normalized_root-mean-square_deviation>
  MAE <- mean(sqrt((y-pred_y)^2))
  MAPE <- mean(sqrt((y-pred_y)^2)/y) # <https://en.wikipedia.org/wiki/Mean_absolute_percentage_error>
  MASE <- mean(sqrt((y-pred_y)^2))/mean(sqrt((diff(y))^2)) # <https://en.wikipedia.org/wiki/Mean_absolute_scaled_error>
  Spearman_rho <- cor.test(y, pred_y, method = "spearman")$estimate
  R2 <- function(y, pred_y){1 - (sum((y-pred_y)^2) / sum((y-mean(y))^2))}
  R_squared <- R2(y, pred_y)
  adj.R2<-function(y, pred_y, k){
    n=length(y);
    1-(1-R2(y, pred_y)^2)*(n-1)/(n-k-1) # k is # of predictors
  }
  if(!is.na(n_features)){
    Adj_R_squared <- adj.R2(y, pred_y, k = n_features)
  }else{
    Adj_R_squared <- NA
  }
  perf<-list()
  perf$MSE<-MSE
  perf$RMSE<-RMSE
  perf$nRMSE<-nRMSE
  perf$MAE<-MAE
  perf$MAPE<-MAPE
  perf$MASE<-MASE
  perf$Spearman_rho <- Spearman_rho
  perf$R_squared<-R_squared
  perf$Adj_R_squared<-Adj_R_squared
  return(invisible(perf))
}

#' @title get.reg.oob.performance
#' @description Get oob performance from forests for training data. Last update: Apr. 7, 2019
#' @param model A object of class \code{rf.out.of.bag} with elements.
#' @param y The responsive variable for regression in training data.
#' @return \code{perf}
#' @author Shi Huang
#' @export
"get.reg.oob.performance" <- function(model, y){
  pred_y<-as.numeric(model$predictions)
  y<-as.numeric(y)
  perf<-get.reg.performance(pred_y, y, n_features=model$num.independent.variables)
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
  perf<-get.reg.performance(pred_y, y, n_features=ncol(newx))
  return(invisible(perf))
}

#' @title get.auroc
#' @description calculates the area under the ROC curve
#' @param y A binary factor vector indicates observed values.
#' @param predictor A predictor. For example, one column in the probability table, indicating
#'            the probability of a given observation being the positive class in the factor y.
#' @param positive_class A class of the factor y.
#' @return The auroc value.
#' @examples
#' y<-factor(c(rep("A", 31), rep("B", 29)))
#' y1<-factor(c(rep("A", 18), rep("B", 20), rep("C", 22)))
#' y0<-factor(rep("A", 60))
#' pred <- c(runif(30, 0.5, 0.9), runif(30, 0, 0.6))
#' prob <-data.frame(A=pred, B=1-pred)
#' positive_class="A"
#' get.auroc(predictor=prob[, positive_class], y, positive_class="A")
#' get.auroc(predictor=prob[, positive_class], y, positive_class="B")
#' get.auroc(predictor=prob[, positive_class], y0, positive_class="A")
#' get.auroc(predictor=prob[, positive_class], y1, positive_class="A")
#' @author Shi Huang
#' @export
get.auroc <- function(predictor, y, positive_class) {
  pos_class<-function(f, positive_class){
    #if(nlevels(f)>2) cat("More than two levels in the factor.")
    idx<-which(levels(f)==positive_class)
    levels(f)[-idx]<-0
    levels(f)[idx]<-1
    factor(f)
  }
  #require(ROCR)
  if(nlevels(factor(y))==1){
    cat("All y values are identical!")
    auroc <- NA
  }else{
    y<-pos_class(y, positive_class=positive_class)
    pred <- prediction(predictor, y)
    auroc  <- performance(pred, "auc")@y.values[[1]]
  }
  return(auroc)
}

#' @title get.auprc
#' @description calculates the area under Precision-recall curve (AUPRC).
#' @param y A binary factor vector indicates observed values.
#' @param predictor A predictor. For example, one column in the probability table, indicating
#'            the probability of a given observation being the positive class in the factor y.
#' @param positive_class A class of the factor y.
#' @details It’s a bit trickier to interpret AUPRC than it is to interpret AUROC.
#' The AUPRC of a random classifier is equal to the fraction of positives (Saito et al.),
#' where the fraction of positives is calculated as (# positive examples / total # examples).
#' That means that _different_ classes have _different_ AUPRC baselines. A class with 12% positives
#' has a baseline AUPRC of 0.12, so obtaining an AUPRC of 0.40 on this class is great. However
#' a class with 98% positives has a baseline AUPRC of 0.98, which means that obtaining an AUPRC
#' of 0.40 on this class is bad.
#' @return auprc
#' @examples
#' y<-factor(c(rep("A", 10), rep("B", 50)))
#' y1<-factor(c(rep("A", 18), rep("B", 20), rep("C", 22)))
#' y0<-factor(rep("A", 60))
#' pred <- c(runif(10, 0.4, 0.9), runif(50, 0, 0.6))
#' prob <-data.frame(A=pred, B=1-pred)
#' positive_class="A"
#' get.auprc(predictor=prob[, positive_class], y, positive_class="A")
#' get.auprc(predictor=prob[, positive_class], y, positive_class="B")
#' get.auprc(predictor=prob[, positive_class], y0, positive_class="A")
#' get.auprc(predictor=prob[, positive_class], y1, positive_class="A")
#' @author Shi Huang
#' @export
get.auprc<- function(predictor, y, positive_class){
  if(nlevels(factor(y))==1){
    cat("All y values are identical!")
    auprc <- NA # if y contains all identical values, no threshold can be used for auprc calculation
  }else{
    df<-data.frame(y, predictor)
    prob_pos<-df[df$y==positive_class, "predictor"]
    prob_neg<-df[df$y!=positive_class,"predictor"]
    pr<-pr.curve(prob_pos, prob_neg, curve=FALSE,  max.compute = T, min.compute = T, rand.compute = T)
    #rel_auprc<-round((pr$auc.integral-pr$min$auc.integral)/(pr$max$auc.integral-pr$min$auc.integral), 3)
    auprc<-pr$auc.integral
    if(pr$auc.integral==pr$rand$auc.integral) cat("Note: ", auprc, "is the auprc of a random classifier!\n")
  }
  auprc
}


#' @title get.mislabel.scores
#' @description Get mislabelled scores of samples
#' @param y A response vector. If a factor, classification is assumed, otherwise regression is assumed.
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
  mm <- model.matrix(~0 + y)
  y.prob.other.max <- apply(y.prob * (1-mm),1,max)
  y.prob.alleged <- apply(y.prob * mm, 1, max)
  result <- cbind(y.prob.alleged, y.prob.other.max, y.prob.alleged - y.prob.other.max)
  rownames(result) <- rownames(y.prob)
  colnames(result) <- c('P(alleged label)','P(second best)','P(alleged label)-P(second best)')
  return(result)
}
