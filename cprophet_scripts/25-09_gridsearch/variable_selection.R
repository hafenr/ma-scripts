createTargetSet <- function(stats, target.var) {
  # Create a training set
  is.target.col <- colnames(stats) == target.var
  is.param.col <- grepl('^param__', colnames(stats))
  trainset <- stats[, is.target.col | is.param.col]
  colnames(trainset) <- sub('param__', '', colnames(trainset))

  trainset
}

fitDecisionTree <- function(stats, target.var) {
  trainset <- createTargetSet(stats, target.var)
  # Decision tree
  fit <- rpart(as.formula(paste(target.var, '~ .')), data=trainset)
  rpart.plot(fit, varlen=0)
  fit
}

fitLASSO <- function(stats, target.var) {
  # LASSO
  trainset.lasso <- createTargetSet(stats, target.var)
  X <- trainset.lasso[, colnames(trainset.lasso) != target.var]
  y <- trainset.lasso[[target.var]]
  is.fac <- sapply(X, is.factor)
  X[is.fac] <- data.frame(do.call(cbind, lapply(X[, is.fac], as.integer)))
  non.changing <- sapply(X, function(col) sd(col) == 0)
  X <- as.matrix(X[, !non.changing])
  X <- scale(X)
  lasso.cv <- cv.glmnet(X, y)
  coef(lasso.cv)
}

fitDecisionTree(stats, 'score_apmw')
