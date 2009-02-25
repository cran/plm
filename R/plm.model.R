plm.within <- function(formula, data, effect){
  pdim <- pdim(data)

  X <- model.matrix(formula, data, part = "first", model = "within", effect = effect)
  y <- pmodel.response(data, part = "first", model = "within", effect = effect)
  if (length(formula) == 2){
    W <- model.matrix(formula, data, part = "second", model = "within", effect = effect)
    if (ncol(W) < ncol(X)) stop("insufficient number of instruments")
  }
  else W <- NULL
  result <- mylm(y, X, W)
  card.fixef <- switch(effect,
                       "individual" = pdim$nT$n,
                       "time"       = pdim$nT$T,
                       "twoways"    = pdim$nT$n + pdim$nT$T - 1
                       )
  df <- df.residual(result) - card.fixef
  result <- list(coefficients = coef(result),
                 vcov         = result$vcov*df.residual(result)/df,
                 residuals    = resid(result),
                 df.residual  = df,
                 formula      = formula,
                 model        = data)
  class(result) <- c("plm", "panelmodel")
  result
}    
    
plm.between <- function(formula, data, effect){
  X <- model.matrix(formula, data, part = "first", model = "between", effect = effect)
  y <- pmodel.response(data, model="between", effect = effect)
  names.X <- colnames(X)
  if (length(formula) == 2){
    W <- model.matrix(formula, data, part = "second",model = "between", effect = effect)
    if (ncol(W) < ncol(X)) stop("insufficient number of instruments")
  }
  else W <-  NULL
  result <- mylm(y, X, W)
  result <- list(coefficients = coef(result),
                 vcov         = result$vcov,
                 residuals    = resid(result),
                 df.residual  = df.residual(result),
                 formula      = formula,
                 model        = data)
  class(result) <- c("plm", "panelmodel")
  result
}    

plm.pooling <- function(formula, data){
  X <- model.matrix(formula, data, part = "first", model = "pooling")
  y <- pmodel.response(data, model = "pooling")
  if (length(formula) == 2){
    W <- model.matrix(formula, data, part = "second", model = "pooling")
    if (ncol(W) < ncol(X)) stop("insufficient number of instruments")
  }
  else W <- NULL
  result <- mylm(y, X, W)
  result <- list(coefficients = coef(result),
                 vcov         = result$vcov,
                 residuals    = resid(result),
                 df.residual  = df.residual(result),
                 formula      = formula,
                 model        = data)
  class(result) <- c("plm","panelmodel")
  result  
}

plm.fd <- function(formula, data){
  X <- model.matrix(formula, data, part = "first", model = "fd");
  y <- pmodel.response(data, model = "fd")
  N <- nrow(X)
  if (length(formula) == 2){
    W <- model.matrix(formula, data, part = "second", model = "fd")
    if (ncol(W) < ncol(X)) stop("insufficient number of instruments")
  }
  else W <- NULL
  result <- mylm(y, X, W)
  result <- list(coefficients = coef(result),
                 vcov         = result$vcov,
                 residuals    = resid(result),
                 df.residual  = df.residual(result),
                 formula      = formula,
                 model        = data)

  class(result) <- c("plm","panelmodel")
  result  
}

plm.random <- function(formula, data, effect, ercomp, ivar){
  has.instruments <- length(formula) == 2
  pdim <- pdim(data)
  balanced <- pdim$balanced
  estec <- ercomp(formula, data, effect, method = ercomp)
  sigma2 <- estec$sigma2
  theta <- estec$theta
  if (effect == "individual") cond <- data[["(id)"]]
  if (effect == "time") cond <- data[["(time)"]]
  if (effect == "twoways"){
    if (has.instruments) stop("Instrumental variable random effect
               estimation not implemented for two-ways panels")
    X <- model.matrix(formula, data, model = "random", effect = effect, theta = theta)
    y <- pmodel.response(data, model = "random", effect = effect, theta = theta) 
    result <- mylm(y, X, NULL)
  }
  else{
    if (length(theta) > 1) theta <- theta[as.character(cond)]
    X <- model.matrix(formula, data, model = "random", effect = effect, theta = theta)
    y <- pmodel.response(data, model = "random", effect = effect, theta = theta)
    if (!has.instruments){
      result <- mylm(y, X, NULL)
    }
    else{
      X <- X/sqrt(sigma2$idios)
      y <- y/sqrt(sigma2$idios)
      Wb <- model.matrix(formula, data, part = "second", model = "Between", effect = effect)
      Ww <- model.matrix(formula, data, part = "second", model = "within", effect = effect)
      if (ncol(Wb) < ncol(X)) stop("Insufficient number of instruments\n")
      if(ivar == "baltagi") W <- cbind(Ww,Wb)
      if (! balanced) sig2one <- sigma2$one[as.character(cond)] else sig2one <- sigma2$one
      if(ivar == "bvk"){
        # to construct the W matrix, we must take into account the fact that 
        # ncol(Wb) >= ncol(Ws)
        W <- Wb/sqrt(sig2one)
        W[,colnames(Ww)] <- W[,colnames(Ww)]+Ww/sqrt(sigma2$idios)
      }
      if (ncol(W) < ncol(X)) stop("insufficient number of instruments")
      result <- mylm(y, X, W)
    }
  }
  result <- list(coefficients = coef(result),
                 vcov         = result$vcov,
                 residuals    = resid(result),
                 df.residual  = df.residual(result),
                 formula      = formula, 
                 model        = data,
                 ercomp       = estec)
  class(result) <- c("plm","panelmodel")
  result
}

plm.ht <- function(formula, data){
  # estimate the within model without instrument

  old.formula <- formula
  formula <- pFormula(formula(formula))
  within <- plm.within(formula, data, effect = "individual")
  fixef <- fixef(within)

  formula <- old.formula
  id <- data[["(id)"]]
  time <- data[["(time)"]]
  pdim <- pdim(data)
  balanced <- pdim$balanced
  T <- pdim$nT$T
  n <- pdim$nT$n
  N <- pdim$nT$N
  Ti <- pdim$Tint$Ti

  # get the typology of the variables
  X <- model.matrix(formula, data, model = "within")
  W <- model.matrix(formula, data, model = "within", part = "second")

  exo.all <- colnames(W)
  tot.all <- colnames(X)
  tot.cst <- attr(X,"timeconst")
  tot.var <- tot.all[!(tot.cst %in% tot.all)]

  exo.cst <- attr(W,"timeconst")
  exo.var <- exo.all[!(exo.all %in% exo.cst)]
  end.cst <- tot.cst[!(tot.cst %in% exo.cst)]
  end.var <- tot.var[!(tot.var %in% exo.var)]
  if (length(end.cst) > length(exo.var)){
    stop(" The number of endogenous time-invariant variables is greater
           than the number of exogenous time varying variables\n")
    }
  
  X <- model.matrix(formula, data, model = "pooling")
  if (length(exo.var) > 0) XV <- X[ , exo.var, drop = FALSE] else XV <- NULL
  if (length(end.var) > 0) NV <- X[ , end.var, drop = FALSE] else NV <- NULL
  if (length(exo.cst) > 0) XC <- X[ , exo.cst, drop = FALSE] else XC <- NULL
  if (length(end.cst) > 0) NC <- X[ , end.cst, drop = FALSE] else NC <- NULL

  sigma2 <- list()
  sigma2$one <- 0
  sigma2$idios <- deviance(within)/(N-n)
  if (length(tot.cst) !=0 ){
    zo <- twosls(fixef[as.character(id)],cbind(XC,NC),cbind(XC,XV))
  }
  else{
    zo <- lm(fixef~1)
  }
  
  ssr <- deviance(zo)/N

  if(balanced){
    sigma2$id <- ssr-sigma2$idios/T
    theta <- 1-sqrt(sigma2$idios/(sigma2$idios+T*sigma2$id))
  }
  else{
    sigma2$id <- ssr-sigma2$idios/T
    theta <- 1-sqrt(sigma2$idios/(sigma2$idios+Ti*sigma2$id))
    theta <- theta[as.character(id)]
  }

  estec <- structure(list(sigma2 = sigma2, theta = theta),
                     class = "ercomp",
                     balanced = balanced,
                     effect = "individual")
  y <- pmodel.response(data, model = "random", theta = theta)
  X <- model.matrix(formula, data, model = "random", theta = theta)
  within.inst <- model.matrix(formula, data, model = "within")
  between.inst <- model.matrix(formula, data, model = "Between",
                               part = "second")[, exo.var, drop = FALSE] 
  W <- cbind(within.inst, XC, between.inst)
  result <- twosls(y,X,W)
  K <- length(data)
  ve <- lev2var(data[,-c(K-1,K)])
  varlist <- list(xv = ve[exo.var],
                  nv = ve[end.var],
                  xc = ve[exo.cst[exo.cst != "(Intercept)"]],
                  nc = ve[end.cst]
                  )
  varlist <- lapply(varlist, function(x){ names(x) <- NULL; x})

  result <- list(coefficients = coef(result),
                 vcov         = vcov(result),
                 residuals    = resid(result),
                 df.residual  = df.residual(result),
                 formula      = formula, 
                 model        = data,
                 varlist      = varlist,
                 ercomp      = estec)
  names(result$coefficients) <- colnames(result$vcov) <-
    rownames(result$vcov) <- colnames(X)
  class(result) <- c("plm","panelmodel")
  result
}

mylm <- function(y, X, W = NULL){
  names.X <- colnames(X)
  if (is.null(W)) result <- lm(y ~ X - 1)
  else result <- twosls(y, X, W)
  if (any(is.na(coef(result)))){
    na.coef <- is.na(coef(result))
    X <- X[, !na.coef, drop = FALSE]
    if (is.null(W)) result <- lm(y ~ X - 1)
    else result <- twosls(y, X, W)
  }
  result$vcov <- vcov(result)
  names(result$coefficients) <- colnames(result$vcov) <-
    rownames(result$vcov) <- colnames(X)
  result
}
