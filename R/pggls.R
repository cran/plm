pggls <- function(formula,data,effect="individual",model,...){
  cl <- match.call()
  model.name <- model
  if (!any(class(data) %in% "pdata.frame")){
    stop("argument data should be a pdata.frame\n")
  }
  if(!(effect %in% names(effect.pggls.list))){
    stop(paste("effect must be one of",oneof(effect.pggls.list)))
  }
  if(!(model.name %in% names(model.pggls.list))){
    stop("model must be one of",oneof(model.pggls.list))
  }

  require(kinship)
  x <- plm.extract(formula=formula,data=data)
  pdim <- attr(x,"pdim")
  pvar <- attr(x,"pvar")
  balanced <- pdim$balanced
  time.names <- pdim$panel.names$time.names
  id.names <- pdim$panel.names$id.names

  pmodel <- list(model.name=model,formula=formula,effect=effect)
  nt <- pdim$Tint$nt
  Ti <- pdim$Tint$Ti
  T <- pdim$nT$T
  n <- pdim$nT$n
  N <- pdim$nT$N
  X <- x$X
  y <- x$y
  id <- x$id
  time <- x$time

  model.res <- switch(model,"random"="pooling",
                      "within"="within")
  m <- plm(formula=formula,data=data,effect=effect,model=model.res)

  coef.names <- colnames(X)
  K <- pdim$K <- ncol(X)
  
  if (effect=="time"){
    cond <- time
    other <- id
    cond.variation <- pvar$time.variation # is this needed?
    other.variation <- pvar$id.variation  # idem
    ncond <- T
    nother <- n
    cond.names <- time.names
    other.names <- id.names
    groupsdim <- nt
  } else {
    cond <- id
    other <- time
    other.variation <- pvar$id.variation  # idem
    cond.variation <- pvar$time.variation # idem
    ncond <- n
    nother <- T
    cond.names <- id.names
    other.names <- time.names
    groupsdim <- Ti
  }
  
  myord <- order(cond,other)
  
  ## reorder data
  resid <- m$residuals[myord]
  X <- m$model$X[myord,]
  y <- m$model[[1]][myord]
  
  ## conditions hav to be reordered as well
  cond<-cond[myord]
  other<-other[myord]
  
  ## drop first time period (see Wooldridge 10.5, eq. 10.61)
  drop1<-FALSE
  ## the function turns out to work irrespective of dropping
  ## one time period or not!! absolutely the same results...
  ## The 'if' parameterization is just for debugging. Set drop1=T
  ## if you want to avoid any problems...
  if(drop1 && model=="within") {
    numeric.t <- as.numeric(other)
    t1 <- which(numeric.t!=min(numeric.t))
    resid <- resid[t1]
    X0 <- X
    y0 <- y
    X <- X[t1,]
    y <- y[t1]
    cond<-cond[t1]
    other<-other[t1]
    nother<-nother-1
    other.names<-other.names[-1]
  }
  
  tres <- array(NA,dim=c(nother,nother,ncond),dimnames=list(other.names,other.names,cond.names))
  lcnd <- levels(cond)
  if(pdim$balanced){
    for (i in 1:ncond){
      ut <- resid[cond==lcnd[i]]
      tres[,,i] <- ut%o%ut
    }
    subOmega <- apply(tres,1:2,mean)
    omega <- bdsmatrix(rep(nother,ncond),rep(subOmega,ncond))
  } else {
    lti <- list()
    for (i in 1:ncond){
      cond.i <- cond==lcnd[i]
      ut <- resid[cond.i]
      names(ut) <- lti[[i]] <- other[cond.i]
      out <- ut%o%ut
      tres[names(ut),names(ut),i] <- out
    }
    subOmega <- apply(tres,1:2,mean,na.rm=TRUE)
    list.cov.blocks <- list()
    for (i in 1:ncond){
      list.cov.blocks[[i]] <- subOmega[lti[[i]],lti[[i]]]
    }
    omega <- bdsmatrix(groupsdim,unlist(list.cov.blocks))
  }
  A <- crossprod(X,solve(omega,X))
  B <- crossprod(X,solve(omega,y))
  vcov <- solve(A) # was: solve(crossprod(X,solve(omega,X)))
  coef <- as.vector(solve(A,B))
  if(drop1 && model=="within"){
    X <- X0
    y <- y0
  }
  residuals <- y-as.vector(crossprod(t(X),coef))
  coef.names <- switch(model,"within"=coef.names,
                       "random"=c("(intercept)",coef.names)
                               )
  df.residual <- nrow(X)-ncol(X)
  fitted.values <- y-residuals
  model <- data.frame(y,X)
  names(model)[[1]] <- deparse(formula[[2]])
  names(coef) <- rownames(vcov) <- colnames(vcov) <- coef.names
  fullGLS <- list(coefficients=coef,residuals=residuals,fitted.values=fitted.values,
                  vcov=vcov,df.residual=df.residual,model=m$model,sigma=subOmega,call=cl)
  fullGLS <- structure(fullGLS,pdim=pdim,pvar=pvar,pmodel=pmodel)
  class(fullGLS)=c("pggls","panelmodel")
  fullGLS

}

summary.pggls <- function(object,...){
  pmodel <- attr(object,"pmodel")
  std.err <- sqrt(diag(object$vcov))
  b <- object$coefficients
  z <- b/std.err
  p <- 2*(1-pnorm(abs(z)))
  CoefTable <- cbind(b,std.err,z,p)
  colnames(CoefTable) <- c("Estimate","Std. Error","z-value","Pr(>|z|)")
  object$CoefTable <- CoefTable
  y <- object$model[,1]
  object$tss <- tss(y)
  object$ssr <- sum(residuals(object)^2)
  object$rsqr <- 1-object$ssr/object$tss
  class(object) <- c("summary.pggls")
  return(object)
}

print.summary.pggls <- function(x,digits=5,length.line=70,...){
  pmodel <- attr(x,"pmodel")
  pdim <- attr(x,"pdim")
  effect <- pmodel$effect
  formula <- pmodel$formula
  model.name <- pmodel$model.name
  centre("Model Description",length.line)
  cat(paste(effect.pggls.list[[effect]],"\n",sep=""))
  cat(paste(model.pggls.list[[model.name]],"\n",sep=""))
  print.form(formula,"Model Formula             : ",length.line)
  centre("Panel Dimensions",length.line)
  print(pdim)
  centre("Residuals",length.line)
  save.digits <- unlist(options(digits=digits))
  on.exit(options(digits=save.digits))
  print(summary(unlist(residuals(x))))
  centre("Coefficients",length.line)
  printCoefmat(x$CoefTable,digits=digits)
  centre("Overall Statistics",length.line)
  cat(paste("Total Sum of Squares       : ",signif(x$tss,digits),"\n",sep=""))
  cat(paste("Residual Sum of Squares    : ",signif(x$ssr,digits),"\n",sep=""))
  cat(paste("Rsq                        : ",signif(x$rsqr,digits),"\n",sep=""))
  cat(paste(trait(length.line),"\n"))
  invisible(x)
}
