pvcm.extract <- function(formula,data,effect){
  pvar <- attr(data,"pvar")
  pdim <- attr(data,"pdim")
  indexes <- attr(data,"indexes")
  time.names <- pdim$panel.names$time.names
  id.names <- pdim$panel.names$id.names
  id.index.name <- indexes$id
  time.index.name <- indexes$time
  id <- data[[id.index.name]]
  time <- data[[time.index.name]]
  name.data <- paste(deparse(substitute(data)))
  formod <- paste(deparse(formula),collapse="")
  fortot <- as.formula(paste(formod,"+",id.index.name,"+",time.index.name))

  cond <- id
  other <- time
  if (effect=="time"){
    cond <- time
    other <- id
  }

  mf <- model.frame(fortot,data=data)
  nr <- nrow(mf)
  if(nr!=nrow(data)){
    for (i in names(mf)){
      if (is.factor(mf[[i]])){
        x <- mf[[i]]
        if (length(unique(x)) < length(levels(x))){
          mf[[i]] <- mf[[i]][,drop=TRUE]
        }
      }
    }
  }
  id <- mf[[id.index.name]]
  time <- mf[[time.index.name]]
  mtx <- terms(formula)
  X <- model.matrix(mtx,mf)[,-1,drop=F]
  y <- model.response(mf)
  yX <- cbind(y,X)
  name.y <- paste(deparse(formula[[2]]))
  colnames(yX)[1] <- name.y
  pdim <- pdim(id,time)
  pvar <- pvar(X,id,time)
  ml <- split(as.data.frame(yX),cond)
  ols <- lapply(ml,lm)
  extr <- list(ml=ml,ols=ols,id=id,time=time)
  extr <- structure(extr,pdim=pdim,pvar=pvar,indexes=indexes,data=name.data)
  extr
}

pvcm <-  function(formula,data,effect="individual",model,...){
  model.name <- model
  cl <- match.call()

  if (!any(class(data) %in% "pdata.frame")){
    stop("argument data should be a pdata.frame\n")
  }
  if(!(effect %in% names(effect.pvcm.list))){
    stop(paste("effect must be one of",oneof(effect.pvcm.list)))
  }

  if(!(model.name %in% names(model.pvcm.list))){
    stop("model must be one of",oneof(model.pvcm.list))
  }
  
  pmodel <- list(model.name=model.name,formula=formula,effect=effect)
  
  extr <- pvcm.extract(formula,data,effect)
  extr <- structure(extr,pmodel=pmodel)
  
  result <- switch(model.name,
                   "within"=pvcm.within(extr,...),
                   "random"=pvcm.random(extr,...)
                   )
  class(result) <- c("pvcm","panelmodel")
  result$call <- cl
  result
}

pvcm.within <- function(extr,...){
  ml <- extr$ml
  ols <- extr$ols
  pdim <- attr(extr,"pdim")
  pmodel <- attr(extr,"pmodel")
  T <- pdim$nT$T
  n <- pdim$nT$n
  N <- pdim$nT$N
  ncond <- n

  if (pmodel$effect=="time"){
    ncond <- T
  }
  coef <- as.data.frame(t(sapply(ols,coefficients)))
  residuals <- lapply(ols,residuals)
  residuals <- unlist(residuals)
  vcov <- lapply(ols,vcov)
  std <- as.data.frame(t(sapply(vcov,function(x) sqrt(diag(x)))))
  names(coef) <- names(std) <- colnames(coef)
  ssr <- sum(residuals^2)
  y <- unlist(lapply(ml,function(x) x[,1]))
  fitted.values <- y-residuals
  tss <- tss(y)
  df.residuals <- N-ncond*ncol(coef)
  nopool <- list(coefficients=coef,residuals=residuals,fitted.values=fitted.values,
                 vcov=vcov,df.residuals=df.residuals,model=ml,std.error=std)
  nopool <- structure(nopool,class="pvcm",pmodel=pmodel,pdim=pdim,pvar=pvar)
  nopool
}

pvcm.random <- function(extr,...){
  ols <- extr$ols
  ml <- extr$ml
  pdim <- attr(extr,"pdim")
  pmodel <- attr(extr,"pmodel")
  T <- pdim$nT$T
  n <- pdim$nT$n
  N <- pdim$nT$N
  ncond <- n

  if (pmodel$effect=="time"){
    ncond <- T
  }
  coefm <- t(sapply(ols,coefficients))
  coef <- lapply(ols,coefficients)
  K <- ncol(coefm)-1
  res <- lapply(ols,residuals)
  coefb <- apply(coefm,2,sum)/ncond
  xpxm1 <- lapply(ml,function(x) solve(crossprod(as.matrix(cbind(1,x[,-1])))))
  D1 <- suml(lapply(oppl(coef,coefb,"-"),function(x) crossprod(t(x))))/(ncond-1)
  sigi <- lapply(res,function(x) sum(x^2)/(length(x)-K-1))
  sigim <- lapply(sigi,function(x) x*matrix(1,K+1,K+1))
  s2xpxm1i <- mapply("*",xpxm1,sigim,SIMPLIFY=F)
  D2 <- matrix(apply(mapply("*",xpxm1,sigi),1,sum)/ncond,ncol=K+1)
  eig <- prod(eigen(D1-D2)$values>=0)
  if (eig){
    Delta <- D1-D2
  }
  else{
    Delta <- D1 # si le Delta precedent nes pas semi defini positif
  }
  Delta <- list(Delta)
  Wn <- mapply("+",s2xpxm1i,rep(Delta,ncond),SIMPLIFY=F)
  Wnm <- sapply(Wn,solve)
  Wn <- lapply(Wn,solve)
  vcovb <- solve(matrix(apply(Wnm,1,sum),K+1,K+1))
  W <- lapply(Wn,function(x) vcovb%*%x)
  beta <- apply(mapply("%*%",W,coef),1,sum)
  df.residuals <- N-ncol(coefm)
  X <- lapply(ml,function(x) as.matrix(cbind(1,x[,-1])))
  y <- lapply(ml,function(x) x[,1])
  tss <- tss(unlist(y))
  haty <- oppl(X,beta,"%*%")
  residuals <- mapply("-",y,haty,SIMPLIFY=FALSE)
  ssr <- sum(unlist(residuals)^2)
  names.res <- lapply(residuals,rownames)
  residuals <- lapply(residuals,as.vector)
  zo <- function(x,y){
    names(x) <- y
    x
  }
  residuals <- mapply(zo,residuals,names.res,SIMPLIFY=FALSE)
  residuals <- unlist(residuals)
  y <- unlist(residuals)
  fitted.values <- y-residuals
  names(beta) <- rownames(vcovb) <- colnames(vcovb) <- colnames(coefm)
  swamy <- list(coefficients=beta,residuals=residuals,fitted.values=fitted.values,
                 vcov=vcovb,df.residuals=df.residuals,model=ml,Delta=Delta[[1]])
  swamy <- structure(swamy,pdim=pdim,pvar=pvar,pmodel=pmodel)
  swamy
}

summary.pvcm <- function(object,...){
  pmodel <- attr(object,"pmodel")
  if (pmodel$model.name=="random"){
    std.err <- sqrt(diag(object$vcov))
    b <- object$coefficients
    z <- b/std.err
    p <- 2*(1-pnorm(abs(z)))
    CoefTable <- cbind(b,std.err,z,p)
    colnames(CoefTable) <- c("Estimate","Std. Error","z-value","Pr(>|z|)")
    object$CoefTable <- CoefTable
  }
  object$ssr <- sum(residuals(object)^2)
  y <-  unlist(lapply(object$model,function(x) x[,1]))
  object$tss <- tss(y)
  object$rsqr <- 1-object$ssr/object$tss
  class(object) <- c("summary.pvcm")
  return(object)
}

print.summary.pvcm <- function(x,digits=5,length.line=70,...){
  pmodel <- attr(x,"pmodel")
  pdim <- attr(x,"pdim")
  effect <- pmodel$effect
  formula <- pmodel$formula
  model.name <- pmodel$model.name
  centre("Model Description",length.line)
  cat(paste(effect.pvcm.list[[effect]],"\n",sep=""))
  cat(paste(model.pvcm.list[[model.name]],"\n",sep=""))
  print.form(formula,"Model Formula             : ",length.line)
  centre("Panel Dimensions",length.line)
  print(pdim)
  centre("Residuals",length.line)
  save.digits <- unlist(options(digits=digits))
  on.exit(options(digits=save.digits))
  print(summary(unlist(residuals(x))))
  if (model.name=="random"){
  centre("Estimated mean of the coefficients",length.line)
    printCoefmat(x$CoefTable,digits=digits)
    centre("Estimated variance of the coefficients",length.line)
    print(x$Delta,digits=digits)
  }
  if (model.name=="within"){
    centre("Coefficients",length.line)
    print(summary(x$coefficients))
  }
  centre("Overall Statistics",length.line)
  cat(paste("Total Sum of Squares       : ",signif(x$tss,digits),"\n",sep=""))
  cat(paste("Sum of Squares Residuals   : ",signif(x$ssr,digits),"\n",sep=""))
  cat(paste("Rsq                        : ",signif(x$rsqr,digits),"\n",sep=""))
  cat(paste(trait(length.line),"\n"))
  invisible(x)
}
