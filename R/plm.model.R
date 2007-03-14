plm.within <- function(y,X,W,id,time,pvar,pdim,pmodel,indexes,cl,...){
  formula <- pmodel$formula
  effect <- pmodel$effect
  rnames <- rownames(X)
  T <- pdim$nT$T
  n <- pdim$nT$n
  N <- pdim$nT$N

  K <- pdim$K <- ncol(X)
  coef.names <- colnames(X)

  cond <- id
  other <- time
  other.variation <- pvar$time.variation
  cond.variation <- pvar$id.variation
  ncond <- n

  if (effect=="time"){
    cond <- time
    other <- id
    cond.variation <- pvar$time.variation
    other.variation <- pvar$id.variation
    ncond <- T
  }

  Kw <- pdim$Kw <- sum(other.variation)
  Kd <- pdim$Kd <- sum(cond.variation & other.variation)

  X.m <- papply(X,mymean,cond)
  X.with <- (X-X.m)[,other.variation,drop=F]
  y.m <- papply(unclass(y),mymean,cond)
  y.with <- y-y.m

  if(effect=="twoways"){
    X.time <- papply(X,mymean,time)
    X.mean <- matrix(rep(apply(X,2,mean),N),ncol=K,byrow=T)
    X.with <- (X-X.m-X.time+X.mean)[,cond.variation & other.variation,drop=F]
    y.time <- papply(unclass(y),mymean,time)
    y.with <- y-y.m-y.time+mean(y)
    df.within <- N-T-n-Kd+1
    coef.within <- other.variation & cond.variation
  }
  else{
    coef.within <- other.variation
    df.within <- N-n-Kw
  }
  if (nrow(X.with)<=ncol(X.with)){
    stop("Within estimation impossible (insufficient number of observations)\n")
  }
  
  if(is.null(W)){
    within <- lm(y.with~X.with-1)
    if (effect=="twoways"){
      FE.time <- y.time-as.vector(X.time[,cond.variation & other.variation,drop=F]%*%within$coef)
    }
  }
  else{
    if (ncol(W)<ncol(X)+1) stop("Insufficient number of instruments\n")
    W.m <- W.m <- papply(W,mymean,cond)
    if(effect=="twoways"){
      W.time <- papply(W,mymean,time)
      W.mean <- matrix(rep(apply(W,2,mean),N),ncol=K,byrow=T)
      W.with <- (W-W.m-W.time+W.mean)[,cond.variation & other.variation]
    }
    else{
      W.with <- W-W.m
    }
    within <- twosls(y.with,X.with,W.with)
  }
  
  within <- plmformat(within,coef.names[coef.within],rnames,
                      df.within,"within",pdim,pmodel,indexes,cl)

  within$FE <- y.m-as.vector(X.m[,coef.within,drop=F]%*%within$coef)
  attr(within$FE,"cm") <- tapply(within$FE,cond,mean)
  within$alpha <- as.vector(mean(y)-apply(X[,coef.within,drop=F],2,mean)%*%within$coef)
  if(effect=="twoways"){
    within$FE.time <- FE.time
  }
  within
}    

plm.between <- function(y,X,W,id,time,pvar,pdim,pmodel,indexes,cl,...){
  if (pmodel$effect=="twoways"){
    effect <- "individual"
  }
  else{
    effect <- pmodel$effect
  }
    
  formula <- pmodel$formula
  effect <- pmodel$effect

  balanced <- pdim$balanced

  T <- pdim$nT$T
  n <- pdim$nT$n
  N <- pdim$nT$N

  K <- pdim$K <- ncol(X)
  coef.names <- colnames(X)

  cond <- id
  other <- time
  other.variation <- pvar$time.variation
  cond.variation <- pvar$id.variation
  ncond <- n

  if (effect=="time"){
    cond <- time
    other <- id
    cond.variation <- pvar$time.variation
    other.variation <- pvar$id.variation
    ncond <- T
  }

  Kb <- pdim$Kb <- sum(cond.variation)

  X.m <- papply(X,mymean,cond)
  X.bet <- attr(X.m,"cm")[,cond.variation,drop=F]
  y.m <- papply(unclass(y),mymean,cond)
  y.bet <- attr(y.m,"cm")

  coef.within <- other.variation
  df.between <- n-Kb-1

  if (nrow(X.bet)<=ncol(X.bet)){
    stop("Between estimation impossible (insufficient number of observations)\n")
  }
  
  if(is.null(W)){
    between <- lm(y.bet~X.bet)
  }
  else{
    if (ncol(W)<ncol(X)+1) stop("Insufficient number of instruments\n")
    W.m <- papply(W,mymean,cond)
    W.bet <- attr(W.m,"cm")
    between <- twosls(y.bet,X.bet,W.bet,TRUE)
  }
  rnames <- rownames(X.bet)
  between <- plmformat(between,c("(intercept)",coef.names[cond.variation]),rnames,
                       df.between,"between",pdim,pmodel,indexes,cl)
  between  
}    

plm.pooling <- function(y,X,W,id,time,pvar,pdim,pmodel,indexes,cl,...){
  formula <- pmodel$formula
  effect <- pmodel$effect
  rnames <- rownames(X)

  K <- pdim$K <- ncol(X)
  N <- nrow(X)
  coef.names <- colnames(X)
  
  if (is.null(W)){
    pooling <- lm(y~X)
  }
  else{
    if (ncol(W)<ncol(X)+1) stop("Insufficient number of instruments\n")
    pooling=twosls(y,X,W,TRUE)
  }
  X <- cbind(1,X)
  colnames(X)[1] <- "(intercept)"
  pooling <- plmformat(pooling,c("(intercept)",coef.names),rnames,
                       N-K-1,"pooling",pdim,pmodel,indexes,cl)
  pooling$model[["X"]] <- cbind(1,pooling$model[["X"]])
  colnames(pooling$model[["X"]])[1] <- "(intercept)"
  pooling
}

plm.fd <- function(y,X,W,id,time,pvar,pdim,pmodel,indexes,cl,...){
  formula <- pmodel$formula
  effect <- pmodel$effect
  rnames <- rownames(X)

  K <- pdim$K <- ncol(X)
  N <- nrow(X)
  coef.names <- colnames(X)
  
  if (is.null(W)){
    pooling <- lm(y~X)
  }
  else{
    pooling=twosls(y,X,W,TRUE)
  }
  X <- cbind(1,X)
  colnames(X)[1] <- "(intercept)"
  pooling <- plmformat(pooling,c("(intercept)",coef.names),rnames,
                       N-K-1,"pooling",pdim,pmodel,indexes,cl)
}


plm.random <- function(y,X,W,id,time,pvar,pdim,pmodel,indexes,cl,...){
  random.method <- pmodel$random.method
  formula <- pmodel$formula
  effect <- pmodel$effect
  rnames <- rownames(X)
  inst.method <- pmodel$inst.method
  balanced <- pdim$balanced
  Ti <- pdim$Tint$Ti
  T <- pdim$nT$T
  n <- pdim$nT$n
  N <- pdim$nT$N
  K <- pdim$K <- ncol(X)
  coef.names <- colnames(X)

  cond <- id
  other <- time
  other.variation <- pvar$time.variation
  cond.variation <- pvar$id.variation
  ncond <- n

  if (effect=="time"){
    cond <- time
    other <- id
    cond.variation <- pvar$time.variation
    other.variation <- pvar$id.variation
    ncond <- T
  }
  
  X.m <- papply(X,mymean,cond)
  y.m <- papply(unclass(y),mymean,cond)
  
  if(effect=="twoways"){
    X.time <- papply(X,mymean,time)
    X.mean <- matrix(rep(apply(X,2,mean),N),ncol=K,byrow=T)
    y.time <- papply(unclass(y),mymean,time)
  }

  if (!balanced) random.method <- "swar"
  estec <-switch(random.method,
                 "walhus"=walhus(y,X,W,id,time,pvar,pdim,pmodel,indexes,cl,...),
                 "amemiya"=amemiya(y,X,W,id,time,pvar,pdim,pmodel,indexes,cl,...),
                 "swar"=swar(y,X,W,id,time,pvar,pdim,pmodel,indexes,cl,...),
                 "nerlove"=nerlove(y,X,W,id,time,pvar,pdim,pmodel,indexes,cl,...),               
                 {
                   warning("theta must be one of walhus, amemiya, swar or nerlove")
                   theta="swar"
                   swar(y,X,W,id,time,pvar,pdim,pmodel,indexes,cl,...)
                 }
                 )
  sigma2 <- estec$sigma2
  if (balanced){
    theta <- estec$theta
  }
  else{
    sigma2$one <- (Ti*sigma2$id+sigma2$idios)[cond]
    theta <- 1-sqrt(sigma2$idios/(sigma2$idios+Ti*sigma2$id))
    theta <- theta[cond]
  }
  if(is.null(W)){
    if(effect=="twoways"){
      X.re=cbind(1-theta$id-theta$time+theta$total,X-theta$id*X.m-theta$time*X.time+theta$total*X.mean)
      y.re=y-theta$id*y.m-theta$time*y.time+theta$total*mean(y)
    }  
    else{
      X.re <- cbind(1-theta,X-theta*X.m)
      y.re <- y-theta*y.m
    }
    colnames(X.re)[1] <- "(intercept)"
    random <- lm(y.re~X.re-1)
  }
  else{
    if(effect=="twoways"){
      warning("Instrumental variable random effect estimation not implemented for two-ways panels")
    }
    else{
      if (ncol(W)<ncol(X)+1) stop("Insufficient number of instruments\n")
      X.re <- cbind(1-theta,X-theta*X.m)/sqrt(sigma2$idios)
      y.re <- (y-theta*y.m)/sqrt(sigma2$idios)
      W.m <- papply(W,mymean,cond)
      W.with <- W-W.m
      
      if(inst.method=="baltagi"){
        W.re <- cbind(W.with,W.m)
      }
      if(inst.method=="bvk"){
        W.re <- cbind(W.with/sqrt(sigma2$idios)+W.m/sqrt(sigma2$one))
      }
      random <- twosls(y.re,X.re,W.re)
    }
  }
  random <- plmformat(random,c("(intercept)",coef.names),rnames,
                      N-K-1,"random",pdim,pmodel,indexes,cl)
  random$theta <- theta
  random$sigma2 <- sigma2
  random
}

plm.ht <- function(y,X,W,id,time,pvar,pdim,pmodel,indexes,cl,...){
  within <- plm.within(y,X,NULL,id,time,pvar,pdim,pmodel,indexes,cl,...)
  formula <- pmodel$formula
  effect <- pmodel$effect
  rnames <- rownames(X)
  balanced <- pdim$balanced
  T <- pdim$nT$T
  n <- pdim$nT$n
  N <- pdim$nT$N
  Ti <- pdim$Tint$Ti

  K <- pdim$K <- ncol(X)
  coef.names <- colnames(X)

  cond <- id
  other <- time
  other.variation <- pvar$time.variation
  cond.variation <- pvar$id.variation
  ncond <- n

  if (effect=="time"){
    cond <- time
    other <- id
    cond.variation <- pvar$time.variation
    other.variation <- pvar$id.variation
    ncond <- T
  }

  Kw <- pdim$Kw <- sum(other.variation)
  Kb <- pdim$Kb <- sum(cond.variation)
  Kd <- pdim$Kd <- sum(cond.variation & other.variation)

  inst.names <- colnames(W[,-1,drop=F])
  
  varcond <- function(x,id){
    tapply(x,id,myvar)
  }
  
  inst.time <- apply(apply(W[,-1,drop=F],2,varcond,id),2,sum)!=0

                                        # x are time varying variables
                                        # z are time independent variables
                                        # 1 are exogenous  variables
                                        # 2 are endogenous variables
  x <- coef.names[other.variation]
  z <- coef.names[!other.variation]
  x1 <- inst.names[inst.time]
  z1 <- inst.names[!inst.time]
  z2 <- z[!z%in%z1]
  x2 <- x[!x%in%x1]

  if (length(z2)>length(x1)){
    stop(" The number of endogenous time-invariant variables is greater than the number of exogenous time varying variables\n")
  }
  if (length(x1)!=0) X1 <- X[,x1,drop=F] else X1 <- NULL
  if (length(x2)!=0) X2 <- X[,x2,drop=F] else X2 <- NULL
  if (length(x1)!=0) Z1 <- X[,z1,drop=F] else Z1 <- NULL
  if (length(x1)!=0) Z2 <- X[,z2,drop=F] else Z2 <- NULL


  sigma2 <- list()
  sigma2$one <- 0
  sigma2$idios <- sum(within$res^2)/(N-n)
  if (length(z)!=0){
    zo <- twosls(within$FE,cbind(Z1,Z2),cbind(Z1,X1),TRUE)
  }
  else{
    zo <- lm(within$FE~1)
  }
  ssr <- sum(zo$residuals^2)/N

  if(balanced){
    sigma2$id <- ssr-sigma2$idios/T
    theta <- 1-sqrt(sigma2$idios/(sigma2$idios+T*sigma2$id))
  }
  else{
    sigma2$id <- ssr-sigma2$idios/T
#    sigma2$one <- (Ti*sigma2$id+sigma2$idios)[cond]
    theta <- 1-sqrt(sigma2$idios/(sigma2$idios+Ti*sigma2$id))
    theta <- theta[cond]
  }
  y.ra <- y-theta*papply(unclass(y),mymean,id)
  
  X.bet <- papply(X,mymean,id)
  X.ra <- cbind(1-theta,X-theta*X.bet)

  within.inst <- cbind(X1,X2)
  if (!is.null(within.inst)) within.inst <- within.inst-papply(within.inst,mymean,id)
  
  between.inst <- X1
  if (!is.null(between.inst)) between.inst <- papply(X1,mymean,id)
  
  W <- cbind(within.inst,Z1,between.inst)
  
  y.ra <- as.vector(y.ra)
  rownames(X.ra) <- names(y.ra) <- rownames(W)<- 1:N
  ht <- twosls(y.ra,X.ra,W)
  K <- pdim$K <- ncol(X.ra)
  ht <- plmformat(ht,c("(intercept)",coef.names),rnames,
                  N-K-1,"ht",pdim,pmodel,indexes,cl)
  ht$theta <- theta
  ht$sigma2 <- sigma2
  ht$varlist=list(x1=x1,x2=x2,z1=z1,z2=z2)
  ht
}
