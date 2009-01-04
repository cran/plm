plm.within <- function(y, X, W, id, time, pvar, pdim, pmodel, indexes, cl, ...){
  effect <- pmodel$effect
  rnames <- rownames(X)
  T <- pdim$nT$T ; n <- pdim$nT$n ; N <- pdim$nT$N ; K <- pdim$K <- ncol(X)
  coef.names <- colnames(X)
  if(effect == "time"){
    cond <- time ; other <- id ; cond.variation <- pvar$time.variation
    other.variation <- pvar$id.variation ; ncond <- T
  }
  else{
    cond <- id ; other <- time ; cond.variation <- pvar$id.variation
    other.variation <- pvar$time.variation ; ncond <- n

  }
  Kw <- pdim$Kw <- sum(other.variation)
  Kd <- pdim$Kd <- sum(cond.variation & other.variation)
  if (Kw == 0) stop("the within model can't be computed because all the variables are individual specific")
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
    df.within <- N-ncond-Kw
  }
  if (nrow(X.with) <= ncol(X.with)) stop("Within estimation impossible (insufficient number of observations)\n")
  
  if(is.null(W)) within <- lm(y.with~X.with-1)
  else{
    if (ncol(W) < ncol(X+1)) stop("Insufficient number of instruments\n")
    W.m <- W.m <- papply(W,mymean,cond)
    if(effect=="twoways"){
      W.time <- papply(W,mymean,time)
                                        #K+1 et non K W contient une constante
      W.mean <- matrix(rep(apply(W,2,mean),N),ncol=K+1,byrow=T)
      W.with <- (W-W.m-W.time+W.mean)[,cond.variation & other.variation]
    }
    else W.with <- W-W.m
    within <- twosls(y.with,X.with,W.with)
  }
  within <- plmformat(within,coef.names[coef.within],rnames,
                      df.within,"within",pdim,pmodel,indexes,cl)
  fixef.id <- y.m-as.vector(X.m[,coef.within,drop=F]%*%within$coef)
  attr(fixef.id,"cm") <- tapply(fixef.id,cond,mean)
  within$alpha <- as.vector(mean(y)-apply(X[,coef.within,drop=F],2,mean)%*%within$coef)
  if(effect=="twoways"){
    fixef.time <- y.time-as.vector(X.time[,cond.variation & other.variation,drop=F]%*%within$coef)
    within$fixef <- list(id=fixef.id,time=fixef.time)
  }
  else within$fixef <- fixef.id
  within
}    

plm.between <- function(y, X, W, id, time, pvar, pdim, pmodel, indexes, cl, ...){
  response.name <- deparse(pmodel$formula[[2]])
  effect <- pmodel$effect
  interc <- pmodel$has.intercept
  balanced <- pdim$balanced
  T <- pdim$nT$T ; n <- pdim$nT$n ; N <- pdim$nT$N ; K <- pdim$K <- ncol(X)
  coef.names <- colnames(X)
  if(effect == "time"){
    cond <- time ; other <- id ; cond.variation <- pvar$time.variation
    other.variation <- pvar$id.variation ; ncond <- T
  }
  else{
    cond <- id ; other <- time ; cond.variation <- pvar$id.variation
    other.variation <- pvar$time.variation ; ncond <- n

  }
  Kb <- pdim$Kb <- sum(cond.variation)
  X.m <- papply(X,mymean,cond)
  X.bet <- attr(X.m,"cm")[,cond.variation,drop=F]
  if (ncol(X.bet) == 0) X.bet <- NULL
  if (is.null(X.bet)) stop("the between model is empty")
  y.m <- papply(unclass(y),mymean,cond)
  y.bet <- attr(y.m,"cm")
  coef.within <- other.variation
  df.between <- n-Kb-1
  if (nrow(X.bet)<=ncol(X.bet)){
    stop("Between estimation impossible (insufficient number of observations)\n")
  }
  if(is.null(W)){
    if (interc) between <- lm(y.bet~X.bet) else between <- lm(y.bet~X.bet-1)
  }
  else{
    if (ncol(W)<ncol(X)+1) stop("Insufficient number of instruments\n")
    W.m <- papply(W,mymean,cond)
    W.bet <- attr(W.m,"cm")
    between <- twosls(y.bet,X.bet,W.bet,interc)
  }
  rnames <- rownames(X.bet)
  between <- plmformat(between,coef.names[cond.variation],rnames,
                       df.between,"between",pdim,pmodel,indexes,cl)
  between  
}    

plm.pooling <- function(y, X, W, id, time, pvar, pdim, pmodel, indexes, cl, ...){
  response.name <- deparse(pmodel$formula[[2]])
  effect <- pmodel$effect
  interc <- pmodel$has.intercept
  rnames <- rownames(X)
  K <- pdim$K <- ncol(X)
  N <- nrow(X)
  coef.names <- colnames(X)
  if (is.null(W)){
    if (interc) pooling <- lm(y~X) else pooling <- lm(y~X-1)
  }
  else{
    if (ncol(W)<ncol(X)+1) stop("Insufficient number of instruments\n")
    pooling <- twosls(y,X,W,interc)
  }
  pooling <- plmformat(pooling, coef.names, rnames,N-K-interc,"pooling",
                       pdim, pmodel, indexes, cl)
  pooling
}

plm.fd <- function(y, X, W, id, time, pvar, pdim, pmodel, indexes, cl, ...){
  interc <- pmodel$has.intercept
  response.name <- deparse(pmodel$formula[[2]])
  effect <- pmodel$effect
  rnames <- rownames(X)
  T <- pdim$nT$T ; n <- pdim$nT$n ; N <- pdim$nT$N ; K <- pdim$K <- ncol(X)
  coef.names <- colnames(X)
  Ti <- pdim$Tint$Ti-1
  nt <- pdim$Tint$nt[-1]
  N <- pdim$nT$N <- N-n
  time.names <- pdim$panel.names$time.names[-1]
  pdim$nT <- list(n=n,T=T,N=N)
  pdim$Tint <- list(Ti=Ti,nt=nt)
  pdim$panel.names$time.names <- time.names[-1]
  K <- pdim$K <- ncol(X)
  N <- nrow(X)
  coef.names <- colnames(X)
  X <- rbind(NA,X[2:N,,drop=FALSE]-X[1:(N-1),,drop=FALSE])
  y <- c(NA,y[2:N]-y[1:(N-1)])
  did <- c(1,as.numeric(id[2:N])-as.numeric(id[1:(N-1)]))
  X <- X[did==0,]
  y <- y[did==0]
  rnames <- rnames[did==0]
  
  if (is.null(W)) if (interc) fd <- lm(y~X) else fd <- lm(y~X-1)
  else{
    W <- rbind(NA,W[2:N,]-W[1:(N-1),])
    W <- W[did==0,]
    fd <- twosls(y,X,W,interc)
  }
  fd <- plmformat(fd,coef.names,rnames,
                  N-K-n,"pooling",pdim,pmodel,indexes,cl)
  fd
}

plm.random <- function(y, X, W, id, time, pvar, pdim, pmodel, indexes, cl, ...){
  response.name <- deparse(pmodel$formula[[2]])
  random.method <- pmodel$random.method
  interc <- pmodel$has.intercept
  effect <- pmodel$effect
  rnames <- rownames(X)
  inst.method <- pmodel$inst.method
  balanced <- pdim$balanced
  Ti <- pdim$Tint$Ti ;  nt <- pdim$Tint$nt ; T <- pdim$nT$T ; n <- pdim$nT$n ;
  N <- pdim$nT$N ; K <- pdim$K <- ncol(X)
  coef.names <- colnames(X)
  if (effect=="time"){
    cond <- time ; other <- id ; cond.variation <- pvar$time.variation ;
    other.variation <- pvar$id.variation ; ncond <- T ; ns <- nt
  }
  else{
    cond <- id ; other <- time ; cond.variation <- pvar$id.variation ;
    other.variation <- pvar$time.variation ; ncond <- n ; ns <- Ti
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
                 "nerlove"=nerlove(y,X,W,id,time,pvar,pdim,pmodel,indexes,cl,...)
                 )

  sigma2 <- estec$sigma2
  if (balanced) theta <- estec$theta
  else{
    sigma2$one <- (ns*sigma2$id+sigma2$idios)[cond]
    theta <- 1-sqrt(sigma2$idios/(sigma2$idios+ns*sigma2$id))
    theta <- theta[cond]
  }
  if(is.null(W)){
    if(effect=="twoways"){
      X.re <- cbind(1-theta$id-theta$time+theta$total,X-theta$id*X.m-theta$time*X.time+theta$total*X.mean)
      y.re <- y-theta$id*y.m-theta$time*y.time+theta$total*mean(y)
    }  
    else{
      if (interc){
        X.re <- cbind(1-theta,X-theta*X.m)
        colnames(X.re)[1] <- "(intercept)"
      }
      else X.re <- X-theta*X.m
      rownames(X.re) <- rownames(X)
      y.re <- y-theta*y.m
    }
    random <- lm(y.re~X.re-1)
  }
  else{
    if(effect=="twoways"){
      warning("Instrumental variable random effect estimation not implemented for two-ways panels")
    }
    else{
      if (ncol(W) < ncol(X)+1) stop("Insufficient number of instruments\n")
      if (interc) X.re <- cbind(1-theta,X-theta*X.m)/sqrt(sigma2$idios)
      else X.re <- (X-theta*X.m)/sqrt(sigma2$idios)
      rownames(X.re) <- rownames(X)
      y.re <- (y-theta*y.m)/sqrt(sigma2$idios)
      W.m <- papply(W,mymean,cond)
      W.with <- W-W.m
      if(inst.method=="baltagi") W.re <- cbind(W.with,W.m)
      if(inst.method=="bvk") W.re <- cbind(W.with/sqrt(sigma2$idios)+W.m/sqrt(sigma2$one))
      random <- twosls(y.re,X.re,W.re,FALSE)
    }
  }
  random <- plmformat(random,coef.names,rnames,
                      N-K-1,"random",pdim,pmodel,indexes,cl)
  random$theta <- theta
  random$sigma2 <- sigma2
  random
}

plm.ht <- function(y, X, W, id, time, pvar, pdim, pmodel, indexes, cl, ...){
  var.effects <- attr(X,"var.effects")
  within <- plm.within(y,X,NULL,id,time,pvar,pdim,pmodel,indexes,cl,...)
  effect <- pmodel$effect
  rnames <- rownames(X)
  interc <- pmodel$has.intercept
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
    zo <- twosls(within$fixef,cbind(Z1,Z2),cbind(Z1,X1),TRUE)
  }
  else{
    zo <- lm(within$fixef~1)
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
  ht <- plmformat(ht,coef.names,rnames,
                  N-K-1,"ht",pdim,pmodel,indexes,cl)
  ht$theta <- theta
  ht$sigma2 <- sigma2
  
  transl <- function(x,y){
    co <- match(x,names(y))
    x[!is.na(co)] <- y[na.omit(co)]
    x
  }

  ht$varlist <- list(x1=transl(x1,var.effects),
                     x2=transl(x2,var.effects),
                     z1=transl(z1,var.effects),
                     z2=transl(z2,var.effects))
  
  ht
}
