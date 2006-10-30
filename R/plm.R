plm.fit <- function(y,X,W,id,time,pvar,pdim,pmodel,indexes,cl,...){
  data <- list()
  data$X$X <- X
  data$y$y <- y
  data$id <- id
  data$time <- time
  random.method <- pmodel$random.method
  model.name <- pmodel$model.name
  formula <- pmodel$formula
  effect <- pmodel$effect
  inst.method <- pmodel$inst.method
  time.variation <- pvar$time.variation
  id.variation <- pvar$id.variation
  balanced <- pdim$balanced

  Ti <- pdim$Tint$Ti
  T <- pdim$nT$T
  n <- pdim$nT$n
  N <- pdim$nT$N

  if (is.null(model.name)) model.name <- "all"
  select.within <- (model.name!="pooling")
  select.random <- (model.name=="random" || model.name=="all")
  select.ht <- (model.name=="ht")
  select.pooling <- (model.name!="ht")

  K <- pdim$K <- ncol(X)
  coef.names <- colnames(X)

  cond <- id
  other <- time
  other.variation <- time.variation
  cond.variation <- id.variation
  ncond <- n

  if (effect=="time"){
    cond <- time
    other <- id
    cond.variation <- time.variation
    other.variation <- id.variation
    ncond <- T
  }

  Kw <- pdim$Kw <- sum(other.variation)
  Kb <- pdim$Kb <- sum(cond.variation)
  Kd <- pdim$Kd <- sum(cond.variation & other.variation)

  res <- list()

  attr(res,"pdim") <- pdim
  attr(res,"pmodel") <- pmodel
  attr(res,"pvar") <- pvar

  if(select.pooling){
    if (is.null(W)) pooling <- lm(y~X) else pooling=twosls(y,X,W,TRUE)
    pooling <- plmformat(pooling,c("(intercept)",coef.names),formula,
                         N-K-1,"pooling",tss(y),pdim,pmodel,indexes,cl)
    res$pooling <- pooling
  }

  if(select.within){
    data$X$X.m <- X.m <- papply(X,mymean,cond)
    X.with <- (X-X.m)[,other.variation,drop=F]
    X.bet <- attr(X.m,"cm")[,cond.variation,drop=F]
    data$y$y.m <- y.m <- papply(unclass(y),mymean,cond)
    y.with <- y-y.m
    y.bet <- attr(y.m,"cm")
    
    if(effect=="twoways"){
      data$X$X.time <- X.time <- papply(X,mymean,time)
      X.mean <- matrix(rep(apply(X,2,mean),N),ncol=K,byrow=T)
      X.with <- (X-X.m-X.time+X.mean)[,cond.variation & other.variation]
      X.bet.time <- attr(X.time,"cm")[,other.variation]
      data$y$y.time <- y.time <- papply(unclass(y),mymean,time)
      y.with <- y-y.m-y.time+mean(y)
      y.bet.time <- attr(y.time,"cm")
      df.within <- N-T-n-Kd+1
      df.between <- n-Kb-1
      df.between.time <- T-Kw-1
      coef.within <- other.variation & cond.variation
    }
    else{
      coef.within <- other.variation
      df.within <- N-n-Kw
      df.between <- n-Kb-1
    }
    
    if (!is.null(W) & model.name!="ht"){
      data$W$W.m <- W.m <- papply(W,mymean,cond)
      W.bet <- attr(W.m,"cm")
      if(effect=="twoways"){
        data$W$W.time <- W.time <- papply(W,mymean,time)
        W.mean <- matrix(rep(apply(W,2,mean),N),ncol=K,byrow=T)
        W.with <- (W-W.m-W.time+W.mean)[,cond.variation & other.variation]
        W.bet.time <- attr(W.time,"cm")
      }
      else{
        W.with <- W-W.m
      }
    }
    
    attr(res,"data") <- data
    
    if(is.null(W) || model.name=="ht"){
      if (nrow(X.bet)<=ncol(X.bet)){
        stop("Between estimation impossible (insufficient number of observations)\n")
      }
      else{
        between <- lm(y.bet~X.bet)
      }
      if (nrow(X.with)<=ncol(X.with)){
        stop("Within estimation impossible (insufficient number of observations)\n")
      }
      else{
        within <- lm(y.with~X.with-1)
      }
      if (effect=="twoways"){
        between.time <- lm(y.bet.time~X.bet.time)
        FE.time <- y.time-as.vector(X.time[,time.variation & id.variation,drop=F]%*%within$coef)
      }
    }
    else{
      if (nrow(X.bet)<=ncol(X.bet)){
        stop("Between estimation impossible (insufficient number of observations)\n")
      }
      else{
        between <- twosls(y.bet,X.bet,W.bet,TRUE)
      }
      if (nrow(X.with)<=ncol(X.with)){
        stop("Within estimation impossible (insufficient number of observations)\n")
      }
      else{
        within <- twosls(y.with,X.with,W.with)
      }
      if (effect=="twoways"){
        between.time <- twosls(y.bet.time,X.bet.time,W.bet.time,TRUE)
      }
    }
    between <- plmformat(between,c("(intercept)",coef.names[cond.variation]),formula,
                         df.between,"between",tss(y.bet),pdim,pmodel,indexes,cl)
    within <- plmformat(within,coef.names[coef.within],formula,
                        df.within,"within",tss(y.with),pdim,pmodel,indexes,cl)
    within$FE <- y.m-as.vector(X.m[,coef.within,drop=F]%*%within$coef)
    attr(within$FE,"cm") <- tapply(within$FE,cond,mean)
    within$alpha <- as.vector(mean(y)-apply(X[,coef.within,drop=F],2,mean)%*%within$coef)
    if(effect=="twoways"){
      within$FE.time <- FE.time
      between.time <- plmformat(between.time,c("(intercept)",coef.names[other.variation]),formula,
                                df.between.time,"between.time",tss(y.bet.time),pdim,pmodel,indexes,cl)
      res$between.time <- between.time
    }
    res$within <- within
    res$between <- between

  }    

  if(select.random){
    estec <-switch(random.method,
                   "walhus"=walhus(res),
                   "amemiya"=amemiya(res),
                   "swar"=swar(res),
                   "nerlove"=nerlove(res),
                   {
                     warning("theta must be one of walhus, amemiya, swar or nerlove")
                     theta="walhus"
                     walhus(res)
                   }
                   )
    sigma2 <- estec$sigma2

    if (balanced){
      theta <- estec$theta
    }
    else{
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
      random <- lm(y.re~X.re-1)
    }
    else{
      if(effect=="twoways"){
        warning("Instrumental variable random effect estimation not implemented for two.ways panels")
      }
      else{
        estec <- swar(res)
        theta <- estec$theta
        sigma2 <- estec$sigma2
        X.re <- cbind(1-theta,X-theta*X.m)/sqrt(sigma2$idios)
        y.re <- (y-theta*y.m)/sqrt(sigma2$idios)
        if(inst.method=="baltagi"){
          W.re <- cbind(W.with,W.m)
        }
        if(inst.method=="bvk"){
          W.re <- cbind(W.with/sqrt(sigma2$idios)+W.m/sqrt(sigma2$one))
        }
        random <- twosls(y.re,X.re,W.re)
      }
    }
    random <- plmformat(random,c("(intercept)",coef.names),formula,
                        N-K-1,"random",tss(y.re),pdim,pmodel,indexes,cl)
    random$theta <- theta
    random$sigma2 <- sigma2
  }

  if(select.ht){
    inst.names <- colnames(W[,-1,drop=F])

    varcond <- function(x,id){
      tapply(x,id,myvar)
    }

    inst.time <- apply(apply(W[,-1,drop=F],2,varcond,id),2,sum)!=0

    # x are time varying variables
    # z are time independent variables
    # 1 are exogenous  variables
    # 2 are endogenous variables
    x <- coef.names[time.variation]
    z <- coef.names[!time.variation]
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
      ssr <- sum(zo$residuals^2)/N
      sigma2$id <- ssr-sigma2$idios/T
      theta <- 1-sqrt(sigma2$idios/(sigma2$idios+T*sigma2$id))
    }
    else{
      sigma2 <- swar(res)$sigma2
      theta <- swar(res)$theta
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
    K <- ncol(X.ra)
    ht <- plmformat(ht,c("(intercept)",coef.names),formula,
                    N-K-1,"ht",tss(y.ra),pdim,pmodel,indexes,cl)
    ht$theta <- theta
    ht$sigma2 <- sigma2
    ht$varlist=list(x1=x1,x2=x2,z1=z1,z2=z2)
    ht
  }
  
  if (model.name=="all"){
    res <- list(pooling=pooling,between=between,within=within,random=random)
    if(effect=="twoways") res$between.time <- between.time
    res <- structure(res,data=data,pdim=pdim,pmodel=pmodel,class="plms")
  }
  else{
    res <- switch(model.name,
                  "pooling"=pooling,
                  "within"=within,
                  "random"=random,
                  "between"=between,
                  "ht"=ht
                  )
    res <- structure(res,data=data,pdim=pdim,pmodel=pmodel,class=c("plm","lm"))
  }
  res
}

