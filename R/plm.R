plm.default <- function(y,X,W,id,time,pvar,pdim,pmodel,...){
  data <- list()
  data$X$X <- X
  data$y$y <- y
  data$id <- id
  data$time <- time

  theta <- pmodel$theta
  model <- pmodel$model
  formula <- pmodel$formula
  effect <- pmodel$effect
  trinst <- pmodel$trinst
  time.variation <- pvar$time.variation
  id.variation <- pvar$id.variation
  balanced <- pdim$balanced

  Ti <- pdim$Tint$Ti
  T <- pdim$nT$T
  n <- pdim$nT$n
  N <- pdim$nT$N

  select.within <- ((is.null(model) || model!="pooling"))
  select.random <- ((is.null(model) || model=="random" ) && trinst!="ht")
  select.ht <- (trinst=="ht")
  select.pooling <- trinst!="ht"
  
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
    pooling <- plmformat(pooling,c("(intercept)",coef.names),formula,N-K-1)
    res$pooling <- pooling
  }
      
  if(select.within){
    data$X$X.m <- X.m <- pmean(X,cond)
    X.with <- (X-X.m)[,other.variation]
    X.bet <- attr(X.m,"cm")[,cond.variation]
    data$y$y.m <- y.m <- pmean(y,cond)
    y.with <- y-y.m
    y.bet <- attr(y.m,"cm")


    
    if(effect=="double"){
      data$X$X.time <- X.time <- pmean(X,time)
      X.mean <- matrix(rep(apply(X,2,mean),N),ncol=K,byrow=T)
      X.with <- (X-X.m-X.time+X.mean)[,cond.variation & other.variation]
      X.bet.time <- attr(X.time,"cm")[,other.variation]
      data$y$y.time <- y.time <- pmean(y,time)
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
    
    if (!is.null(W) & trinst!="ht"){
      data$W$W.m <- W.m <- pmean(W,cond)
      W.bet <- attr(W.m,"cm")
      if(effect=="double"){
        data$W$W.time <- W.time <- pmean(W,time)
        W.mean <- matrix(rep(apply(W,2,mean),N),ncol=K,byrow=T)
        W.with <- (W-W.m-W.time+W.mean)[,cond.variation & other.variation]
        W.bet.time <- attr(W.time,"cm")
      }
      else{
        W.with <- W-W.m
      }
    }
    
    attr(res,"data") <- data
    
    if(is.null(W) || trinst=="ht"){

      between <- lm(y.bet~X.bet)
      within <- lm(y.with~X.with-1)
      if (effect=="double"){
        between.time <- lm(y.bet.time~X.bet.time)
        FE.time <- y.time-as.vector(X.time[,time.variation & id.variation,drop=F]%*%within$coef)
      }
    }
    else{
      between <- twosls(y.bet,X.bet,W.bet,TRUE)

      between <- plmformat(between,c("(intercept)",coef.names[cond.variation]),formula,df.between)
      within <- twosls(y.with,X.with,W.with)

      if (effect=="double"){
        between.time <- twosls(y.bet.time,X.bet.time,W.bet.time,TRUE)
      }
    }
    between <- plmformat(between,c("(intercept)",coef.names[cond.variation]),formula,df.between)
    within <- plmformat(within,coef.names[coef.within],formula,df.within)
    within$FE <- y.m-as.vector(X.m[,coef.within,drop=F]%*%within$coef)
    attr(within$FE,"cm") <- tapply(within$FE,cond,mean)
    within$alpha <- mean(y)-apply(X[,coef.within,drop=F],2,mean)%*%within$coef
    if(effect=="double"){
      within$FE.time <- FE.time
      between.time <- plmformat(between.time,c("(intercept)",coef.names[other.variation]),formula,df.between.time)
      res$between.time <- between.time
    }
    res$within <- within
    res$between <- between
  }    
    
  if(select.random){

    estec <-switch(theta,
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

    if (balanced){
      thetahat <- estec$theta
    }
    else{
      thetai <- 1-sqrt(estec$s2eta/(estec$s2eta+Ti*estec$s2mu))
      thetahat <- thetai[cond]
    }
    
    if(is.null(W)){
      if(effect=="double"){
        X.re=cbind(1-estec$theta1-estec$theta2+estec$theta3,X-estec$theta1*X.m-estec$theta2*X.time+estec$theta3*X.mean)
        y.re=y-estec$theta1*y.m-estec$theta2*y.time+estec$theta3*mean(y)
      }  
      else{
        X.re <- cbind(1-thetahat,X-thetahat*X.m)
        y.re <- y-thetahat*y.m
      }
      random <- lm(y.re~X.re-1)
    }
    else{
      if(effect=="double"){
        warning("Instrumental variable random effect estimation not implemented for two.ways panels")
      }
      else{
        estec <- swar(res)
        thetahat <- estec$theta
        s21 <- estec$s21
        s2eta <- estec$s2eta
        X.re <- cbind(1-thetahat,X-thetahat*X.m)/sqrt(s2eta)
        y.re <- (y-thetahat*y.m)/sqrt(s2eta)
        if(trinst=="baltagi"){
          W.re <- cbind(W.with,W.m)
        }
        if(trinst=="bvk"){
          W.re <- cbind(W.with/sqrt(s2eta)+W.m/sqrt(s21))
        }
        random <- twosls(y.re,X.re,W.re)
      }
    }
    random <- plmformat(random,c("(intercept)",coef.names),formula,N-K-1)
    random$thetahat <- thetahat
    random$theta <- theta
  }

  if(select.ht){
    inst.names <- colnames(W[,-1])
    inst.time <- apply(apply(W[,-1,drop=F],2,varcond,id),2,sum)!=0
    inst.id <- apply(apply(W[,-1,drop=F],2,varcond,time),2,sum)!=0
    x1 <- inst.names[inst.time]
    z1 <- inst.names[!inst.time]
    x <- coef.names[time.variation]
    z <- coef.names[!time.variation]
    x2 <- x[!x%in%x1]
    z2 <- z[!z%in%z1]
    X1 <- X[,x1]
    X2 <- X[,x2]
    Z1 <- X[,z1]
    Z2 <- X[,z2]
    s2eta <- sum(within$res^2)/(N-n)
    V <- cbind(Z1,Z2)
    Vi <- lm(V~Z1+X1)$fit
    zo <- twosls(within$FE,V,cbind(Z1,X1),TRUE)
    ssr <- sum(zo$residuals^2)/N
    sigmau <- ssr-s2eta/T
    tetha <- 1-sqrt(s2eta/(s2eta+T*sigmau))
    y.ra <- y-tetha*pmean(y,id)
    X.bet <- pmean(X,id)
    X.ra <- cbind(1-tetha,X-tetha*X.bet)
    Z <- cbind(X1,X2)
    W <- cbind(Z-pmean(Z,id),Z1,pmean(X1,id))
    y.ra <- as.vector(y.ra)
    rownames(X.ra) <- names(y.ra) <- rownames(W)<- 1:N
    z <- twosls(y.ra,X.ra,W)
    K <- ncol(X.ra)
    random <- plmformat(z,c("(intercept)",coef.names),formula,N-K-1)
    random
  }
  
  if (is.null(model) & trinst!="ht"){
    res <- list(pooling=pooling,between=between,within=within,random=random)
    if(effect=="double") res$between.time <- between.time
    res <- structure(res,data=data,pdim=pdim,pmodel=pmodel,class="plms")
  }
  if(!is.null(model)){
    res <- switch(model,
                  "pooling"=pooling,
                  "within"=within,
                  "random"=random,
                  "between"=between
                  )
    res <- structure(res,data=data,pdim=pdim,pmodel=pmodel,class="plm")
  }
  if(trinst=="ht"){
    res <- structure(random,data=data,pdim=pdim,pmodel=pmodel,class="plm")
  }
  res
}
