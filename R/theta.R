swar <- function(y,X,W,id,time,pvar,pdim,pmodel,indexes,cl,...){
  effect <- pmodel$effect
  twoways <- pmodel$effect=="twoways"
  balanced <- pdim$balanced
  T <- pdim$nT$T
  n <- pdim$nT$n
  Kb <- pdim$Kb
  K <- pdim$K
  N <- pdim$nT$N
  Ti <- pdim$Tint$Ti
  within <- plm.within(y,X,W,id,time,pvar,pdim,pmodel,indexes,cl,...)
  Kw <- attr(within,"pdim")$Kw
  if (twoways){
    pmodel.time <- pmodel.id <- pmodel
    pmodel.time$effect <- "time"
    pmodel.id$effect <- "id"
    between.id <- plm.between(y,X,W,id,time,pvar,pdim,pmodel.id,indexes,cl,...)
    between.time <- plm.between(y,X,W,id,time,pvar,pdim,pmodel.time,indexes,cl,...)
    Kb <- attr(between.id,"pdim")$Kb
  }
  else{
    between <- plm.between(y,X,W,id,time,pvar,pdim,pmodel,indexes,cl,...)
    Kb <- attr(between,"pdim")$Kb
  }

  sigma2 <- list()
  if(!twoways){
    if(balanced){
      ssrbet <- T*sum(between$residuals^2)
      sigma2$one <- ssrbet/(n-Kb-1)
      sigma2$idios <- sum(within$residuals^2)/(N-n-Kw)
      sigma2$id <- (sigma2$one-sigma2$idios)/T
      theta <- 1-sqrt(sigma2$idios/sigma2$one)
      z <- list(sigma2=sigma2,theta=theta)
    }
    else{
      data.name <- within$call$data
      if (effect=="individual"){
        condvar <- eval(data.name)[[indexes$id]]
      }
      else{
        condvar <- eval(data.name)[[indexes$time]]
      }
      X.m <- papply(cbind(1,X),mymean,condvar)
      X.sum <- attr(papply(cbind(1,X),mysum,condvar),"cm")
#      tr <- sum(diag(solve(crossprod(cbind(1,X.m)))%*%crossprod(X.sum))
      tr <- sum(diag(solve(crossprod(X.m))%*%crossprod(X.sum)))
      sigma2$idios <- sum(within$residuals^2)/(N-n-Kw)
      ssrbet <- sum(between$residuals^2*Ti)
      sigma2$id <- (ssrbet-(n-Kb-1)*sigma2$idios)/(N-tr)
      z <- list(sigma2=sigma2)
    }
  }
  else{
    if(balanced){
      theta=list()
      sigma2$idios <- sum(within$residuals^2)/((n-1)*(T-1)-Kw)
      lambda2 <- T*sum(between.id$residuals^2)/(n-Kb-1)
      lambda3 <- n*sum(between.time$residuals^2)/(T-Kb-1)
      lambda4 <- lambda2+lambda3-sigma2$idios
      sigma2$id <- (lambda2-sigma2$idios)/T
      sigma2$time <- (lambda3-sigma2$idios)/n
      theta$id <- 1-sqrt(sigma2$idios/lambda2)
      theta$time <- 1-sqrt(sigma2$idios/lambda3)
      theta$total <- theta$id+theta$time+sqrt(sigma2$idios/lambda4)-1
      if (sigma2$time<0) theta$time <- theta$total <- sigma2$time <- 0
      z <- list(theta=theta,sigma2=sigma2)
    }
    else{
      stop("twoway random effect model not implemented for unbalanced panels")
    }
  }
  z
}

walhus <- function(y,X,W,id,time,pvar,pdim,pmodel,indexes,cl,...){
  twoways <- pmodel$effect=="twoways"
  effect <- pmodel$effect
  balanced <- pdim$balanced
  pooling <- plm.pooling(y,X,W,id,time,pvar,pdim,pmodel,indexes,cl,...)
  T <- pdim$nT$T
  n <- pdim$nT$n
  sigma2=list()
  if(!balanced){stop("walhus not implemented for unbalanced panels\n")}
  if(!twoways){
    if (effect=="individual") condvar <- id else condvar <- time
    resi <- pooling$residuals
    sigma2$one <- T*sum(tapply(resi,condvar,mean)^2)/n
    sigma2$idios <- sum((resi-tapply(resi,condvar,mean)[as.character(condvar)])^2)/(n*(T-1))
    sigma2$id <- (sigma2$one-sigma2$idios)/T
    theta <- 1-sqrt(sigma2$idios/sigma2$one)
    z <- list(theta=theta,sigma2=sigma2)
  }
  else{
    theta=list()
    sigma2$idios <- sum((pooling$residuals-tapply(pooling$residuals,id,mean)[as.character(id)]-
                         tapply(pooling$residuals,time,mean)[as.character(time)])^2)/((n-1)*(T-1))
    lambda2 <- sum(tapply(pooling$residuals,id,mean)^2)*T/(n-1)
    lambda3 <- sum(tapply(pooling$residuals,time,mean)^2)*n/(T-1)
    lambda4 <- lambda2+lambda3-sigma2$idios
    sigma2$id <- (lambda2-sigma2$idios)/T
    sigma2$time <- (lambda3-sigma2$idios)/n
    theta$id <- 1-sqrt(sigma2$idios/lambda2)
    theta$time <- 1-sqrt(sigma2$idios/lambda3)
    theta$total <- theta$id+theta$time+sqrt(sigma2$idios/lambda4)-1
    if (sigma2$time<0) theta$time <- theta$total <- sigma2$time <- 0    
    z <- list(theta=theta,sigma2=sigma2)
  }
  z
}

amemiya <- function(y,X,W,id,time,pvar,pdim,pmodel,indexes,cl,...){
  twoways <- pmodel$effect=="twoways"
  effect <- pmodel$effect
  balanced <- pdim$balanced
  within <- plm.within(y,X,W,id,time,pvar,pdim,pmodel,indexes,cl,...)
  T <- pdim$nT$T
  n <- pdim$nT$n
  Kw <- sum(pvar$time.variation)
  Kb <- sum(pvar$id.variation)
  sigma2 <- list()
  if(!balanced)stop("amemiya variance decomposition not implemented for unbalanced panels")
  if(!twoways){
    if (effect=="individual") condvar <- id else condvar <- time
    uest <- within$res+within$FE-within$alpha
    sigma2$one <- T/n*sum(tapply(uest,condvar,mean)^2)
    sigma2$idios <- sum(within$residuals^2)/(n*(T-1)-Kw)
    sigma2$id <- max((sigma2$one-sigma2$idios)/T,0)
    theta <- max(1-sqrt(sigma2$idios/sigma2$one),0)
    z <- list(theta=theta,sigma2=sigma2)
  }
  else{
    theta=list()
    uest <- as.vector(y-within$alpha-X%*%within$coef)
    sigma2$idios <- sum((uest-tapply(uest,id,mean)[as.character(id)]-
                         tapply(uest,time,mean)[as.character(time)])^2)/((n-1)*(T-1))
    lambda2 <- sum(tapply(uest,id,mean)^2)*T/(n-1)
    lambda3 <- sum(tapply(uest,time,mean)^2)*n/(T-1)
    lambda4 <- lambda2+lambda3-sigma2$idios
    sigma2$id <- (lambda2-sigma2$idios)/T
    sigma2$time <- (lambda3-sigma2$idios)/n
    theta$id <- 1-sqrt(sigma2$idios/lambda2)
    theta$time <- 1-sqrt(sigma2$idios/lambda3)
    theta$total <- theta$id+theta$time+sqrt(sigma2$idios/lambda4)-1
    if (sigma2$time<0) theta$time <- theta$total <- sigma2$time <- 0
    z <- list(theta=theta,sigma2=sigma2)
  }
  z
}

nerlove <- function(y,X,W,id,time,pvar,pdim,pmodel,indexes,cl,...){
  within <- plm.within(y,X,W,id,time,pvar,pdim,pmodel,indexes,cl,...)
  twoways <- pmodel$effect=="twoways"
  effect <- pmodel$effect
  balanced <- pdim$balanced
  n <- pdim$nT$n
  N <- pdim$nT$N
  Kb <- pdim$Kb
  Kw <- pdim$Kw
  sigma2=list()
  if(!twoways & balanced){
    sigma2$idios <- sum(within$residuals^2)/N
    sigma2$id <- sum((within$FE-mean(within$FE))^2)/(n-1)
    sigma2$one <- T*sigma2$id+sigma2$idios
    theta <- 1-sqrt(sigma2$idios/sigma2$one)
    z <- list(theta=theta,sigma2=sigma2)
  }
  else stop("nerlove variance decomposition only implemented for balanced oneway panels")
  z
}

