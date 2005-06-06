swar <-    function(res){
  data <- attr(res,"data")
  pdim <- attr(res,"pdim")
  pmodel <- attr(res,"pmodel")
  double <- pmodel$effect=="double"
  balanced <- pdim$balanced
  effect <- pmodel$effect
  within <- res$within
  T <- pdim$nT$T
  n <- pdim$nT$n
  Kb <- pdim$Kb
  Kw <- pdim$Kw
  N <- pdim$nT$N
  Ti <- pdim$Tint$Ti
  id <- data$id
  time <- data$time
  if(!double){
    between <- res$between
    if(balanced){
      ssrbet <- T*sum(between$residuals^2)
      s21 <- ssrbet/(n-Kb-1)
      s2eta <- sum(within$residuals^2)/(N-n-Kw)
      s2mu <- (s21-s2eta)/T
      theta <- 1-sqrt(s2eta/s21)
      z <- list(s21=s21,s2eta=s2eta,s2mu=s2mu,theta=theta)
    }
    else{
      if (effect=="individual") condvar <- id else condvar <- time
      X <- data$X$X
      X.m <- data$X$X.m
      X.sum <- attr(psum(cbind(1,X),condvar),"cm")
      tr <- sum(diag(solve(crossprod(cbind(1,X.m)))%*%crossprod(X.sum)))
      s2eta <- sum(within$residuals^2)/(N-n-Kw)
      ssrbet <- sum(between$residuals^2*Ti)
      s2mu <- (ssrbet-(n-Kb-1)*s2eta)/(N-tr)
      z <- list(s2eta=s2eta,s2mu=s2mu)
    }
  }
  else{
    if(balanced){
      between <- res$between
      between.time <- res$between.time
      s2eta <- sum(within$residuals^2)/((n-1)*(T-1)-Kw)
      lambda2 <- T*sum(between$residuals^2)/(n-Kb-1)
      lambda3 <- n*sum(between.time$residuals^2)/(T-Kb-1)
      lambda4 <- lambda2+lambda3-s2eta
      theta1 <- 1-sqrt(s2eta/lambda2)
      theta2 <- 1-sqrt(s2eta/lambda3)
      theta3 <- theta1+theta2+sqrt(s2eta/lambda4)-1
      if (theta2<0) theta2 <- theta3 <- 0
      z <- list(theta1=theta1,theta2=theta2,theta3=theta3)
    }
    else{
      stop("twoway random effect model not implemented for unbalanced panels")
    }
  }
  z
}



walhus <- function(res){
  data <- attr(res,"data")
  pdim <- attr(res,"pdim")
  pmodel <- attr(res,"pmodel")
  double <- pmodel$effect=="double"
  effect <- pmodel$effect
  balanced <- pdim$balanced
  pooling <- res$pooling
  id <- data$id
  time <- data$time
  T <- pdim$nT$T
  n <- pdim$nT$n
  if(!balanced){stop("walhus not implemented for unbalanced panels\n")}
  if(!double){
    if (effect=="individual") condvar <- id else condvar <- time
    s21 <- T/n*sum(tapply(pooling$residuals,condvar,mean)^2)
    resi <- pooling$residuals
    s2eta <- sum((resi-tapply(resi,id,mean)[id])^2)/(n*T+3)
    s21 <- T*sum(tapply(resi,id,mean)^2)/(n-3)
    s2mu <- (s21-s2eta)/T
    theta <- 1-sqrt(s2eta/s21)
    z <- list(theta=theta,s2eta=s2eta,s21=s21,s2mu=s2mu)
  }
  else{
    s2eta <- sum((pooling$residuals-tapply(pooling$residuals,id,mean)[as.character(id)]-tapply(pooling$residuals,time,mean)[as.character(time)])^2)/((n-1)*(T-1))
    lambda2 <- sum(tapply(pooling$residuals,id,mean)^2)*T/(n-1)
    lambda3 <- sum(tapply(pooling$residuals,time,mean)^2)*n/(T-1)
    lambda4 <- lambda2+lambda3-s2eta
    theta1 <- 1-sqrt(s2eta/lambda2)
    theta2 <- 1-sqrt(s2eta/lambda3)
    theta3 <- theta1+theta2+sqrt(s2eta/lambda4)-1
    if (theta2<0) theta2<-theta3<-0
    z <- list(theta1=theta1,theta2=theta2,theta3=theta3)
  }
  z
}

amemiya <- function(res){
  data <- attr(res,"data")
  pdim <- attr(res,"pdim")
  pmodel <- attr(res,"pmodel")
  double <- pmodel$effect=="double"
  effect <- pmodel$effect
  balanced <- pdim$balanced
  within <- res$within
  id <- data$id
  time <- data$time
  T <- pdim$nT$T
  n <- pdim$nT$n
  if(!balanced)stop("amemiya variance decomposition not implemented for unbalanced panels")
  if(!double){
    if (effect=="individual") condvar <- id else condvar <- time
    uest <- within$res+within$FE-as.vector(within$alpha)
    uest <- data$y$y-data$X$X%*%within$coef-as.vector(within$alpha)
    s21 <- T^2/(n*(T-1)-2)*sum(tapply(uest,condvar,mean)^2)
    s2eta <- sum((within$residuals-tapply(within$residuals,condvar,mean)[as.character(condvar)])^2)/(n*(T-1)-2)
    s2mu <- (s21-s2eta)/T
    theta <- 1-sqrt(s2eta/s21)
    z <- list(theta=theta,s2eta=s2eta,s21=s21,s2mu=s2mu)
  }
  else{
    y <- data$y$y
    X <- data$X$X
    uest <- as.vector(y-within$alpha-X%*%within$coef)
    s2eta <- sum((uest-tapply(uest,id,mean)[as.character(id)]-tapply(uest,time,mean)[as.character(time)])^2)/((n-1)*(T-1))
    lambda2 <- sum(tapply(uest,id,mean)^2)*T/(n-1)
    lambda3 <- sum(tapply(uest,time,mean)^2)*n/(T-1)
    lambda4 <- lambda2+lambda3-s2eta
    theta1 <- 1-sqrt(s2eta/lambda2)
    theta2 <- 1-sqrt(s2eta/lambda3)
    theta3 <- theta1+theta2+sqrt(s2eta/lambda4)-1
    if (theta2<0) theta2<-theta3<-0
    z <- list(theta1=theta1,theta2=theta2,theta3=theta3)
  }
  z
}

nerlove <- function(res){
  data <- attr(res,"data")
  pdim <- attr(res,"pdim")
  pmodel <- attr(res,"pmodel")
  double <- pmodel$effect=="double"
  effect <- pmodel$effect
  balanced <- pdim$balanced
  within <- res$within
  n <- pdim$nT$n
  N <- pdim$nT$N
  if(!double & balanced){
      s2eta <- sum(within$residuals^2)/N
      s2mu <- sum((within$FE-mean(within$FE))^2)/(n-1)
      s21 <- T*s2mu+s2eta
      theta <- 1-sqrt(s2eta/s21)
      z <- list(theta=theta,s2eta=s2eta,s21=s21,s2mu=s2mu)
    }
    else stop("nerlove variance decomposition only implemented for balanced oneway panels")
  z
}
