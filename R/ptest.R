phtest <- function(x,...){
  UseMethod("phtest")
}

phtest.plms <- function(x,...){
  data.name <- paste(deparse(substitute(x)))
  within <- x$within
  random <- x$random
  phtest(within,random,data.name)
}

phtest.panelmodel <- function(x,x2,data=NULL,...){
  # x is assumed to have a greater variance than x2
  coef.wi <- coefficients(x)
  coef.re <- coefficients(x2)
  vcov.wi <- vcov(x)
  vcov.re <- vcov(x2)
  names.wi <- names(coef.wi)
  names.re <- names(coef.re)
  coef.h <- names.re[names.re%in%names.wi]
  dbeta <- coef.wi[coef.h]-coef.re[coef.h]
  df <- length(dbeta)
  dvcov <- vcov.re[coef.h,coef.h]-vcov.wi[coef.h,coef.h]
  stat <- t(dbeta)%*%solve(dvcov)%*%dbeta
  pval <- (1-pchisq(stat,df=df))
  names(stat) <- "chi2"
  parameter <- df
  names(parameter) <- "df"
  if (is.null(data)){
    data.name <- paste(deparse(substitute(x)), "and",  deparse(substitute(x2)))
  }
  else{
    data.name <- data
  }

  res <- list(statistic = stat,
                p.value = pval,
                parameter = parameter,
                method = "Hausman Test",
                data.name = data.name)
  class(res) <- "htest"
  return(res)
}


plmtest <- function(x,...){
  UseMethod("plmtest")
}

plmtest.plms <- function(x,effect="id",type="bp", ...){
  data.name <- paste(deparse(substitute(x)))
  pooling <- x$pooling
  data.name <- pooling$call$data
  class(data)
  id.name <- attr(pooling,"indexes")$id
  time.name <- attr(pooling,"indexes")$time
  id <- eval(data.name)[[id.name]]
  time <- eval(data.name)[[time.name]]
  n <- attr(pooling,"pdim")$nT$n
  T <- attr(pooling,"pdim")$nT$T
  balanced <- attr(x,"pdim")$pdim$balanced
  res <- pooling$res
  plmtest(res,n=n,T=T,balanced=balanced,id=id,time=time,effect=effect,type=type,data=data.name)
}

plmtest.plm <- function(x,effect="id",type="bp", ...){
  data.name <- paste(deparse(substitute(x)))
  id.name <- attr(x,"indexes")$id
  time.name <- attr(x,"indexes")$time
  id <- eval(data.name)[[id.name]]
  time <- eval(data.name)[[time.name]]
  n <- attr(x,"pdim")$nT$n
  T <- attr(x,"pdim")$nT$T
  balanced <- attr(x,"pdim")$pdim$balanced
  res <- x$res
  plmtest(res,n=n,T=T,balanced=balanced,id=id,time=time,effect=effect,type=type,data=data.name)
}

plmtest.default <-  function(x,n=NULL,T=NULL,balanced=NULL,id=NULL,time=NULL,effect="id",type="bp",data=NULL, ...){
  if(effect=="id"){
    if(type != "honda" & type != "bp"){
      warning("type must be one of honda or bp, bp used")
      type="bp"
    }
    stat <-  sqrt(n*T/(2*(T-1)))*(crossprod(tapply(x,id,mean))*T^2/sum(x^2)-1)
    stat <- switch(type,honda=stat,bp=stat^2)
    names(stat) <- switch(type,honda="normal",bp="chi2")
    parameter <- switch(type,honda=NULL,bp=1)
    pval <- switch(type,honda=(1-pnorm(abs(stat)))/2,bp=1-pchisq(stat,df=1))
  }
  if(effect=="time"){
    if(type != "honda" & type != "bp"){
      warning("type must be one of honda or bp, bp used")
      type="bp"
    }
    stat <- sqrt(n*T/(2*(n-1)))*(crossprod(tapply(x,time,mean))*n^2/sum(x^2)-1)
    stat <- switch(type,honda=stat,bp=stat^2)
    names(stat) <- switch(type,honda="normal",bp="chi2")
    parameter <- switch(type,honda=NULL,bp=1)
    pval <- switch(type,honda=(1-pnorm(abs(stat)))/2,bp=1-pchisq(stat,df=1))
  }
  if(effect=="twoways"){
    if(type != "honda" & type != "bp" & type !="ghm" & type != "kw"){
      warning("type must be one of honda or bp, bp used")
      type="bp"
    }
    stat1 <-  sqrt(n*T/(2*(T-1)))*(crossprod(tapply(x,id,mean))*T^2/sum(x^2)-1)
    stat2 <-  sqrt(n*T/(2*(n-1)))*(crossprod(tapply(x,time,mean))*n^2/sum(x^2)-1)
    stat <- switch(type,
                   ghm=max(0,stat1)^2+max(0,stat2)^2,
                   bp=stat1^2+stat2^2,
                   honda=(stat1+stat2)/sqrt(2),
                   kw=sqrt((T-1)/(n+T-2))*stat1+sqrt((n-1)/(n+T-2))*stat2)
    parameter <- 2
    names(stat) <- switch(type,ghm="chi2",honda="normal",bp="chi2",kw="normal")
    pval <- switch(type,ghm=1-pchisq(stat,df=2),honda=(1-pnorm(abs(stat)))/2,bp=1-pchisq(stat,df=2),kw=(1-pnorm(abs(stat))))
  }

  method.type <- switch(type,
                   honda="Honda",
                   bp="Breush-Pagan",
                   ghm="Gourierroux, Holly and Monfort",
                   kw="King and Wu")
  method.effect <- switch(effect,
                          id="individual effects",
                          time="time effects",
                          twoways="two-ways effects")
  method=paste("Lagrange Multiplier Test - ",method.effect," (",method.type,")\n",sep="")
    
  if(type=="honda"){
    res <- list(statistic = stat,
                p.value = pval,
                method = method,
                data.name = data)
  }
  else{
    names(parameter) <- "df"
    res <- list(statistic = stat,
                p.value = pval,
                method = method,
                parameter=parameter,
                data.name = data)
  }
  class(res) <- "htest"
  res
}

pFtest <- function(x,...){
  UseMethod("pFtest")
}

pFtest.plms <- function(x,...){
  data.name <- paste(deparse(substitute(x)))
  within <- x$within
  pooling <- x$pooling
  pFtest(within,pooling,data.name)
}

pFtest.plm <- function(x,z,data=NULL, ...){
  if (is.null(data)){
    data.name <- paste(deparse(substitute(x)), "and",  deparse(substitute(z)))
  }
  else{
    data.name <- data
  }
  within <- x
  pooling <- z
  df1 <- df.residual(pooling)-df.residual(within)
  df2 <- df.residual(within)
  ssrp <- sum(residuals(pooling)^2)
  ssrw <- sum(residuals(within)^2)
  stat <- (ssrp-ssrw)/ssrw/df1*df2
  names(stat) <- "F"
  parameter <- c(df1,df2)
  names(parameter) <- c("df1","df2")
  pval <- 1-pf(stat,df1,df2)
  res <- list(statistic = stat,
              p.value = pval,
              method = "F test for effects",
              parameter=parameter,
              data.name=data.name)
  class(res) <- "htest"
  res
}
  
Ftest <- function(x,...){
  pdim <- attr(x,"pdim")
  model.name <- attr(x,"pmodel")$model
  df1 <- x$df.residual
  df2 <- switch(model.name,
                "within"=pdim$Kw,
                "between"=pdim$Kb,
                "pooling"=pdim$K,
                "random"=pdim$K,
                "ht"=pdim$K
                )
  stat <- (x$tss-x$ssr)/x$ssr*df1/df2
  names(stat) <- "F"
  pval <- 1-pf(stat,df1,df2)
  parameter <- c(df1,df2)
  names(parameter) <- c("df1","df2")
  res <- list(statistic = stat,
              parameter = parameter,
              p.value = pval,
              method = "F test",
              data.name = "data.name")
  class(res) <- "htest"
  res
}

pwaldtest <- function(x,...){
  pdim <- attr(x,"pdim")
  df <- switch(x$model.name,
                "within"=pdim$Kw,
                "between"=pdim$Kb,
                "pooling"=pdim$K,
                "random"=pdim$K,
                "ht"=pdim$K
                )
  if (names(coefficients(x))[1]=="(intercept)"){
    coef <- coefficients(x)[-1]
    vcv <- vcov(x)[-1,-1]
  }
  else{
    coef <- coefficients(x)
    vcv <- vcov(x)
  }
  parameter <- length(coef)
  stat <- coef%*%solve(vcv)%*%coef
  names(stat) <- "chisq"
  names(parameter) <- "df"
  pval <- 1-pchisq(stat,df)
  res <- list(statistic = stat,
              p.value = pval,
              parameter = parameter,
              parameter.name = "df",
              method = "Wald Test",
              data.name = "data.name")
  class(res) <- "htest"
  res
}

pooltest <- function(x,...){
  UseMethod("pooltest")
}


pooltest.plms <- function(x,nopool,effect=FALSE,...){
  if(effect){
    z <- x$within
  }
  else{
    z <- x$pooling
  }
  pooltest(z,nopool)
}

pooltest.plm <- function(x,nopool,...){
  rss <- sum(residuals(x)^2)
  uss <- sum(unlist(residuals(nopool))^2)
  dlr <- df.residual(x)
  dlu <- df.residual(nopool)
  df1 <- dlr-dlu
  df2 <- dlu
  stat <- (rss-uss)/uss*df2/df1
  pval <- 1-pf(stat,df1,df2)
  parameter <- c(df1,df2)
  names(parameter) <- c("df1","df2")
  names(stat)="F"
  data.name <- paste(deparse(substitute(plms)))
  res <- list(statistic = stat,
              parameter=parameter,
              p.value = pval,
              data.name=data.name,
              null.value = "stability",
              method = "F statistic")
  class(res) <- "htest"
  res
}
