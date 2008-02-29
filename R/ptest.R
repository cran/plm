phtest <- function(x,...){
  UseMethod("phtest")
}

phtest.formula <- function(x, data, model = c("within","random"), effect = "individual", index=NULL, ...){
  if(length(model)!=2) stop("two models should be indicated")
  
  for (i in 1:2){
    model.name <- model[i]
    if(!(model.name %in% names(model.plm.list))){
      stop("model must be one of ",oneof(model.plm.list))
    }
  }
  form <- x
  x <- plm(form,data,model=model[1],effect=effect,index=index)
  x2 <- plm(form,data,model=model[2],effect=effect,index=index)
  phtest(x,x2)
}

phtest.panelmodel <- function(x, x2, ...){
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
  stat <- abs(t(dbeta)%*%solve(dvcov)%*%dbeta)
#  pval <- (1-pchisq(stat,df=df))
  pval <- pchisq(stat,df=df,lower.tail=FALSE)
  names(stat) <- "chisq"
  parameter <- df
  names(parameter) <- "df"
#  data.name <- paste(deparse(substitute(x)), "and",  deparse(substitute(x2)))
  data.name <- paste(deparse(x$call$formula))
  alternative <- "one model is inconsistent"
#  null.value <- "both models are consistent"
  res <- list(statistic = stat,
              p.value = pval,
              parameter = parameter,
              method = "Hausman Test",
              data.name = data.name,
 #             null.value=null.value,
              alternative=alternative)
  class(res) <- "htest"
  return(res)
}


plmtest <- function(x,...){
  UseMethod("plmtest")
}

plmtest.plm <- function(x, effect = "individual", type = "honda", ...){
  data.name <- paste(deparse(x$call$formula))
  if (is.null(x$call$model) || x$call$model!="pooling") stop("the argument must be a pooling model")
  id <- x$index$id
  time <- x$index$time
  pdim <- attr(x,"pdim")
  
  n <- pdim$nT$n
  T <- pdim$nT$T
  balanced <- pdim$balanced
  res <- resid(x)
  plmtest(res,n=n,T=T,balanced=balanced,
          id=id,time=time,effect=effect,
          type=type,data=data.name)
}

plmtest.formula <- function(x, data, effect = "individual", type = "honda", index=NULL, ...){
  data.name <- paste(deparse(x))
  data <- plm.data(data,index=index)
  x <- plm(x,data,model="pooling")

  id <- x$index$id
  time <- x$index$time
  pdim <- pdim(data)
  
  n <- pdim$nT$n
  T <- pdim$nT$T
  balanced <- pdim$balanced
  res <- resid(x)
  plmtest(res,n=n,T=T,balanced=balanced,
          id=id,time=time,effect=effect,
          type=type,data=data.name)
}


plmtest.default <-  function(x, n = NULL, T = NULL, balanced = NULL,
                             id = NULL, time = NULL, effect = "individual",
                             type= "bp" , data = NULL, ...){
  if(effect=="individual"){
    if(type != "honda" & type != "bp"){
      warning("type must be one of honda or bp, bp used")
      type="bp"
    }
    stat <-  sqrt(n*T/(2*(T-1)))*(crossprod(tapply(x,id,mean))*T^2/sum(x^2)-1)
    stat <- switch(type,honda=stat,bp=stat^2)
    names(stat) <- switch(type,honda="normal",bp="chisq")
    parameter <- switch(type,honda=NULL,bp=1)
#    pval <- switch(type,honda=(1-pnorm(abs(stat)))/2,bp=1-pchisq(stat,df=1))
    pval <- switch(type,
                   honda=pnorm(abs(stat),lower.tail=FALSE)/2,
                   bp=pchisq(stat,df=1,lower.tail=FALSE))
  }
  if(effect=="time"){
    if(type != "honda" & type != "bp"){
      warning("type must be one of honda or bp, bp used")
      type="bp"
    }
    stat <- sqrt(n*T/(2*(n-1)))*(crossprod(tapply(x,time,mean))*n^2/sum(x^2)-1)
    stat <- switch(type,honda=stat,bp=stat^2)
    names(stat) <- switch(type,honda="normal",bp="chisq")
    parameter <- switch(type,honda=NULL,bp=1)
#    pval <- switch(type,honda=(1-pnorm(abs(stat)))/2,bp=1-pchisq(stat,df=1))
    pval <- switch(type,
                   honda=pnorm(abs(stat),lower.tail=FALSE)/2,
                   bp=pchisq(stat,df=1,lower.tail=FALSE))
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
    names(stat) <- switch(type,ghm="chisq",honda="normal",bp="chisq",kw="normal")
#    pval <- switch(type,ghm=1-pchisq(stat,df=2),honda=(1-pnorm(abs(stat)))/2,bp=1-pchisq(stat,df=2),kw=(1-pnorm(abs(stat))))
    pval <- switch(type,ghm=pchisq(stat,df=2,lower.tail=FALSE),
                   honda=pnorm(abs(stat),lower.tail=FALSE)/2,
                   bp=pchisq(stat,df=2,lower.tail=FALSE),
                   kw=pnorm(abs(stat),lower.tail=FALSE))
  }

  method.type <- switch(type,
                   honda="Honda",
                   bp="Breusch-Pagan",
                   ghm="Gourieroux, Holly and Monfort",
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
#  res$null.value <- "no effects"
  res$alternative <- "significant effects"
  class(res) <- "htest"
  res
}

pFtest <- function(x,...){
  UseMethod("pFtest")
}

pFtest.formula <- function(x, data, effect = "individual", index=NULL,...){
  form <- x
  x <- plm(form,data,effect=effect,model="within",index=index)
  z <- plm(form,data,effect=effect,model="pooling",index=index)
  pFtest(x,z,...)
}
  
  

pFtest.plm <- function(x, z, ...){
#  data.name <- paste(deparse(substitute(x)), "and",  deparse(substitute(z)))
  data.name <-  paste(deparse(x$call$formula))
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
#  pval <- 1-pf(stat,df1,df2)
  pval <- pf(stat,df1,df2,lower.tail=FALSE)
#  null.value <- "no effects"
  alternative <- "significant effects"
  res <- list(statistic = stat,
              p.value = pval,
              method = "F test for effects",
              parameter=parameter,
              data.name=data.name,
#              null.value=null.value,
              alternative=alternative)
  class(res) <- "htest"
  res
}
  
Ftest <- function(x, ...){
  pdim <- attr(x,"pdim")
  model.name <- attr(x,"pmodel")$model
  data.name <- x$call$formula
  df1 <- x$df.residual
  df2 <- switch(model.name,
                "within"=pdim$Kw,
                "between"=pdim$Kb,
                "pooling"=pdim$K,
                "random"=pdim$K,
                "ht"=pdim$K,
                "fd"=pdim$K
                )
  stat <- (x$tss-x$ssr)/x$ssr*df1/df2
  names(stat) <- "F"
  pval <- pf(stat,df1,df2,lower.tail=FALSE)
  parameter <- c(df1,df2)
  names(parameter) <- c("df1","df2")
  res <- list(statistic = stat,
              parameter = parameter,
              p.value = pval,
              method = "F test"
              )
  class(res) <- "htest"
  res
}

pwaldtest <- function(x, ...){
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
#  pval <- 1-pchisq(stat,df)
  pval <- pchisq(stat,df,lower.tail=FALSE)
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

pooltest.formula <- function(x, data, effect = "individual", model = "within", index=NULL, ...){
  plm.model <- plm(x,data,effect=effect,model=model,index=index)
  pvcm.model <- pvcm(x,data,effect=effect,model="within",index=index)
  print(summary(plm.model))
  print(summary(pvcm.model))
  pooltest(plm.model,pvcm.model)
}
  

pooltest.plm <- function(x, z, ...){
  data.name <- paste(deparse(x$call$formula))
  rss <- sum(residuals(x)^2)
  uss <- sum(unlist(residuals(z))^2)
  dlr <- df.residual(x)
  dlu <- df.residual(z)
  df1 <- dlr-dlu
  df2 <- dlu
  stat <- (rss-uss)/uss*df2/df1
#  pval <- 1-pf(stat,df1,df2)
  pval <- pf(stat,df1,df2,lower.tail=FALSE)
  parameter <- c(df1,df2)
  names(parameter) <- c("df1","df2")
  names(stat)="F"
#  x.name <- paste(deparse(substitute(x)))
#  z.name <- paste(deparse(substitute(z)))
#  data.name <- paste(x.name,z.name,sep=" and ",collapse="")
  res <- list(statistic = stat,
              parameter=parameter,
              p.value = pval,
              data.name=data.name,
#              null.value = "stability",
              alternative = "unstability",
              method = "F statistic")
  class(res) <- "htest"
  res
}
