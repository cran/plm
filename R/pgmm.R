X3D <- function(form,data){
  tf <- terms(form)
  if (attr(tf,"intercept")==1){
    if (attr(tf,"response")==1){
      form <- update(form,.~.-1)
    }
    else{
      form <- update(form,~.-1)
    }
  }
  pdim <- attr(data,"pdim")
  id.names <- pdim$panel.names$id.names
  time.names <- pdim$panel.names$time.names
  n <- pdim$nT$n
  T <- pdim$nT$T
  id.name <- attr(data,"indexes")$id
  time.name <- attr(data,"indexes")$time
  mf <- model.frame(form,data,na.action=NULL)
  X <- as.data.frame(model.matrix(form,mf))
  rselect <- rownames(X)
  indexes <- data[rselect,c(id.name,time.name)]
  K <- ncol(X)
  tot <- as.vector(t(outer(id.names,time.names,paste)))
  totg <- paste(indexes[,1],indexes[,2])
  rownames(X) <- totg
  X <- X[tot,,drop=F]
  rownames(X) <- tot
  array(as.vector(as.matrix(X)),
        dim=c(T,n,K),
        dimnames=(list(time.names,id.names,colnames(X))))
}

pgmm.extract <- function(formula,data,effect,instruments,
                 inst.transformation,
                 lags.endog,first.period,last.period,
                 model,verbal){

  indexes <- attr(data,"indexes")
  id.index.name <- indexes$id
  time.index.name <- indexes$time

  if (is.character(formula)){
    pure.ar.model <- TRUE
    fortot <-  as.formula(paste(formula,"~",id.index.name,"+",time.index.name))
    y.name <- formula
    x.formula <- X <- names.X <- NULL
    K <- 0
  }
  else{
    pure.ar.model <- FALSE
    fortot <-  as.formula(paste(formula,"+",id.index.name,"+",time.index.name))
    y.name <- paste(deparse(formula[[2]]))
    x.formula <- paste(deparse(formula[[3]]),collapse="")
    x.formula <- as.formula(paste("~ ",x.formula,collapse="",sep=""))
  }
  mf <- model.frame(fortot,data=data,na.action=NULL)
  id <- mf[[id.index.name]]
  time <- mf[[time.index.name]]
  first.diff.model <- ifelse(effect %in% c("individual","twoways"),
                             TRUE,FALSE)
  t.start <- 1+(first.diff.model)+lags.endog
  pdim <- attr(data,"pdim")
  pvar <- attr(data,"pvar")
  time.var <- attr(data,"indexes")$time
  y.names <- switch(as.character(lags.endog),
                   "0"=NULL,
                   "1"=c(paste("lag(",y.name,")",collapse="",sep="")),
                   "2"=c(paste("lag(",y.name,")",collapse="",sep=""),
                     paste("lag(",y.name,",2)",collapse="",sep=""))
                 )
  y.formula <- as.formula(paste("~ ",y.name,collapse="",sep=""))
  if (!pure.ar.model){
    X <- X3D(x.formula,data)
    K <- dim(X)[3]
    names.X <- dimnames(X)[[3]]
  }
    
  y <- X3D(y.formula,data)
##   K <- dim(X)[3]
##   T <- dim(X)[1]
##   n <- dim(X)[2]
  T <- dim(y)[1]
  n <- dim(y)[2]
  Ky <- lags.endog

  names.time <- dimnames(y)[[1]]
  names.id <- dimnames(y)[[2]]

# Creation of a list containing available periods for each individual
  pdim <- pdim(id,time)
  if (!pdim$balanced){
    stop("pgmm not implemented for unbalanced panel data")
    ti <- list()
    for (i in pdim$panel.names$id.names){
      ti[[i]] <- time[id==i][-(1:(t.start-1))]
      ti[[i]] <- as.character(ti[[i]])
    }
  }
  
# Matrix Z of instruments (3D array)
  
  if (is.null(instruments)){
    if (lags.endog>0){
      names.inst <- c(y.name,names.X)
    }
    else{
      names.inst <- names.X
    }
    Z <- array(NA,dim=c(T,n,K+(Ky>0)),
               dimnames=list(names.time,names.id,
                 names.inst))
    if (lags.endog>0){
      Z[,,1] <- y[,,]
    }
    if (!pure.ar.model)  Z[,,((Ky>0)+1):(K+(Ky>0))] <- X
  }
  else{
    if (length(instruments)!=2 & class(instruments)!="formula"){
      stop("instruments must be a one side formula\n")
    }
    Z <- X3D(instruments,data)
    names.inst <- dimnames(Z)[[3]]
  }
  Ki <- length(names.inst)

#  inclusion of time dummies (3D array) and modification of the names of variables and instruments
  
  if (effect=="twoways"){
    t.formula <- as.formula(paste("~ ",time.var,sep="",collapse=""))
    time <- X3D(t.formula,data)
    time <- time[,,t.start:T]
    t.names <- dimnames(time)[[3]]
    Kt <- T-t.start+1
    var.names <- c(y.names,names.X,t.names)
    names.inst2 <- c(names.inst,"(intercept)")
  } 
 else{
    time <- NULL
    Kt <- 0
    t.names <- NULL
    var.names <- c(y.names,names.X)
    names.inst2 <- names.inst
  }

# Matrix of variables (3D array) the first is the dependent variable
  
  yX <- array(Z,dim=c(T,n,K+Ky+Kt+1),
              dimnames=list(names.time,names.id,
                c(y.name,var.names)))
  yX[,,1] <- y
  if (lags.endog>0){
    yX[2:T,,2] <- y[1:(T-1),,]
  }
  if (lags.endog>1){
    yX[3:T,,3] <- y[1:(T-2),,]
  }
  if (!pure.ar.model)  yX[,,(Ky+2):(K+Ky+1)] <- X
  if (effect=="twoways"){
    yX[,,(K+Ky+2):(K+Ky+Kt+1)] <- time
  }
  if (effect %in% c("individual","twoways")){
    yX[2:T,,] <- yX[2:T,,]-yX[1:(T-1),,]
  }
  yX <- yX[t.start:T,,]

# Construction of the vectors of first/last period for instruments and instruments transformation
  
  lfp <- length(first.period)
  llp <- length(last.period)
  lit <- length(inst.transformation)
  if (lfp!=Ki){
    first.period <- switch(as.character(lfp),
                           "1"=rep(first.period,Ki),
                           "2"=c(first.period[1],rep(first.period[2],Ki-1)),
                           stop("the first.period argument is irrelevant")
                           )
  }
  if (llp!=Ki){
    last.period <- switch(as.character(llp),
                           "1"=rep(last.period,Ki),
                           "2"=c(last.period[1],rep(last.period[2],Ki-1)),
                           stop("the last.period argument is irrelevant")
                           )
  }
  if (lit!=Ki){
    inst.transformation <- switch(as.character(lit),
                          "1"=rep(inst.transformation,Ki),
                          "2"=c(inst.transformation[1],rep(inst.transformation[2],Ki-1)),
                          stop("the inst.transformation argument is irrelevant")
                          )
  }
  names(first.period) <- names(last.period) <- names(inst.transformation) <- names.inst

# J is an array containing the first, last period and the number of period for each instrument and each year
# W is a list containing lists (one for each instrument) of instruments values for each period

  J <- array(0,dim=c(T,length(names.inst2),3),
             dimnames=list(names.time,names.inst2,
               c("start","end","n")))
  for (ni in names.inst2){
    for (t in 1:T){
      if (ni=="(intercept)"){
        J[t,ni,"start"] <- 1
        J[t,ni,"end"] <- 1
      }
      else{
        J[t,ni,"start"] <- max(1,t+first.period[ni])
        J[t,ni,"end"] <- max(1,min(t+last.period[ni],T))
      }
    }
  }
  J[,,"n"] <- J[,,"end"]-J[,,"start"]+1
  
  W <- list()
  for (ni in names.inst){
    wi <- list()
    for (t in t.start:T){
      andebut <- J[t,ni,"start"]
      anfin <- J[t,ni,"end"]
      if (inst.transformation[ni]=="l"){
        wi[[names.time[t]]] <- Z[andebut:anfin,,ni]
      }
      if (inst.transformation[ni]=="d"){
        wi[[names.time[t]]] <- Z[andebut:anfin,,ni]-Z[(andebut-1):(anfin-1),,ni]
      }
    }
    W[[ni]] <- wi
  }
  J <- J[t.start:T,,,drop=F]
  Jtotal <- sum(J[,,"n"])
  ra <- J[,,"n",drop=F]
  BJ <- BJ(apply(J[,,"n",drop=F],1,sum))

# omega is a matrix where each column is an individual. Each column is a vector containing all instruments for
# all periods (the W_i block diagonal matrix in vector form)

  omega <- c()
  for (an in names.time[t.start:T]){
    for (vk in names(W)){
      omega <- rbind(omega,W[[vk]][[an]])
     }
    if (effect=="twoways") omega <- rbind(omega,1)
  }
  omega[is.na(omega)] <- 0

  
  yX[is.na(yX)] <- 0
  
  time.span <- names.time[t.start:length(names.time)]
  row.names.omega <- rep(time.span,apply(J[,,"n",drop=F],1,sum))
  rownames(omega) <- row.names.omega
  
# In case of unbalanced panel, blocks of omega corresponding to missing
# periods and unavailable instruments for present periods are
# fixed to zero
  
  if (!pdim$balanced){
    for (i in pdim$panel.names$id.names){
      null.row.omega <- !(row.names.omega %in% ti[[i]])
      omega[null.row.omega,i] <- 0
      null.period.id <- !(time.span %in% ti[[i]])
      yX[null.period.id,i,] <- 0
    }
  }
  omega[is.na(omega)] <- 0
  WyX <- matrix(0,nrow=Jtotal,ncol=K+Ky+Kt+1)
  if(verbal) cat("Moments creation : loop on individuals ")
  for (i in names.id){
    WyX1i <- as.matrix(t(BJ)%*%yX[,i,])
    WyX2i <- matrix(omega[,i],ncol=K+Ky+Kt+1,nrow=Jtotal)
    WyXi <- WyX1i*WyX2i
    WyX <- WyX+WyXi
    if(verbal) cat(paste(i," "))
  }
  if (verbal) cat("\n")

  extract <- list(omega=omega,yX=yX,WyX=WyX,J=J,t.start=t.start,K=K,Ky=Ky,Kt=Kt)
  attr(extract,"pdim") <- pdim
  extract
}


pgmm.one.step <- function(omega,yX,WyX,J,K,Ky,Kt,verbal,effect){
  names.yX <- dimnames(yX)
  dim.yX <- dim(yX)
  names.id <- names.yX[[2]]
  names.time <- names.yX[[1]]
  var.names <- names.yX[[3]][-1]
  n <- dim.yX[2]
  T <- dim.yX[1]

  Jtotal <- sum(J[,,"n"])
  BJ <- BJ(apply(J[,,"n",drop=F],1,sum))
 
  if (effect %in% c("individual","twoways")){
    G=matrix(0,T,T)
    for (i in 1:(T-1)){
      G[i,i]=2
      G[i,i+1]=-1
      G[i+1,i]=-1
    }
    G[T,T]=2
  }
  else{
    G <- diag(rep(1,T),nrow=T)
    rownames(G) <- colnames(G) <- time.span
  }
  
  V <- (t(BJ)%*%G%*%BJ)

  WDW <- matrix(0,nrow=Jtotal,ncol=Jtotal)
  if(verbal) cat("First Step : loop on individuals ")
  for (i in names.id){
    omegai <- omega[,i]
    WDW1i <- outer(omegai,omegai)
    WDWi <- WDW1i*V
    WDW <- WDW+WDWi
    if(verbal) cat(paste(i," "))
  }
  if (verbal) cat("\n")
  WDWm1 <- as.matrix(solve(WDW))
  Wy <- WyX[,1,drop=F]
  WX <- WyX[,-1,drop=F]
  
  vcov <- solve(t(WX)%*%WDWm1%*%WX)
  coefficients <- vcov%*%(t(WX)%*%WDWm1%*%Wy)
  dim(coefficients) <- NULL
  names(coefficients) <- var.names
  yX2 <- yX
  for (k in var.names){
    yX2[,,k] <- coefficients[k]*yX[,,k]
  }
  residuals <- yX2[,,1]-apply(yX2[,,-1],c(1,2),sum)
  fitted.values <- yX2[,,1]-residuals
  dvW <- apply(crossprod(BJ,residuals)*omega,1,sum)

  stat <- t(dvW)%*%WDWm1%*%dvW
  dim(stat) <- NULL
  p <- nrow(omega)
  Ktot <- K+Ky+Kt
  parameter <- p-Ktot

  res1 <- list(coefficients=coefficients,
               residuals=residuals,
               fitted.values=fitted.values,
               vcov=vcov,
               df.residual=df.residual,
               model=yX,
               omega=omega,
               WDWm1=WDWm1,
               WyX=WyX,
               J=J)
}

pgmm.two.steps <- function(omega,yX,WyX,J,K,Ky,Kt,verbal,effect){
  res.one.step <- pgmm.one.step(omega,yX,WyX,J,K,Ky,Kt,verbal,effect)
  names.yX <- dimnames(yX)
  dim.yX <- dim(yX)
  names.id <- names.yX[[2]]
  names.time <- names.yX[[1]]
  var.names <- names.yX[[3]][-1]
  n <- dim.yX[2]
  T <- dim.yX[1]
  Jtotal <- sum(J[,,"n"])
  BJ <- BJ(apply(J[,,"n",drop=F],1,sum))
  WDW <- matrix(0,nrow=Jtotal,ncol=Jtotal)

  if(verbal) cat("Second Step : loop on individuals ")
  for (i in names.id){
    WDW1i <- crossprod(t(omega[,i]))
    vipvi <- crossprod(t(res.one.step$residuals[,i]))
    WDW2i <- (t(BJ)%*%vipvi%*%BJ)
    WDWi <- WDW1i*WDW2i
    WDW <- WDW+WDWi
    if(verbal) cat(paste(i," "))
  }
  
  if (verbal) cat("\n")
  Wy <- WyX[,1,drop=F]
  WX <- WyX[,-1,drop=F]

  WDWm1 <- as.matrix(solve(WDW))
  vcov <- solve(t(WX)%*%WDWm1%*%WX)
  coefficients <- vcov%*%(t(WX)%*%WDWm1%*%Wy)
  dim(coefficients) <- NULL
  names(coefficients) <- var.names
  
  yX2 <- yX
  for (k in var.names){
    yX2[,,k] <- coefficients[k]*yX[,,k]
  }
  
  residuals <- yX2[,,1]-apply(yX2[,,-1],c(1,2),sum)
  fitted.values <- yX2[,,1]-residuals
  df.residual <- n*T-K-Ky-Kt
  dvW <- apply(crossprod(BJ,residuals)*omega,1,sum)
  
  stat <- t(dvW)%*%WDWm1%*%dvW
  dim(stat) <- NULL
  p <- nrow(omega)
  Ktot <- K+Ky+Kt
  parameter <- p-Ktot
  
  res2 <- list(coefficients=coefficients,
               residuals=residuals,
               fitted.values=fitted.values,
               vcov=vcov,
               df.residual=df.residual,
               model=yX,
               omega=omega,
               WDWm1=WDWm1,
               WyX=WyX,
               J=J)
}

pgmm <- function(formula,data,effect="individual",model="twosteps",
                 instruments=NULL,inst.transformation="l",
                 lags.endog=1,first.period=-99,last.period=-1,
                 ...){
  model.name <- model
  verbal=FALSE
  if (!any(class(data) %in% "pdata.frame")){
    stop("argument data should be a pdata.frame\n")
  }
  if(!(effect %in% names(effect.pgmm.list))){
    stop(paste("effect must be one of",oneof(effect.pgmm.list)))
  }
  if(!(model.name %in% names(model.pgmm.list))){
    stop("model must be one of",oneof(model.pgmm.list))
  }
  
  pmodel <- list(formula=formula,effect=effect,instruments=instruments,
                 model.name=model)
  if (is.null(effect)) effect <- "none"
  cl <- match.call()
  z <- pgmm.extract(formula,data,effect,instruments,
                    inst.transformation,
                    lags.endog,first.period,last.period,
                    model,verbal)
  pdim <- attr(z,"pdim")
  omega <- z$omega
  yX <- z$yX
  J <- z$J
  t.start <- z$t.start
  K <- z$K
  Ky <- z$Ky
  Kt <- z$Kt
  WyX <- z$WyX

  if(model=="onestep"){
    z <- pgmm.one.step(omega,yX,WyX,J,K,Ky,Kt,verbal,effect)
  }
  if (model=="twosteps"){
    z <- pgmm.two.steps(omega,yX,WyX,J,K,Ky,Kt,verbal,effect)
  }
  z$K <- list(K=K,Ky=Ky,Kt=Kt)
  z$call <- cl
  z <- structure(z,class="pgmm",pdim=pdim,pmodel=pmodel)
  z
}

print.pgmm <- function(x,digits=5,...){
  print(x$coefficients)
}

summary.pgmm <- function(object,...){
  rowsel <- object$K$K+object$K$Ky
  vcov <- object$vcov
  std.err <- sqrt(diag(vcov))
  b <- object$coefficients
  z <- b/std.err
  p <- 2*(1-pnorm(abs(z)))
  CoefTable <- cbind(b,std.err,z,p)
  colnames(CoefTable) <- c("Estimate","Std. Error","z-value","Pr(>|z|)")
  object$CoefTable <- CoefTable[1:rowsel,,drop=F]
  object$sargan <- sargan(object)
  if (object$call$effect=="twoways") object$wald.td <- wald.td(object)
  class(object) <- "summary.pgmm"
  object
}

print.summary.pgmm <- function(x,digits=5,length.line=70,...){
  pdim <- attr(x,"pdim")
  pmodel <- attr(x,"pmodel")
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
  print(summary(as.vector(residuals(x))))

  centre("Model Description",length.line)
  print(x$CoefTable)
  centre("Specification tests",length.line)

  cat("Sargan Test                   : ",names(x$sargan$statistic),
      "(",x$sargan$parameter,") = ",x$sargan$statistic,
      " (p.value=",x$sargan$p.value,")\n",sep="")
  if (x$call$effect=="twoways"){
    cat("Wald test for time dummies    : ",names(x$wald.td$statistic),
        "(",x$wald.td$parameter,") = ",x$wald.td$statistic,
        " (p.value=",x$wald.td$p.value,")\n",sep="")
  }
  invisible(x)
}

sargan <- function(x){
  BJ <- BJ(apply(x$J[,,"n",drop=F],1,sum))
  dvW <- apply(crossprod(BJ,x$residuals)*x$omega,1,sum)
  stat <- t(dvW)%*%x$WDWm1%*%dvW
  dim(stat) <- NULL
  p <- nrow(z$omega)
  Ktot <- z$K$K+z$K$Ky+z$K$Kt
  parameter <- p-Ktot
  names(parameter) <- "df"
  names(stat) <- "chi2"
  method <- "Sargan test"
  pval <- 1-pchisq(stat,df=parameter)
  sargan <- list(statistic = stat,
                 p.value = pval,
                 parameter = parameter,
                 method = "Sargan Test")
  class(sargan) <- "htest"
  sargan
}

wald.td <- function(x){
  Ktot <- length(x$coefficients)
  coef <- x$coefficients[(Ktot-x$K$Kt+1):Ktot]
  vcov <- x$vcov[(Ktot-x$K$Kt+1):Ktot,(Ktot-x$K$Kt+1):Ktot]
  stat <- t(coef)%*%solve(vcov)%*%coef
  names(stat) <- "chi2"
  parameter <- x$K$Kt
  pval <- 1-pchisq(stat,df=parameter)
  wald <- list(statistic = stat,
               p.value = pval,
               parameter = parameter,
               method = "Wald test")
  class(wald) <- "htest"
  wald
}

BJ <- function(x){
  require(Matrix)
  MJ <- list()
  for (i in 1:length(x)){
    MJ[[i]] <- matrix(rep(1,x[i]),nrow=1)
  }
  bdiag(MJ)
}
    
