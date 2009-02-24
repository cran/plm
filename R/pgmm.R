pgmm.formula <- function(formula,data,effect,model,instruments,gmm.inst,lag.gmm,transformation){
  if (is.null(instruments)){
    if (!is.list(gmm.inst)){
      var.gmm <- attr(terms(gmm.inst),"term.labels")
      J <- length(var.gmm)
      var.tot <- attr(formula,"var")
      lag.tot <- attr(formula,"lag")
      log.tot <- attr(formula,"log")
      diff.tot <- attr(formula,"diff")
      var.inst <- var.tot[! var.tot %in% var.gmm]
      if (length(var.inst)>0){
        lag.inst <- lag.tot[!var.tot %in% var.gmm]
        log.inst <- log.tot[!var.tot %in% var.gmm]
        diff.inst <- diff.tot[!var.tot %in% var.gmm]
        form.inst <- as.formula(paste("~",paste(var.inst,collapse="+"),sep=""))
        instruments <- dynformula(form.inst,lag.inst,diff.inst,log.inst)
      }
    }
  }
  if (!is.list(lag.gmm)){
      var.gmm <- attr(terms(gmm.inst),"term.labels")
      J <- length(var.gmm)
    lag.gmm <- rep(list(lag.gmm),J)
  }
  max.lag.gmm <- max(sapply(lag.gmm,function(x) x[1]))
  max.lag.model <- max(sapply(attr(formula,"lag"),max))+1
  time.lost <- max(max.lag.model,max.lag.gmm)
  
  gmm <- list(gmm.inst=gmm.inst,lag.gmm=lag.gmm,instruments=instruments)
  if (transformation=="ld"){
    lag.gmm <- rep(list(c(1,1)),J)
    gmm.inst.level <- dynformula(gmm.inst,diff=T)
    gmm.level <- list(gmm.inst=gmm.inst.level,lag.gmm=lag.gmm,instruments=instruments)
    z <- list(gmm.diff=gmm,gmm.level=gmm.level,time.lost=time.lost)
  }
  else{
    z <- list(gmm=gmm,time.lost=time.lost)
  }
  z
}

extract.data <- function(formula,data,time.name,id.name){
  if (length(formula)==3){
    formula <- update(formula,.~.-1)
    exo.name <- deparse(formula[[2]])
  }
  else{
    formula <- update(formula,~.-1)
  }
  indexes <- attr(data,"indexes")
  yX <- model.frame(formula,data,na.action=NULL)
  yX$time <- data[[time.name]]
  yX <- split(yX,data[[id.name]])
  yX <- lapply(yX,function(x){ rownames(x) <- x[["time"]];return(x)})
  if (length(formula)==3){
    yX <- lapply(yX,function(x) cbind(x[,1],model.matrix(formula,x)))
    yX <- lapply(yX,function(x){colnames(x)[1] <- exo.name;x})
  }
  else{
    yX <- lapply(yX,function(x) model.matrix(formula,x))
  }    
}

pgmm <- function(formula, data, effect = "individual",
                 model = "twosteps",
                 instruments = NULL, gmm.inst, lag.gmm, transformation = "d", fsm = NULL, index = NULL, ...){

  model.name <- model
  data.name <- paste(deparse(substitute(data)))
  new.data.name <- "mydata"
  data2 <- data2plm.data(data,index)
  data <- data2$data
  id.name <- data2$id.name
  time.name <- data2$time.name
  for (i in 1:length(data)){
    attr(data[[i]],"data") <- new.data.name
    attr(data[[i]],"class") <- c("pserie",attr(data[[i]],"class"))
  }
  indexes <- list(id=id.name,time=time.name)
  class(indexes) <- "indexes"
  pdim <- pdim(data)
  data <- structure(data,pdim=pdim,indexes=indexes)
  nframe <- length(sys.calls())
  assign(new.data.name,data,env=sys.frame(which=nframe))
  if(is.null(fsm)){
    fsm <- switch(transformation,
                  "d"="G",
                  "ld"="full"
                  )
  }
  cl <- match.call()
  cl$formula <- formula
  cl$gmm.inst <- gmm.inst

  if (is.null(effect)) effect <- "none"
  time.names <- pdim$panel.names$time.names
  id.names <- pdim$panel.names$id.names
  z <- pgmm.formula(formula,data,effect,model,instruments,gmm.inst,lag.gmm,transformation)

  pmodel <- list(formula=formula,effect=effect,instruments=instruments,
                 model.name=model,transformation=transformation,time.lost=z$time.lost)
  if (transformation=="ld"){
    gmm.desc <- z$gmm.diff
    gmm.desc.level <- z$gmm.level
  }
  else{
    gmm.desc <- z$gmm
  }
  time.lost <- z$time.lost
  data.name <- deparse(substitute(data))

  T <- pdim$nT$T
  ti <- split(data[[time.name]],data[[id.name]])

  yX <- extract.data(formula,data,time.name,id.name)
  yX <- lapply(yX,function(x) if(time.lost==1) x else x[-c(1:(time.lost-1)),])
  W <- extract.data(gmm.desc$gmm.inst,data,time.name,id.name)
  J <- makeJ(time.names,gmm.desc$gmm.inst,gmm.desc$lag.gmm,time.lost)
  W <- lapply(W,momatrix,J,time.names)
  if (is.null(gmm.desc$instruments)){
    In <- NULL
  }
  else{
    In <- extract.data(gmm.desc$instruments,data,time.name,id.name)
    In <- lapply(In,function(x) if(time.lost==1) x else x[-c(1:(time.lost-1)),])
  }
  if (transformation=="ld"){
    Wl <- extract.data(gmm.desc.level$gmm.inst,data,time.name,id.name)
    Wl <- lapply(Wl,function(x) x[-1,,drop=FALSE])
    Jl <- makeJ(time.names[-1],gmm.desc.level$gmm.inst,gmm.desc.level$lag.gmm,time.lost-1)
    Wl <- lapply(Wl,momatrix,Jl,time.names[-1])
    Wl <- lapply(Wl,function(x){prems <- which(time.names==rownames(x)[1]);z <- rbind(0,x);rownames(z)[1] <- time.names[prems-1];z})
  }
  else{
    Wl=NULL
  }
  if (effect!="individual"){
    if (transformation=="ld"){
      time.dummies <- cbind(1,diag(1,T)[,-(1:(time.lost))])
      dimnames(time.dummies) <- list(time.names,c("(intercept)",time.names[(time.lost+1):T]))
    }
    else{
      time.dummies <- diag(1,T)[,-(1:(time.lost))]
      dimnames(time.dummies) <- list(time.names,time.names[(time.lost+1):T])
    }
  }
  else{
    time.dummies <- NULL
  }

  if (transformation=="ld"){
    result <- pgmm.sys(yX,W,Wl,In,time.dummies,model,fsm,cl)
  }
  else{
    result <- pgmm.diff(yX,W,In,time.dummies,model,fsm,cl)
  }
  structure(result,class=c("pgmm","panelmodel"),pdim=pdim,pmodel=pmodel)
}


pgmm.diff <- function(yX,W,In,time.dummies,model,fsm,cl){
  if(!is.null(time.dummies)){
    yX <- lapply(yX,function(x) cbind(x,time.dummies[rownames(x),]))
    W <- lapply(W,function(x) cbind(x,time.dummies[rownames(x),]))
  }
  yX <- lapply(yX,diff)
  if (!is.null(In)){
    In <- lapply(In,diff)
    W <- mapply(cbind,W,In,SIMPLIFY=FALSE)
  }
  Vi <- lapply(W,function(x) crossprod(t(crossprod(x,FSM(dim(x)[1],fsm))),x))
  A1 <- solve(suml(Vi))*length(W)
  WyXi <- mapply(crossprod,W,yX,SIMPLIFY=FALSE)
  Wyi <- lapply(WyXi,function(x) x[,1])
  WXi <- lapply(WyXi,function(x) x[,-1])
  Wy <- suml(Wyi)
  WX <- suml(WXi)
  var.names <- colnames(yX[[1]])
  B1 <- solve(t(WX)%*%A1%*%WX)
  rownames(B1) <- colnames(B1) <- var.names[-1]
  coefficients <- B1%*%(t(WX)%*%A1%*%Wy)
  dim(coefficients) <- NULL
  names(coefficients) <- var.names[-1]
  residuals <- lapply(yX,function(x) as.vector(x[,1]-crossprod(t(x[,-1]),coefficients)))
  outresid <- lapply(residuals,function(x) outer(x,x))
  A2 <- mapply(crossprod,W,outresid,SIMPLIFY=FALSE)
  A2 <- mapply("%*%",A2,W,SIMPLIFY=FALSE)
  A2 <- solve(suml(A2))
  B2 <- solve(t(WX)%*%A2%*%WX)
  rownames(B2) <- colnames(B2) <- var.names[-1]
  if (model=="twosteps"){
    coef1s <- coefficients
    coefficients <- B2%*%(t(WX)%*%A2%*%Wy)
    dim(coefficients) <- NULL
    names(coefficients) <- var.names[-1]
    vcov <- B2
  }
  else{
    vcov <- B1
  }
  residuals <- lapply(yX,function(x){
    nz <- rownames(x)
    z <- as.vector(x[,1]-crossprod(t(x[,-1]),coefficients))
    names(z) <- nz
    z
  }
                      )
  fitted.values <- mapply(function(x,y) x[,1]-y,yX,residuals)
  n <- apply(sapply(yX,dim),1,sum)[1]
  K <- length(attr(terms(as.formula(cl$formula)),"term.labels"))
  Kt <- dim(yX[[1]])[2]-K-1
  p <- ncol(W[[1]])
  Ky <- attr(as.formula(cl$formula),"lag")[[1]][2]
  if(is.na(Ky)) Ky <- 0
  K <- list(K=K-Ky,Ky=Ky,Kt=Kt)
  if (model=="twosteps") coefficients <- list(coef1s,coefficients)

  list(coefficients=coefficients,residuals=residuals,vcov=vcov,
       fitted.values=fitted.values,
       df.residual=df.residual,
       model=yX,W=W,K=K,A1=A1,A2=A2,call=cl)
}




pgmm.sys <- function(yX,W,Wl,In,time.dummies,model,fsm,cl){
  if(!is.null(time.dummies)){
    yX <- lapply(yX,function(x) cbind(x,time.dummies[rownames(x),]))
    Wl <- lapply(Wl,function(x) cbind(x,time.dummies[rownames(x),]))
  }
  else{
    yX <- lapply(yX,function(x){x <- cbind(x,1);colnames(x)[dim(x)[2]] <- "(intercept)";x})
    Wl <- lapply(Wl,function(x){x <- cbind(x,1);colnames(x)[dim(x)[2]] <- "(intercept)";x})
  }
    
  if (!is.null(In)){
    Inl <- In
    In <- lapply(In,diff)
    W <- mapply(cbind,W,In,SIMPLIFY=FALSE)
    Wl <- mapply(cbind,Wl,Inl,SIMPLIFY=FALSE)
  }
  var.names <- colnames(yX[[1]])
  yXl <- yX
  yX <- lapply(yX,diff)
  pi <- lapply(Wl,nrow)
  F <- lapply(pi,FSM,fsm)
  WS <- mapply(bdiag,W,Wl,SIMPLIFY=FALSE)
  yXS <- mapply(rbind,yX,yXl,SIMPLIFY=FALSE)
  WyXi <- mapply(crossprod,WS,yXS,SIMPLIFY=FALSE)
  Wyi <- lapply(WyXi,function(x) x[,1])
  WXi <- lapply(WyXi,function(x) x[,-1])
  Wy <- suml(Wyi)
  WX <- suml(WXi)
  Vi <- mapply(function(x,y) crossprod(t(crossprod(x,y)),x),WS,F,SIMPLIFY=FALSE)
  A1 <- solve(suml(Vi))*length(WS)
  B1 <- solve(t(WX)%*%A1%*%WX)
  coefficients <- B1%*%(t(WX)%*%A1%*%Wy)
  dim(coefficients) <- NULL
  names(coefficients) <- var.names[-1]
  residuals <- lapply(yXl,function(x){nx <- rownames(x);z <- as.vector(x[,1]-crossprod(t(x[,-1]),coefficients));names(z) <- nx;z})
  resid <- lapply(residuals,function(x) c(diff(x),x))
  outresid <- lapply(resid,function(x) outer(x,x))
  Vi <- mapply(function(x,y) crossprod(t(crossprod(x,y)),x),WS,outresid,SIMPLIFY=FALSE)
  A2 <- solve(suml(Vi))
  B2 <- solve(t(WX)%*%A2%*%WX)
  vcov <- B1
  if (model=="twosteps"){
    coef1s <- coefficients
    coefficients <- B2%*%(t(WX)%*%A2%*%Wy)
    dim(coefficients) <- NULL
    names(coefficients) <- var.names[-1]
    vcov <- B2
    residuals <- lapply(yXl,function(x){nx <- rownames(x);z <- as.vector(x[,1]-crossprod(t(x[,-1]),coefficients));names(z) <- nx;z})
  }
  fitted.values <- mapply(function(x,y) x[,1]-y,yXl,residuals)
  n <- apply(sapply(yX,dim),1,sum)[1]
  K <- length(attr(terms(as.formula(cl$formula)),"term.labels"))
  Kt <- dim(yX[[1]])[2]-K-1
  p <- ncol(W[[1]])
  dim(coefficients) <- NULL
  names(coefficients) <- rownames(vcov) <- colnames(vcov) <- var.names[-1]
  Ky <- attr(as.formula(cl$formula),"lag")[[1]][2]
  if(is.na(Ky)) Ky <- 0
  K <- list(K=K-Ky,Ky=Ky,Kt=Kt)
  if (model=="twosteps") coefficients <- list(coef1s,coefficients)
  list(coefficients=coefficients,residuals=residuals,
       fitted.values=fitted.values,vcov=vcov,
       df.residual=df.residual,
       model=yXl,W=WS,K=K,A1=A1,A2=A2,call=cl,Wd=W,Wl=Wl)
}


makeJ <- function(time.names,gmminst,lag.gmm,time.lost){
  T <- length(time.names)
  names.gmm <- attr(terms(gmminst),"term.label")
  J <- array(0,dim=c(T,length(names.gmm),3),
             dimnames=list(time.names,names.gmm,
               c("start","end","n")))
  first.period <- sapply(lag.gmm,max)
  last.period <- sapply(lag.gmm,min)
  names(first.period) <- names(last.period) <- names.gmm
  for (ni in names.gmm){
    for (t in 1:T){
      J[t,ni,"start"] <- max(1,t-first.period[ni])
      J[t,ni,"end"] <- max(1,min(t-last.period[ni],T))
    }
  }
  
  J[,,"n"] <- J[,,"end",drop=FALSE]-J[,,"start",drop=FALSE]+1
  if (time.lost!=0){
    J <- J[-(1:time.lost),,,drop=FALSE]
  }
  J
}

momatrix <- function(x,J,ttot){
  names.gmm <- dimnames(J)[[2]]
  z <- matrix(0,nrow=length(ttot),ncol=length(names.gmm),dimnames=list(ttot,names.gmm))
  z[rownames(x),] <- x
  t.kept <- dimnames(J)[[1]]
  t.drop <- length(ttot)-length(t.kept)
  start <- which(ttot==rownames(x)[1])
  cnames <- c()
  for (y in t.kept){
    my <- c()
    for (ng in names.gmm){
      my <- c(my,z[seq(J[y,ng,1],J[y,ng,2]),ng])
    }
    cnames <- c(cnames,names(my))
    if (y==t.kept[1]){
      maty <- matrix(my,nrow=1)
    }
    else{
      lgn <- c(rep(0,ncol(maty)),my)
      maty <- cbind(maty,matrix(0,nrow=nrow(maty),ncol=length(my)))
      maty <- rbind(maty,lgn)
    }
  }
  rownames(maty) <- t.kept
  maty <- maty[rownames(x)[(t.drop+1):dim(x)[1]],]
  maty
}

gg <- function(name,lags,diff){
  lags <- switch(length(lags),
                 "1"=c(0,lags),
                 "2"=sort(lags),
                 stop("lags should be of length 1 or 2\n")
                 )
  lag.string <- ifelse(diff,"diff","lag")
  chlag <- c()
  if (lags[2]!=0){
    lags <- lags[1]:lags[2]
    for (i in lags){
      if (i==0){
        if (diff) chlag <- c(chlag,paste("diff(",name,")")) else chlag <- c(chlag,name)
      }
      else{
        ichar <- paste(i)
        chlag <- c(chlag,paste(lag.string,"(",name,",",i,")",sep=""))
      }
    }
    ret <- paste(chlag,collapse="+")
  }
  else{
    if (diff) chlag <- paste("diff(",name,")") else chlag <- name
    ret <- chlag
  }
  ret
}   

print.dynformula <- function(x,...){
  attr(x,"lag") <- attr(x,"var") <- attr(x,"log") <- attr(x,"diff") <- NULL
  print.formula(x)
}

G <- function(t){
  G <- matrix(0,t,t)
  for (i in 1:(t-1)){
    G[i,i] <- 2
    G[i,i+1] <- -1
    G[i+1,i] <- -1
  }
  G[t,t] <- 2
  G
}


FD <- function(t){
  FD <- Id(t)[-1,]
  for (i in 1:(t-1)){
    FD[i,i] <- -1
  }
  FD
}

Id <- function(t){
  diag(rep(1,t))
}

FSM <- function(t,fsm){
  switch(fsm,
         "I"=Id(t),
         "G"=G(t),
         "GI"=bdiag(G(t-1),diag(1,t)),
         "full"=rbind(cbind(G(t-1),FD(t)),cbind(t(FD(t)),Id(t)))
         )
}

summary.pgmm <- function(object,robust=FALSE,...){
  model.name <- attr(object,"pmodel")$model
  transformation <- attr(object,"pmodel")$transformation
  if (robust){
    vv <- vcovHC(object)
  }
  else{
    vv <- vcov(object)
  }
  rowsel <- object$K$K+object$K$Ky
  std.err <- sqrt(diag(vv))
  b <- coef(object)
  z <- b/std.err
#  p <- 2*(1-pnorm(abs(z)))
  p <- 2*pnorm(abs(z),lower.tail=FALSE)
  CoefTable <- cbind(b,std.err,z,p)
  colnames(CoefTable) <- c("Estimate","Std. Error","z-value","Pr(>|z|)")
  object$CoefTable <- CoefTable[1:rowsel,,drop=FALSE]
  object$sargan <- sargan(object)
  object$m1 <- mtest(object,1,vv)
  object$m2 <- mtest(object,2,vv)
  object$wald.coef <- wald(object,"param",vv)
  if (describe(object, "effect") == "twoways") object$wald.td <- wald(object,"time",vv)
  class(object) <- "summary.pgmm"
  object
}

print.summary.pgmm <- function(x,digits=max(3, getOption("digits") - 2), width = getOption("width"),...){
  transformation <- attr(x,"pmodel")$transformation
  pdim <- attr(x,"pdim")
  pmodel <- attr(x,"pmodel")
  effect <- pmodel$effect
  formula <- pmodel$formula
  model.name <- pmodel$model.name
  cat(paste(effect.pgmm.list[effect]," ",sep=""))
  cat(paste(model.pgmm.list[model.name],"\n",sep=""))
  cat("\nCall:\n")
  print(x$call)
  cat("\n")
  print(pdim)
  ntot <- apply(sapply(x$model,dim),1,sum)[1]
  cat("\nNumber of Observations Used: ",ntot,"\n")
  
  cat("\nResiduals\n")
  print(summary(unlist(residuals(x))))
  cat("\nCoefficients\n")
  printCoefmat(x$CoefTable,digits=digits)

  cat("\nSargan Test: ",names(x$sargan$statistic),
      "(",x$sargan$parameter,") = ",x$sargan$statistic,
      " (p.value=",format.pval(x$sargan$p.value,digits=digits),")\n",sep="")

  cat("Autocorrelation test (1): ",names(x$m1$statistic),
      " = ",x$m1$statistic,
      " (p.value=",format.pval(x$m1$p.value,digits=digits),")\n",sep="")
  
  cat("Autocorrelation test (2): ",names(x$m2$statistic),
      " = ",x$m2$statistic,
      " (p.value=",format.pval(x$m2$p.value,digits=digits),")\n",sep="")
  cat("Wald test for coefficients: ",names(x$wald.coef$statistic),
      "(",x$wald.coef$parameter,") = ",x$wald.coef$statistic,
      " (p.value=",format.pval(x$wald.coef$p.value,digits=digits),")\n",sep="")
  
  
  if (describe(x, "effect") == "twoways"){
    cat("Wald test for time dummies: ",names(x$wald.td$statistic),
        "(",x$wald.td$parameter,") = ",x$wald.td$statistic,
        " (p.value=",format.pval(x$wald.td$p.value,digits=digits),")\n",sep="")
  }
  invisible(x)
}


bdiag <- function(...){
  if (nargs() == 1)
    x <- as.list(...)
  else
    x <- list(...)
  n <- length(x)
  if(n==0) return(NULL)
  x <- lapply(x, function(y) if(length(y)) as.matrix(y) else
              stop("Zero-length component in x"))
  d <- array(unlist(lapply(x, dim)), c(2, n))
  rr <- d[1,]
  cc <- d[2,]
  rsum <- sum(rr)
  csum <- sum(cc)
  out <- array(0, c(rsum, csum))
  ind <- array(0, c(4, n))
  rcum <- cumsum(rr)
  ccum <- cumsum(cc)
  ind[1,-1] <- rcum[-n]
  ind[2,] <- rcum
  ind[3,-1] <- ccum[-n]
  ind[4,] <- ccum
  imat <- array(1:(rsum * csum), c(rsum, csum))
  iuse <- apply(ind, 2, function(y, imat) imat[(y[1]+1):y[2],
                                               (y[3]+1):y[4]], imat=imat)
  iuse <- as.vector(unlist(iuse))
  out[iuse] <- unlist(x)
  return(out)
} 


sargan <- function(object){
  transformation <- attr(object,"pmodel")$transformation
  model.name <- attr(object,"pmodel")$model
  Ktot <- object$K$K+object$K$Ky+object$K$Kt
  if (transformation=="ld"){
#    WS <- mapply(bdiag,object$W[[1]],object$W[[2]])
    resS <- lapply(object$residuals,function(x) c(diff(x),x))
    z <- suml(mapply(function(x,y) t(x)%*%y,object$W,resS,SIMPLIFY=FALSE))
    p <- ncol(object$W[[1]])
  }
  else{
    z <- suml(mapply(function(x,y) t(x)%*%y,object$W,object$residuals,SIMPLIFY=FALSE))
    p <- ncol(object$W[[1]])
  }
  if (model.name=="onestep"){
    A <- object$A1
  }
  else{
    A <- object$A2
  }
  stat <- as.numeric(crossprod(z,t(crossprod(z,A))))
  parameter <- p-Ktot
  names(parameter) <- "df"
  names(stat) <- "chisq"
  method <- "Sargan test"
#  pval <- 1-pchisq(stat,df=parameter)
  pval <- pchisq(stat,df=parameter,lower.tail=FALSE)
  sargan <- list(statistic = stat,
                 p.value = pval,
                 parameter = parameter,
                 method = "Sargan Test")
  class(sargan) <- "htest"
  sargan
}

wald <- function(x,param="coef",vcov=NULL){
  myvcov <- vcov
  if (is.null(vcov)){
    vv <- vcov(x)
  }
  else if (is.function(vcov)){
    vv <- myvcov(x)
  }
  else{
    vv <- myvcov
  }
  transformation <- attr(x,"pmodel")$transformation
  model.name <- attr(x,"pmodel")$model.name
  if (model.name=="onestep"){
    coefficients <- x$coefficients
  }
  else{
    coefficients <- x$coefficients[[2]]
  }
  Ktot <- length(coefficients)
  if (param=="time"){
    start <- switch(transformation,
                    "d"=(Ktot-x$K$Kt+1),
                    "ld"=(Ktot-x$K$Kt+2)
                    )
    end <- Ktot
  }
  else{
    end <- Ktot-x$K$Kt
    start <- 1
  }
  coef <- coefficients[start:end]
  vv <- vv[start:end,start:end]
  stat <- t(coef)%*%solve(vv)%*%coef
  names(stat) <- "chisq"
  parameter <- length(coef)
#  pval <- 1-pchisq(stat,df=parameter)
  pval <- pchisq(stat,df=parameter,lower.tail=FALSE)
  wald <- list(statistic = stat,
               p.value = pval,
               parameter = parameter,
               method = "Wald test")
  class(wald) <- "htest"
  wald
}

mtest <- function(object,order=1,vcov=NULL){
  myvcov <- vcov
  if (is.null(vcov)){
    vv <- vcov(object)
  }
  else if (is.function(vcov)){
    vv <- myvcov(object)
  }
  else{
    vv <- myvcov
  }
  time.names <- attr(object,"pdim")$panel.names$time.names
  time.lost <- attr(object,"pmodel")$time.lost
  model.name <- attr(object,"pmodel")$model
  time.names <- time.names[-(1:time.lost)]
  transformation <- attr(object,"pmodel")$transformation
  resid <- object$residuals
  if (transformation=="ld"){
    resid <- lapply(resid,diff)
    X <- lapply(object$model,function(x) rbind(diff(x[,-1]),x[,-1]))
    W <- object$W
    A2 <- object$A2
  }
  else{
    X <- lapply(object$model,function(x) x[,-1])
    W <- object$W
  }
  if (model.name=="onestep"){
    A <- object$A1
  }
  else{
    A <- object$A2
    
  }
  Eb <- lapply(resid,function(x){
    z <- rep(0,length(time.names))
    names(z) <- time.names
    z[names(x)] <- x
    return(z)
  })
  El <- lapply(resid,function(x){
    nx <- names(x)
    z <- c(rep(0,order),x[1:(length(x)-order)])
    names(z) <- nx
    z
  })
  Elb <- lapply(El,function(x){
    z <- rep(0,length(time.names))
    names(z) <- time.names
    z[names(x)] <- x
    z
  })

  if(transformation=="ld"){
    Eb <- lapply(Eb,function(x) c(x,rep(0,length(x)+1)))
    El <- lapply(El,function(x) c(x,rep(0,length(x)+1)))
    Elb <- lapply(Elb,function(x) c(x,rep(0,length(x)+1)))
    resid <- lapply(resid,function(x) c(x,rep(0,length(x)+1)))
  }
  EVE <- suml(mapply(function(x,y) t(y)%*%x%*%t(x)%*%y,Eb,Elb,SIMPLIFY=FALSE))
  EX <- suml(mapply(crossprod,El,X,SIMPLIFY=FALSE))
  XZ <- suml(mapply(crossprod,W,X,SIMPLIFY=FALSE))
  ZVE <- suml(mapply(function(x,y,z) t(x)%*%y%*%t(y)%*%z,W,resid,El,SIMPLIFY=FALSE))
                                        #  denom <- EVE-2*EX%*%object$B2%*%t(XZ)%*%A2%*%ZVE+EX%*%object$vcov%*%t(EX)
  denom <- EVE-2*EX%*%vcov(object)%*%t(XZ)%*%A%*%ZVE+EX%*%vv%*%t(EX)
  num <- suml(mapply(crossprod,Eb,Elb,SIMPLIFY=FALSE))
  stat <- num/sqrt(denom)
  names(stat) <- "normal"
#  pval <- 1-pnorm(abs(stat))
  pval <- pnorm(abs(stat),lower.tail=FALSE)
  mtest <- list(statistic = stat,
                p.value = pval,
                method = paste("Autocorrelation test of degree",order))
  class(mtest) <- "htest"
  mtest
}



coef.pgmm <- function(object,...){
  model.name <- attr(object,"pmodel")$model
  if(model.name=="onestep"){
    coefficients <- object$coefficients
  }
  else{
    coefficients <- object$coefficients[[2]]
  }
  coefficients
}


dynformula <- function(formula,lag.form=NULL,diff.form=NULL,log.form=NULL){
  endog <- attr(terms(formula),"term.labels")
  is.int <- attr(terms(formula),"intercept")
  if(length(formula)==3){
    exo <- deparse(formula[[2]])
    lhs <- TRUE
  }
  else{
    exo <- NULL
    lhs <- FALSE
  }
  
  K <- length(endog)

  if (is.null(lag.form)){
    lag.form <- rep(list(0),K+lhs)
  }
  else{
    if (!is.list(lag.form)){
      lag.form <- list(lag.form)
    }
    if (!is.null(names(lag.form))){
      nam <- names(lag.form)
      olag.form <- lag.form

      unnamed <- nam%in%""
      sum.unnamed <- sum(unnamed)
      if (sum.unnamed>0){
        if (sum.unnamed!=1) stop("Only one unnamed element is adminited\n")
        else default <- lag.form[[which(unnamed)]]
        lag.form <- rep(list(default),K+lhs)
        names(lag.form) <- c(exo,endog)
        lag.form[nam] <- olag.form
      }
      else{
        lag.form <- rep(list(0),K+lhs)
        names(lag.form) <- c(exo,endog)
        lag.form[nam] <- olag.form
      }
    }
    else{
      if (length(lag.form)==1){
        lag.form <- rep(lag.form,c(K+lhs))
      }
      else if (!length(lag.form) %in% c(K,K+lhs)) stop("irrelevant length for lag.form\n")
    }
  }
  if (is.null(diff.form)){
    diff.form <- rep(list(FALSE),K+lhs)
  }
  else{
    if (!is.list(diff.form)) diff.form <- list(diff.form)
    if (!is.null(names(diff.form))){
      nam <- names(diff.form)
      odiff.form <- diff.form
      unnamed <- nam%in%""
      sum.unnamed <- sum(unnamed)
      if (sum.unnamed>0){
        if (sum.unnamed!=1) stop("Only one unnamed element is adminited\n")
        else default <- diff.form[[which(unnamed)]]
        diff.form <- rep(list(default),K+lhs)
        names(diff.form) <- c(exo,endog)
        diff.form[nam] <- odiff.form
      }
      else{
        diff.form <- rep(list(FALSE),K+lhs)
        names(diff.form) <- c(exo,endog)
        diff.form[nam] <- odiff.form
      }
    }
    else{
      if (length(diff.form)==1){
        diff.form <- rep(diff.form,c(K+lhs))
      }
      else if (!length(diff.form) %in% c(K,K+lhs)) stop("irrelevant length for diff.form\n")
    }
  }

  if (is.null(log.form)){
    log.form <- rep(list(FALSE),K+lhs)
  }
  else{
    if (!is.list(log.form)) log.form <- list(log.form)
    if (!is.null(names(log.form))){
      nam <- names(log.form)
      olog.form <- log.form
      unnamed <- nam%in%""
      sum.unnamed <- sum(unnamed)
      if (sum.unnamed>0){
        if (sum.unnamed!=1) stop("Only one unnamed element is adminited\n")
        else default <- log.form[[which(unnamed)]]
        log.form <- rep(list(default),K+lhs)
        names(log.form) <- c(exo,endog)
        log.form[nam] <- olog.form
      }
      else{
        log.form <- rep(list(FALSE),K+lhs)
        names(log.form) <- c(exo,endog)
        log.form[nam] <- olog.form
      }
    }
    else{
      if (length(log.form)==1){
        log.form <- rep(log.form,c(K+lhs))
      }
      else if (!length(log.form) %in% c(K,K+lhs)) stop("irrelevant length for log.form\n")
    }
  }

  chendog <- c()
  if (lhs){
    if (log.form[[1]]) exo <- paste("log(",exo,")",sep="")
    if (length(lag.form)==K){
      lag.form <- c(0,lag.form)
    }
    if (length(lag.form[[1]])==1 && lag.form[[1]]!=0){
      lag.form[[1]] <- c(1,lag.form[[1]])
    }
    if (length(lag.form[[1]])!=1){
      chendog <- c(chendog,gg(exo,lag.form[[1]],diff.form[[1]]))
    }
  }

  j <- 1*lhs
  for (i in endog){
    j <- j+1
    if (log.form[[j]]) i <- paste("log(",i,")",sep="")
    chendog <- c(chendog,gg(i,lag.form[[j]],diff.form[[j]]))
  }
  chendog <- paste(chendog,collapse="+")
  if (!is.null(exo)){
    if (diff.form[[1]]) exo <- paste("diff(",exo,")",sep="")
    formod <- as.formula(paste(exo,"~",chendog,sep=""))
  }
  else{
    formod <- as.formula(paste("~",chendog,sep=""))
  }
  if (is.int==0) formod <- update(formod,.~.-1)
  structure(formod,class=c("dynformula","formula"),lag=lag.form,
            diff=diff.form,log=log.form,var=c(exo,endog))
}
