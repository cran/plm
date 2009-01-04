pvcm <-  function(formula, data,subset ,na.action, effect = c("individual","time"),
                  model = c("within","random"), index = NULL, ...){
  has.intercept <- (attr(terms(formula),'intercept') == 1)
  effect <- match.arg(effect)
  name.y <- paste(deparse(formula[[2]]))
  model.name <- match.arg(model)
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
  attr(data,"indexes") <- indexes
  nframe <- length(sys.calls())
  assign(new.data.name,data,env=sys.frame(which=nframe))

  cl <- match.call()
  mf <- match.call(expand.dots=FALSE)
  m <- match(c("formula","data","subset","weights","na.action","offset"),names(mf),0)
  mf <- mf[c(1,m)]
  mf$drop.unused.levels <- TRUE
  mf$data <- as.name(new.data.name)
  mf[[1]] <- as.name("model.frame")
  mindexes <- mf
  mindexes[["formula"]] <- formula(paste("~",id.name,"+",time.name,sep="",collapse=""))
  mf$formula <- formula
  mf <- eval(mf,sys.frame(which=nframe))
  mindexes <- eval(mindexes,sys.frame(which=nframe))
#  y <- model.response(mf,"numeric")
  int.row.names <- intersect(attr(mf,"row.names"),attr(mindexes,"row.names"))
  mf <- mf[int.row.names,]
  mindexes <- mindexes[int.row.names,]
  attr(mf,"row.names") <- attr(mindexes,"row.names") <- int.row.names
  terms.mf <- attr(mf,"terms")
  id <- mindexes[[id.name]]
  time <- mindexes[[time.name]]
  if (effect=="time"){ cond <- time ; other <- id}
  else{ cond <- id ; other <- time}
  pdim <- pdim(id,time)
  ml <- split(as.data.frame(mf),cond)
  nr <- sapply(ml,function(x) dim(x)[1])>0
  ml <- ml[nr]
  ols <- lapply(ml,
                function(x){
                  X <- model.matrix(terms.mf,x)
                  y <- model.response(x)
                  r <- lm(y~X-1)
                  nc <- colnames(model.frame(r)$X)
                  names(r$coefficients) <- nc
                  r
                }
                )
  extr <- list(ml=ml,ols=ols,id=id,time=time)
  pmodel <- list(model.name = model.name, formula = formula,
                 effect = effect, has.intercept = has.intercept)
  extr <- structure(extr, pdim = pdim, indexes = indexes,
                    data = data.name, pmodel = pmodel)
  extr$call <- cl
  result <- switch(model.name,
                   "within" = pvcm.within(extr,...),
                   "random" = pvcm.random(extr,...)
                   )
  class(result) <- c("pvcm","panelmodel")
  result$call <- cl
  result
}

pvcm.within <- function(extr,...){
  ml <- extr$ml
  ols <- extr$ols
  pdim <- attr(extr,"pdim")
  pmodel <- attr(extr,"pmodel")
  T <- pdim$nT$T
  n <- pdim$nT$n
  N <- pdim$nT$N
  interc <- pmodel$has.intercept
  if (pmodel$effect=="time") ncond <- T else ncond <- n
  coef <- as.data.frame(t(sapply(ols,coefficients)))
  residuals <- unlist(lapply(ols,residuals))
  vcov <- lapply(ols,vcov)
  std <- as.data.frame(t(sapply(vcov,function(x) sqrt(diag(x)))))
  names(coef) <- names(std) <- colnames(coef)
  ssr <- sum(residuals^2)
  y <- unlist(lapply(ml,function(x) x[,1]))
  fitted.values <- y-residuals
  tss <- tss(y)
  df.residuals <- N-ncond*ncol(coef)
  nopool <- list(coefficients = coef, residuals = residuals, fitted.values = fitted.values,
                 vcov = vcov, df.residuals = df.residuals, model = ml, std.error = std)
  nopool <- structure(nopool, class = "pvcm", pmodel = pmodel, pdim = pdim)
  nopool
}

pvcm.random <- function(extr, ...){
  interc <- attr(extr,"pmodel")$has.intercept
  ols <- extr$ols
  ml <- extr$ml
  pdim <- attr(extr,"pdim")
  pmodel <- attr(extr,"pmodel")
  T <- pdim$nT$T ; n <- pdim$nT$n ; N <- pdim$nT$N
  if (pmodel$effect == "time") ncond <- T else ncond <- n
  coefm <- t(sapply(ols,coefficients))
  coef <- lapply(ols,coefficients)
  K <- ncol(coefm)-interc
  res <- lapply(ols,residuals)
  coefb <- apply(coefm,2,sum)/ncond
  xpxm1 <- lapply(ml,function(x){
    if (interc) x <- cbind(1,x[,-1]) else x <- x[-1];solve(crossprod(as.matrix(x)))}
                  )
  D1 <- suml(lapply(oppl(coef,coefb,"-"),function(x) crossprod(t(x))))/(ncond-1)
  sigi <- lapply(res,function(x) sum(x^2)/(length(x)-K-interc))
  sigim <- lapply(sigi,function(x) x*matrix(1,K+interc,K+interc))
  s2xpxm1i <- mapply("*",xpxm1,sigim,SIMPLIFY=F)
  D2 <- matrix(apply(mapply("*",xpxm1,sigi),1,sum)/ncond,ncol=K+interc)
  eig <- prod(eigen(D1-D2)$values>=0)
  if (eig) Delta <- D1-D2 else Delta <- D1  # si le Delta precedent nes pas semi defini positif
  Delta <- list(Delta)
  Wn <- mapply("+",s2xpxm1i,rep(Delta,ncond),SIMPLIFY=F)
  Wnm <- sapply(Wn,solve)
  Wn <- lapply(Wn,solve)
  vcovb <- solve(matrix(apply(Wnm,1,sum),K+interc,K+interc))
  W <- lapply(Wn,function(x) vcovb%*%x)
  beta <- apply(mapply("%*%",W,coef),1,sum)
  df.residuals <- N-ncol(coefm)
  X <- lapply(ml,function(x){
    if (interc) x <- as.matrix(cbind(1,x[,-1])) else x <- as.matrix(x[,-1]);x}
              )
  y <- lapply(ml,function(x) x[,1])
  tss <- tss(unlist(y))
  haty <- oppl(X,beta,"%*%")
  residuals <- mapply("-",y,haty,SIMPLIFY=FALSE)
  ssr <- sum(unlist(residuals)^2)
  names.res <- lapply(residuals,rownames)
  residuals <- lapply(residuals,as.vector)
  zo <- function(x,y){names(x) <- y;x}
  residuals <- unlist(mapply(zo,residuals,names.res,SIMPLIFY=FALSE))
  y <- unlist(y)
  fitted.values <- y-residuals
  names(beta) <- rownames(vcovb) <- colnames(vcovb) <- colnames(coefm)
  swamy <- list(coefficients = beta, residuals = residuals, fitted.values = fitted.values,
                 vcov = vcovb, df.residuals = df.residuals, model = ml, Delta = Delta[[1]])
  swamy <- structure(swamy, pdim = pdim, pmodel = pmodel)
  swamy
}

summary.pvcm <- function(object,...){
  pmodel <- attr(object,"pmodel")
  if (pmodel$model.name=="random"){
    std.err <- sqrt(diag(object$vcov))
    b <- object$coefficients
    z <- b/std.err
    p <- 2*(1-pnorm(abs(z)))
    CoefTable <- cbind(b,std.err,z,p)
    colnames(CoefTable) <- c("Estimate","Std. Error","z-value","Pr(>|z|)")
    object$CoefTable <- CoefTable
  }
  object$ssr <- sum(residuals(object)^2)
  y <-  unlist(lapply(object$model,function(x) x[,1]))
  object$tss <- tss(y)
  object$rsqr <- 1-object$ssr/object$tss
  class(object) <- c("summary.pvcm")
  return(object)
}

print.summary.pvcm <- function(x,digits=max(3, getOption("digits") - 2), width = getOption("width"),...){
  pmodel <- attr(x,"pmodel")
  pdim <- attr(x,"pdim")
  effect <- pmodel$effect
  formula <- pmodel$formula
  model.name <- pmodel$model.name
  cat(paste(effect.pvcm.list[effect]," ",sep=""))
  cat(paste(model.pvcm.list[model.name],"\n",sep=""))
  cat("\nCall:\n")
  print(x$call)
  cat("\n")
  print(pdim)
  cat("\nResiduals:\n")
  print(summary(unlist(residuals(x))))
  if (model.name=="random"){
  cat("\nEstimated mean of the coefficients:\n")
    printCoefmat(x$CoefTable,digits=digits)
    cat("\nEstimated variance of the coefficients:\n")
    print(x$Delta,digits=digits)
  }
  if (model.name=="within"){
    cat("\nCoefficients:\n")
    print(summary(x$coefficients))
  }
  cat("\n")
  cat(paste("Total Sum of Squares: ",signif(x$tss,digits),"\n",sep=""))
  cat(paste("Residual Sum of Squares: ",signif(x$ssr,digits),"\n",sep=""))
  cat(paste("Multiple R-Squared: ",signif(x$rsqr,digits),"\n",sep=""))
  invisible(x)
}
