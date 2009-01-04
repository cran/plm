papply <- function(x, ...){
  UseMethod("papply")
}

papply.pserie <- function(x,func,effect="individual",...){
  na.x <- is.na(x)
  data.name <- attr(x,"data")
  sc <- sapply(sys.calls(),function(x) as.character(x[[1]]))
  leframe <- which(!is.na(match(sc,c("pgmm","pgls","pgmm","pvcm"))))
  leframe <- ifelse(length(leframe)==1,leframe,1)
  
  data <- get(data.name,sys.frame(which=leframe))
  indexes <- attr(data,"indexes")
  cond <- switch(effect,
                 "individual"= data[[indexes$id]],
                 "time"= data[[indexes$time]],
                 stop("effect must be individual or time")
                 )
  cm <- tapply(x,cond,func)
  Cm <- cm[as.character(cond)]
  Cm[na.x] <- NA
  attr(Cm,"cm") <- cm
  class(Cm) <- "pserie"
  attr(Cm,"data") <- data.name
  Cm
}

papply.matrix <- function(x,func,cond,...){
  na.x <- is.na(x)
  cm <- apply(x,2,tapply,cond,func)
  Cm <- cm[as.character(cond),,drop=F]
  Cm[na.x] <- NA
  attr(Cm,"cm") <- cm
  Cm
}

papply.default <- function(x,func,cond,...){
  na.x <- is.na(x)
  cm <- tapply(x,cond,func)
  Cm <- cm[as.character(cond)]
  Cm[na.x] <- NA
  attr(Cm,"cm") <- cm
  Cm
}

Between <- function(x,...){
  UseMethod("Between")
}

Between.pserie <- function(x,effect="individual", ...){
  if (is.numeric(x)){
    res <- papply(x,mymean,effect)
  }
  else{
    stop("The Between function only applies to numeric vectors\n")
  }
  attr(res,"cm") <- NULL
  res
}

Between.default <- function(x,cond, ...){
  if (is.numeric(x)){
    res <- papply(x,mymean,cond)
  }
  else{
    stop("The Between function only applies to numeric vectors\n")
  }
  attr(res,"cm") <- NULL
  res
}

Between.matrix <- function(x,cond, ...){
  if (is.numeric(x)){
    res <- papply(x,mymean,cond)
  }
  else{
    stop("The Between function only applies to numeric vectors\n")
  }
  attr(res,"cm") <- NULL
  res
}

between <- function(x,...){
  UseMethod("between")
}

between.pserie <- function(x,effect="individual", ...){
  data.name <- attr(x,"data")
  classx <- class(x)
  sc <- sapply(sys.calls(),function(x) as.character(x[[1]]))
  leframe <- which(!is.na(match(sc,c("pgmm","pgls","pgmm","pvcm"))))
  leframe <- ifelse(length(leframe)==1,leframe,1)

  data <- get(data.name,sys.frame(which=leframe))
  indexes <- attr(data,"indexes")
  cond <- switch(effect,
                 "individual"= data[[indexes$id]],
                 "time"= data[[indexes$time]],
                 stop("effect must be individual or time")
               )
  if (is.numeric(x)){
    res <- tapply(x,cond,mymean)
  }
  else{
    stop("The between function only applies to numeric vectors\n")
  }
  res
}

between.matrix <- function(x,cond, ...){
  if (is.numeric(x)){
    res <- apply(x,2,tapply,cond,mymean)
  }
  else{
    stop("The between function only applies to numeric vectors\n")
  }
  res
}

between.default <- function(x,cond, ...){
  if (is.numeric(x)){
    res <- tapply(x,cond,mymean)
  }
  else{
    stop("The between function only applies to numeric vectors\n")
  }
  res
}

within <- function(x,...){
  UseMethod("within")
}

within.pserie <- function(x,effect="individual", ...){
  if (!is.numeric(x)) stop("the within function only applies to numric vectors")
  res <- switch(effect,
                "individual"= x-Between(x,effect),
                "time"=  x-Between(x,effect),
                "twoways"= x-Between(x,"individual")-Between(x,"time")+mean(x),
                stop("effect must be individual, time or twoways")
                )
  res
}

within.default <- function(x,cond, ...){
  res <- x-papply(x,mymean,cond)
  attr(res,"cm") <- NULL
  res
}

within.matrix <- function(x,cond, ...){
  if (is.numeric(x)){
    res <- x-papply(x,mymean,cond)
  }
  else{
    stop("The within function only applies to numeric vectors\n")
  }
  attr(res,"cm") <- NULL
  res
}

diff.pserie <- function(x,lag=0,...){
  if (!is.numeric(x)) stop("diff meaningfull only for numeric vectors\n")
  xlagt <- lag(x,k=lag)
  xlagtm1 <- lag(x,k=lag+1)
  xdiff <- xlagt-xlagtm1
  xdiff
}

lag.pserie <- function(x,k=1,...){
  
  N <- length(x)
  data.name <- attr(x,"data")
  classx <- class(x)
  if (is.factor(x)){
    levelsx <- levels(x)
    nlevelsx <- length(levelsx)
  }
  sc <- sapply(sys.calls(),function(x) as.character(x[[1]]))
  leframe <- which(!is.na(match(sc,c("pgmm","pggls","plm","pvcm"))))
  if (length(leframe>1)) leframe <- leframe[length(leframe)]
  leframe <- ifelse(length(leframe)==1,leframe,1)

  # il faut evaluer  dans le frame en position 1,
  # celui qui correspond à la fonction pgmm
  # print(ls(sys.frame(which=1)))
  data <- get(data.name,sys.frame(which=leframe))
  id.name <- attr(data,"indexes")$id
  time.name <- attr(data,"indexes")$time
  id <- data[[id.name]]
  time <- data[[time.name]]
  id <- as.numeric(id)
  res <- c(rep(NA,k),x[1:(N-k)])
  lagid <- id-c(rep(NA,k),id[1:(N-k)])
  res[lagid!=0]=NA
  if (is.factor(x)){
    res <- factor(res,levels=1:nlevelsx,labels=levelsx)
  }
  attr(res,"data") <- data.name
  class(res) <- classx
  res
}


myvar <- function(x){
  if(any(is.na(x))) x <- x[!is.na(x)]
  n <- length(x)
  z <- switch(as.character(n),
              "0"=NA,
              "1"=0,
              var(x))
  z
}

mymean <- function(x){
  if(any(is.na(x))) x <- x[!is.na(x)]
  n <- length(x)
  z <- switch(as.character(n),
              "0"=NA,
              mean(x))
  z
}

mysum <- function(x){
  if(any(is.na(x))) x <- x[!is.na(x)]
  n <- length(x)
  z <- switch(as.character(n),
              "0"=NA,
              sum(x))
  z
}

twosls <- function(y,X,W,intercept=FALSE){
  Xhat <- lm(X~W)$fit
  if(!is.matrix(Xhat)){
    Xhat <- matrix(Xhat,ncol=1)
    colnames(Xhat) <- colnames(X)
  }
  if(intercept){
    model <- lm(y~Xhat)
    residuals <- y-as.vector(cbind(1,X)%*%model$coef)
  }
  else{
    model <- lm(y~Xhat-1)
    residuals <- y-as.vector(as.matrix(X)%*%model$coef)
  }
  model$residuals <- as.vector(residuals)
  model
}
  
tss <- function(x){
  n <- length(x)
  sum(x^2)-n*mean(x)^2
}


sumres <- function(x){
  sr <- summary(residuals(x))
  srm <- sr["Mean"]
  if (abs(srm)<1e-10){
    sr <- sr[c(1:3,5:6)]
  }
  sr
}

fixef.plm <- function(object, effect = c("individual","time"), ...){
  effect <- match.arg(effect)
  pmodel <- attr(object,"pmodel")
  model.name <- pmodel$model
  if (model.name!="within"){
    stop("fixef function only implemented for within models\n")
  }
  else{
    if (pmodel$effect!="twoways"){
      fixef <- attr(object$fixef,"cm")
    }
    else{
      if (is.null(effect)){
        stop("effect should be indicated\n")
      }
      else{
        if (effect=="individual"){
          fixef <- attr(object$fixef$id,"cm")
        }
        else{
          fixef <- attr(object$fixef$time,"cm")
        }
      }
    }
  }
  if (pmodel$effect!="twoways"){
    bet <- update(object,model="between")
    xb <- model.matrix(bet)
    sigma2 <- sum(residuals(bet)^2)/df.residual(bet)
    # the between model contains a constant, the within one don't
    vcov <- vcov(object)[names(coef(bet))[-1],names(coef(bet))[-1]]
    T <- attr(bet,"pdim")$nT$T
    s2 <- sum(resid(object)^2)/df.residual(object)
#    sefixef <- sqrt(apply(xb,1,function(x) t(x)%*%vcov%*%x))
    sefixef <- sqrt(s2/T+apply(xb[,rownames(vcov),drop=FALSE],1,function(x) t(x)%*%vcov%*%x))
    intercept <- object$alpha
    fixef <- structure(fixef-intercept,se=sefixef,intercept=intercept,class="fixef")
  }
  else{
    class(fixef) <- "fixef"
  }
  fixef
}

print.fixef <- function(x,digits= max(3, getOption("digits") - 2),width=getOption("width"),...){
  if (!is.null(attr(x,"intercept"))){
#    if (diff) intercept <- attr(x,"intercept")
    attr(x,"se") <- attr(x,"intercept") <- attr(x,"class") <- NULL
  }
  else{
    attr(x,"class") <- NULL
  }
  print.default(x)
}

summary.fixef <- function(object,...){
  if (is.null(attr(object,"intercept"))) stop("summary method not defined for two-ways effects")
  se <- attr(object,"se")
  alpha <- attr(object,"intercept")
  zvalue <- (object)/se
  res <- cbind(object,se,zvalue,(1-pnorm(abs(zvalue)))*2)
  colnames(res) <- c("Estimate","Std. Error","t-value","Pr(>|t|)")
  class(res) <- "summary.fixef"
  res
}

print.summary.fixef <- function(x,digits= max(3, getOption("digits") - 2),width=getOption("width"),...){
  printCoefmat(x,digits=digits)
}
  

suml <- function(x){
  n <- length(x)
  if (!is.null(dim(x[[1]]))){
    d <- dim(x[[1]])
    s <- matrix(0,d[1],d[2])
    for (i in 1:n){
      s <- s+x[[i]]
    }
  }
  else{
    s <- rep(0,length(x[[n]]))
    for (i in 1:n){
      s <- s+x[[i]]
    }
  }
  s
}

oppl <- function(x,y,func){
  n <- length(x)
  z <- list()
  if (!is.list(y)){
    for (i in 1:n){
      t <- paste("\"",func,"\"","(x[[i]],y)",sep="")
      z[[i]] <- eval(parse(text=t))
    }
  }
  else{
    for (i in 1:n){
      t <- paste("\"",func,"\"","(x[[i]],y[[i]])",sep="")
      z[[i]] <- eval(parse(text=t))
    }
  }
  z
}

is.one.side.formula <- function(x){
  class(x)=="formula" && length(x)==2
}

rbindl <- function(x){
  n <- length(x)
  d <- dim(x[[1]])
  s <- c()
  for (i in 1:n){
    s <- rbind(s,x[[i]])
  }
}
