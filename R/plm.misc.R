extractcond <- function(x,...){
  UseMethod("extractcond")
}

extractcond.matrix <- function(x,cond=NULL,effect="individual",data=NULL){
  if(!is.null(data)){
    cond.index.name <- ifelse(effect=="individual",attr(data,"indexes")$id,attr(data,"indexes")$time)
    cond <- data[[cond.index.name]]
  }
  if (is.null(data) & is.null(cond)) stop("argument data or cond missing\n")
  cond
}

extractcond.default <- function(x,cond=NULL,effect="individual",data=NULL){
  if (is.null(cond)){
    if (is.null(data)){
      if (is.null(attr(x,"data"))){ stop("argument data missing with no default\n")}
      else{
        data.name <- attr(x,"data")
        assign("data",eval(parse(text=data.name)))
      }
    }
    cond.index.name <- ifelse(effect=="individual",attr(data,"indexes")$id,attr(data,"indexes")$time)
    cond <- data[[cond.index.name]]
  }
  cond
}

csvect <- function(x,id){
  if (any(is.na(x))){
    z <- na.omit(data.frame(x,id))
    res <- tapply(z$x,z$id,sum)
  }
  else res <- tapply(x,id,sum)
  res
}

papply <- function(x,cond,func){
  na.x <- is.na(x)
  if(!is.matrix(x) || ncol(x)==1){
    cm <- tapply(x,cond,func)
    Cm <- cm[as.character(cond)]
    Cm[na.x] <- NA
    Cm <- as.matrix(Cm)
    cm <- as.matrix(cm)
  }
  else{
    cm <- apply(x,2,tapply,cond,func)
    Cm <- cm[as.character(cond),]
    Cm[na.x] <- NA
  }
  attr(Cm,"cm") <- cm
  Cm
}

pmean <- function(x,cond=NULL,effect="individual",data=NULL){
  cond <- extractcond(x=x,cond=cond,effect=effect,data=data)
  papply(x,cond,mymean)
}

psum <- function(x,cond=NULL,effect="individual",data=NULL){
  cond <- extractcond(x=x,cond=cond,effect=effect,data=data)
  papply(x,cond,mysum)
}

pdiff <- function(x,lag=0){
  data <- attr(x,"data")
  xlagt <- plag(x,order=lag)
  xlagtm1 <- plag(x,order=lag-1)
  xdiff <- xlagt-xlagtm1
  attr(xdiff,"data") <- data
  xdiff
}

plag <- function(x,order=1){
  nT <- length(x)
  data.name <- attr(x,"data")
  assign("data",eval(parse(text=data.name)))
  z <- attr(data,"indexes")
  time.name <- z$time
  id.name <- z$id
  time <- eval(parse(text=paste(data.name,"$",time.name)))
  id <- eval(parse(text=paste(data.name,"$",id.name)))
  res <- c(rep(NA,order),x[1:(nT-order)])
  lagid <- id-c(rep(NA,order),id[1:(nT-order)])
  res[lagid!=0]=NA
  attr(res,"data") <- data.name
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



varcond <- function(x,id){
 tapply(x,id,myvar)
}

twosls <- function(y,X,W,intercept=FALSE){
  Xhat <- lm(X~W)$fit
  if(intercept){
    model <- lm(y~Xhat)
    residuals <- y-cbind(1,X)%*%model$coef
  }
  else{
    model <- lm(y~Xhat-1)
    residuals <- y-as.matrix(X)%*%model$coef
  }
  model$residuals <- as.vector(residuals)
  model
}
  


  
