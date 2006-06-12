nopool <- function(y,...){
  UseMethod("nopool")
}
nopool.formula <- function(y,data=data,effect="individual",...){
  formula <- y
  indexes <- attr(data,"indexes")
  id.index.name <- indexes$id
  time.index.name <- indexes$time
  id <- data[[id.index.name]]
  time <- data[[time.index.name]]
  if (effect=="individual"){
    cond <- id
    other <- time
  }
  else{
    cond <- time
    other <- id
  }

  name.data <- paste(deparse(substitute(data)))
  formod <- paste(deparse(formula),collapse="")
  fortot <- as.formula(paste(formod,"+",id.index.name,"+",time.index.name))
  mf <- model.frame(fortot,data=data)
  
  mtx <- terms(formula)
  X <- model.matrix(mtx,mf)[,-1,drop=F]
  y <- model.response(mf)
  res <- nopool(y,X,cond)
  res
}


nopool.default <- function(y,X,cond,...){
  K <- ncol(X)
  N <- nrow(X)
  namesX <- c("(intercept)",colnames(X))
  coef <- list()
  residuals <- c()
  coefficients <- c()
  std.error <- c()
  valcond <- names(table(cond))
  n <- length(valcond)
  for (i in 1:n){
    suby <- y[cond==valcond[i]]
    subX <- X[cond==valcond[i],]
    z <- lm(suby~subX)
    residuals <- c(residuals,z$residuals)
    coefficients <- rbind(coefficients,z$coefficients)
    std.error <- rbind(std.error,sqrt(diag(vcov(z))))
    
  }
  coefficients <- as.data.frame(coefficients)
  std.error <- as.data.frame(std.error)
  names(coefficients) <- names(std.error) <- namesX
  ssr <- sum(residuals^2)
  df.residuals <- N-n*(K+1)
  result <- list(residuals=residuals,df.residuals=df.residuals,
         ssr=ssr,coefficients=coefficients,std.error=std.error)
  class(result) <- "nopool"
  result
}


plot.nopool <- function(x,...){
  coef <- x$coefficients
  K <- length(coef)
  nc <- K%/%2+1
  nr <- 2
  oldpar <- par(mfrow=c(nr,nc))
  for(i in 1:K){
    hist(coef[[i]],probability=T,main=paste(colnames(coef))[i],xlab=NULL,ylab=NULL)
  }
}

print.nopool <- function(x,digits=3,...){
  print(x$coefficients)
}


summary.nopool <- function(object,...){
  zz <- summary(object$coefficients)
  class(zz) <- c("summary.nopool","table")
  return(zz)
}

print.summary.nopool <- function(x,digits=3,...){
  print.table(x)
}

pooltest <- function(plms,nopool=NULL,effect=FALSE){
  if(is.null(nopool)){
    if (is.null(plms$nopool)) stop("error : no nopool object\n")
    else nopool <- plms$nopool
  }
  if(effect){
    rss <- plms$within$ssr
    uss <- nopool$ssr
    dlr <- plms$within$df.residual
  }
  else{
    rss <- plms$pooling$ssr
    uss <- nopool$ssr
    dlr <- plms$pooling$df.residual
  }
    
  dlu <- nopool$df.residuals
  df1 <- dlr-dlu
  df2 <- dlu
  stat <- (rss-uss)/uss*df2/df1
  pval <- 1-pf(stat,df1,df2)
  names(stat)="F"
  data.name <- paste(deparse(substitute(plms)))
  res <- list(statistic = stat,
              p.value = pval,
              data.name=data.name,
              null.value = "stability",
              alternative = "no stability",
              method = "F statistic")
  class(res) <- "htest"
  res
}
