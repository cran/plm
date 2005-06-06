plmformat <- function(object,names,formula,df.residual,...){
  sobject <- summary(object)
  coefficients <- object$coefficients
  residuals <- object$residuals
  attr(residuals,"cm") <- NULL
  cov.unscaled <- as.matrix(sobject$cov.unscaled)
  ssr <- sum(residuals^2)
  res <- list(coefficients=coefficients,residuals=residuals,cov.unscaled=cov.unscaled,ssr=ssr,formula=formula,df.residual=df.residual)
  names(res$coefficients) <- rownames(res$cov.unscaled) <- colnames(res$cov.unscaled) <- names
  class(res) <- c("plm")
  res
}

print.plm <- function(x,digits=3,...){
  cat("\nModel Formula: ")
  print(x$formula)
  cat("\nCoefficients:\n")
  print(x$coefficients,digits=digits)
  cat("\n")
  invisible(x)
}

summary.plm <- function(object,...){
  vcov <- object$cov.unscaled*object$ssr/object$df.residual
  std.err <- sqrt(diag(vcov))
  b <- object$coefficients
  z <- b/std.err
  p <- 2*(1-pnorm(abs(z)))
  CoefTable <- cbind(b,std.err,z,p)
  colnames(CoefTable) <- c("Estimate","Std. Error","z-value","Pr(>|z|)")
  object$CoefTable <- CoefTable
  class(object) <- "summary.plm"
  return(object)
}

print.summary.plm <- function(x,digits=3,...){
  save.digits <- unlist(options(digits=digits))
  on.exit(options(digits=save.digits))
  cat("\nModel formula: ")
  print(x$formula)
  cat("\nResiduals:\n")
  print(summary(residuals(x)))
  cat("\n")
  printCoefmat(x$CoefTable,digits=digits)
  invisible(x)
}

print.plms <- function(x,digits=3,...){
  for (i in 1:length(x)){
    cat(paste("Model ",names(x)[i]," :\n"))
    print(x[[i]],digits=digits)
  }
  invisible(x)
}

summary.plms <- function(object,...){
  tab.random=summary(object$random)$CoefTable[,1:2]
  tab.dim=dim(tab.random)
  tab.within=summary(object$within)$CoefTable
  names.within=rownames(tab.within)
  tab.within=tab.within[,1:2]
  tab.between=summary(object$between)$CoefTable[,1:2]
  names.random=rownames(tab.random)
  names.between=rownames(tab.between)
  miss.within=names.random[!names.random%in%names.within]
  miss.between=names.random[!names.random%in%names.between]
  CoefTable=cbind(matrix(NA,tab.dim[1],tab.dim[2]),matrix(NA,tab.dim[1],tab.dim[2]),tab.random)
  CoefTable[names.between,1:2]=tab.between
  CoefTable[names.within,3:4]=tab.within
  colnames(CoefTable)=c("between","bse","within","wse","random","rse")
  object$CoefTable=CoefTable
  class(object)="summary.plms"
  object
}

print.summary.plms <- function(x,digits=3,...){
  printCoefmat(x$CoefTable,digits=digits,na.print=".")
  invisible(x)
}

