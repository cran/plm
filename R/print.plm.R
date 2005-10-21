length.line <- 80

vcov.plm <- function(object,...){
  object$cov.unscaled*object$ssr/object$df.residual
}

formula.plm <- function(x,...){
  x$formula
}


## coefficients.plm <- function(object,...){
##   object$coefficient
## }

## residuals.plm <- function(object,...){
##   object$residuals
## }

plmformat2 <- function(object,names,formula,
                      df.residual,model.name,tss,pdim,pmodel,indexes,cl,...){
  sobject <- summary(object)
  coefficients <- object$coefficients
  residuals <- object$residuals
  attr(residuals,"cm") <- NULL
  cov.unscaled <- as.matrix(sobject$cov.unscaled)
  ssr <- sum(residuals^2)
  cl$model.name=model.name
  res <- list(coefficients=coefficients,residuals=residuals,cov.unscaled=cov.unscaled,
              ssr=ssr,formula=formula,df.residual=df.residual,model.name=model.name,tss=tss,
              terms=object$terms,call=cl)
  names(res$coefficients) <- rownames(res$cov.unscaled) <- colnames(res$cov.unscaled) <- names
  attr(res,"pdim") <- pdim
  attr(res,"pmodel") <- pmodel
  attr(res,"indexes") <- indexes
  class(res) <- c("plm")
  res
}


plmformat <- function(object,names,formula,
                      df.residual,model.name,tss,pdim,pmodel,indexes,cl,...){
  sobject <- summary(object)
  residuals <- as.vector(object$residuals)
  attr(residuals,"cm") <- NULL
  cov.unscaled <- as.matrix(sobject$cov.unscaled)
  ssr <- sum(residuals^2)
  cl$model.name=model.name

  object$cov.unscaled <- cov.unscaled
  object$ssr <- ssr
  object$df.residual <- df.residual
  object$model.name <- model.name
  object$tss <- tss
  object$call <- cl
  object$formula <- formula
  object$residuals <- residuals
  ##   res <- list(coefficients=coefficients,residuals=residuals,cov.unscaled=cov.unscaled,
##               ssr=ssr,formula=formula,df.residual=df.residual,model.name=model.name,tss=tss,
##               terms=object$terms,call=cl)

  names(object$coefficients) <- rownames(object$cov.unscaled) <- colnames(object$cov.unscaled) <- names
  attr(object,"pdim") <- pdim
  attr(object,"pmodel") <- pmodel
  attr(object,"indexes") <- indexes
  class(object) <- c("plm","lm")
  object
}

print.plm <- function(x,digits=5,...){
  cat("\nModel Formula: ")
  print(x$formula)
  cat("\nCoefficients:\n")
  print(x$coefficients,digits=digits)
  cat("\n")
  invisible(x)
}

summary.plm <- function(object,...){
  balanced <- attr(object,"pdim")$balanced
  vcov <- object$cov.unscaled*object$ssr/object$df.residual
  std.err <- sqrt(diag(vcov))
  b <- object$coefficients
  z <- b/std.err
  p <- 2*(1-pnorm(abs(z)))
  CoefTable <- cbind(b,std.err,z,p)
  colnames(CoefTable) <- c("Estimate","Std. Error","z-value","Pr(>|z|)")
  object$CoefTable <- CoefTable
  effect <- attr(object,"pmodel")$effect
  if (object$model.name=="random" || object$model.name=="ht"){
    sigma2 <- object$sigma2
    if (effect=="twoways"){
      sigma2 <- unlist(sigma2)
      sigma2Table <- cbind(var=sigma2,std.dev=sqrt(sigma2),share=sigma2/sum(sigma2))
      rownames(sigma2Table) <- c("idiosyncratic","individual","time")
    }
    else{
      if (balanced) sigma2 <- unlist(sigma2[-1]) else sigma2 <- unlist(sigma2)
      sigma2Table <- cbind(var=sigma2,std.dev=sqrt(sigma2),share=sigma2/sum(sigma2))
      rownames(sigma2Table) <- c("idiosyncratic","individual")
    }
  object$sigma2Table <- sigma2Table
  }
  object$rsqr <- 1-object$ssr/object$tss
  object$fstatistic <- Ftest(object)
  class(object) <- c("summary.plm","plm")
  return(object)
}

print.summary.plm <- function(x,digits=5,length.line=70,...){
  pmodel <- attr(x,"pmodel")
  pdim <- attr(x,"pdim")
  effect <- pmodel$effect
  endog <- pmodel$endog
  centre("Model Description",length.line)
  cat(paste(effect.list[[effect]],"\n",sep=""))
  cat(paste(model.list[[x$model.name]]," Model",sep=""))

  if (x$model.name=="random"){
    random.method <- attr(x,"pmodel")$random.method
    cat(paste(" (",random.method.list[[random.method]],"'s transformation)\n",sep=""))
  }
  else{
    cat("\n")
  }

  if (!is.null(endog$instruments)){
    inst.method <- attr(x,"pmodel")$inst.method
    if (x$model.name!="ht"){
      cat(paste("Instrumental variable estimation (",
                inst.method.list[[inst.method]],"'s transformation)\n",sep=""))
    }
  }
  
  print.form(x$formula,"Model Formula             : ",length.line)

  if (!is.null(endog$instruments)){
    if (!is.null(endog$endog)){
      print.form(endog$endog,"Endogenous Variables    : ",length.line)
    }
    print.form(endog$instruments,"Instrumental Variables    : ",length.line)
  }

  if (x$model.name=="ht"){
    cat("Time--Varying Variables    \n")
    names.x1 <- paste(x$varlist$x1,collapse=",")
    names.x2 <- paste(x$varlist$x2,collapse=",")
    names.z1 <- paste(x$varlist$z1,collapse=",")
    names.z2 <- paste(x$varlist$z2,collapse=",")
    cat(paste("    exogenous variables   : ",names.x1,"\n"))
    cat(paste("    endogenous variables  : ",names.x2,"\n"))
    cat("Time--Invariant Variables  \n")
    cat(paste("    exogenous variables   : ",names.z1,"\n"))
    cat(paste("    endogenous variables  : ",names.z2,"\n"))

  }
  
  centre("Panel Dimensions",length.line)
  print(pdim)
  if (x$model.name=="random" || x$model.name=="ht"){
    centre("Effects",length.line)
    printCoefmat(x$sigma2Table)
    print.theta(x,digits)
  }
  centre("Residuals",length.line)
  save.digits <- unlist(options(digits=digits))
  on.exit(options(digits=save.digits))
  print(summary(residuals(x)))
  centre("Coefficients",length.line)
  printCoefmat(x$CoefTable,digits=digits)

  centre("Overall Statistics",length.line)
  cat(paste("Total Sum of Squares       : ",signif(x$tss,digits),"\n",sep=""))
  cat(paste("Sum of Squares Residuals   : ",signif(x$ssr,digits),"\n",sep=""))
  cat(paste("Rsq                        : ",signif(x$rsq,digits),"\n",sep=""))
  cat(paste("F                          : ",signif(x$fstatistic$statistic),"\n",sep=""))
  cat(paste("P(F>0)                     : ",signif(x$fstatistic$p.value),"\n",sep=""))

  cat(paste(trait(length.line),"\n"))
  invisible(x)
}

print.plms <- function(x,digits=5,...){
  for (i in 1:length(x)){
    cat(paste("Model ",names(x)[i]," :\n"))
    print(x[[i]],digits=digits)
  }
  invisible(x)
}

summary.plms <- function(object,...){
  tab.random <- summary(object$random)$CoefTable[,1:2]
  tab.dim <- dim(tab.random)
  tab.within <- summary(object$within)$CoefTable[,1:2]
  names.within <- rownames(tab.within)
  names.random <- rownames(tab.random)
  miss.within <- names.random[!names.random%in%names.within]
  CoefTable <- cbind(matrix(NA,tab.dim[1],tab.dim[2]),tab.random)
  CoefTable[names.within,1:2] <- tab.within
  colnames(CoefTable) <- c("within","wse","random","rse")
  object$CoefTable <- CoefTable
  object$phtest <- phtest(object)
  object$pFtest <- pFtest(object)
  object$plmtest <- plmtest(object)
  class(object) <- "summary.plms"
  object
}

print.summary.plms <- function(x,digits=5,length.line=70,...){
  pmodel <- attr(x,"pmodel")
  pdim <- attr(x,"pdim")
  effect <- pmodel$effect
  endog <- pmodel$endog

  centre("Model Description",length.line)
  
  cat(paste(effect.list[[effect]],"\n",sep=""))
  cat("\n")

  print.form(pmodel$formula,"Model Formula        : ",length.line)

  if (!is.null(endog$instruments)){
    if (!is.null(endog$endog)){
      print.form(endog$endog,"Endogenous Variables : ",length.line)
    }
    print.form(endog$instruments,"Instrument Variables : ",length.line)
  }
  
  centre("Panel Dimensions",length.line)
  print(pdim)

  centre("Coefficients",length.line)
  printCoefmat(x$CoefTable,digits=digits,na.print=".")
  invisible(x)
  centre("Tests",length.line)
  phtest <- x$phtest
  cat("Hausman Test                   : ",names(phtest$statistic),
      "(",phtest$parameter,") = ",x$phtest$statistic,
      " (p.value=",phtest$p.value,")\n",sep="")
  pFtest <- x$pFtest
  cat("F Test                         : ",names(pFtest$statistic),
      "(",pFtest$parameter[[1]],",",pFtest$parameter[[2]],") = ",
      pFtest$statistic," (p.value=",pFtest$p.value,")\n",sep="")
  plmtest <- x$plmtest
  cat("Lagrange Multiplier Test       : ",names(plmtest$statistic),
      "(",plmtest$parameter,") = ",plmtest$statistic," (p.value=",plmtest$p.value,")\n",sep="")

  cat(paste(trait(length.line),"\n"))
  plmtest <- x$plmtest

}

print.theta <- function(x,digits){
  effect <- attr(x,"pmodel")$effect
  pdim <- attr(x,"pdim")
  if (effect!="twoways"){
    if (pdim$balanced){
      cat(paste("theta   : ",signif(x$theta,digits)," \n"))
    }
    else{
      cat("theta  : \n")
      print(summary(x$theta))
    }
  }
  else{
    if(pdim$balanced){
      cat(paste("theta  : ",signif(x$theta$id,digits)," (id) ",signif(x$theta$time,digits)," (time) ",signif(x$theta$total,digits)," (total)\n",sep=""))
    }
  }
}

trait <- function(n=length.line,char="_"){
  marge=""
  for (i in 1:n){
    marge <- paste(marge,char,sep="")
  }
}

centre <- function(x,total=length.line,char="_"){
  cat(paste(trait(total,char),"\n"))
  n <- nchar(x)
  marge.gauche <- trunc((total-n-2)/2)
  marge.droite <- total-n-2-marge.gauche
  marge.gauche <- trait(marge.gauche,char=char)
  marge.droite <- trait(marge.droite,char=char)
  cat(paste(marge.gauche," ",x," ",marge.droite,"\n",sep=""))

}

print.form <- function(x,title,length.line){
  length.line <- length.line-nchar(title)-10
  x <- deparse(x,width.cutoff=length.line)
  n <- length(x)
  nt <- nchar(title)
  cat(paste(title,x[1],"\n",sep=""))
  if (n>1){
    for (i in 2:n){
      cat(paste(trait(nt," "),x[i],"\n",sep=""))
    }
  }
}

