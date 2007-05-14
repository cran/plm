summary.plm <- function(object,...){
  balanced <- attr(object,"pdim")$balanced
  model.name <- attr(object,"pmodel")$model
  effect <- attr(object,"pmodel")$effect
  std.err <- sqrt(diag(vcov(object)))
  b <- coefficients(object)
  z <- b/std.err
  p <- 2*(1-pnorm(abs(z)))
  CoefTable <- cbind(b,std.err,z,p)
  colnames(CoefTable) <- c("Estimate","Std. Error","z-value","Pr(>|z|)")
  object$CoefTable <- CoefTable
  if (model.name=="random" || model.name=="ht"){
    sigma2 <- object$sigma2
    if (effect=="twoways"){
      sigma2 <- unlist(sigma2)
      sigma2Table <- cbind(var=sigma2,std.dev=sqrt(sigma2),share=sigma2/sum(sigma2))
      rownames(sigma2Table) <- c("idiosyncratic","individual","time")
    }
    else{
      if (balanced) sigma2 <- unlist(sigma2[-1]) else sigma2 <- unlist(sigma2[1:2])
      sigma2Table <- cbind(var=sigma2,std.dev=sqrt(sigma2),share=sigma2/sum(sigma2))
      rownames(sigma2Table) <- c("idiosyncratic","individual")
    }
    object$sigma2Table <- sigma2Table
  }
  object$ssr <- sum(residuals(object)^2)
  object$tss <- tss(object$model[[1]])
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
  instruments <- pmodel$instruments
  model.name <- pmodel$model
  centre("Model Description",length.line)
  cat(paste(effect.plm.list[[effect]],"\n",sep=""))
  cat(paste(model.plm.list[[model.name]]," Model",sep=""))

  if (model.name=="random"){
    random.method <- attr(x,"pmodel")$random.method
    cat(paste(" (",random.method.list[[random.method]],"'s transformation)\n",sep=""))
  }
  else{
    cat("\n")
  }

  if (!is.null(instruments)){
    inst.method <- attr(x,"pmodel")$inst.method
    if (model.name!="ht"){
      cat(paste("Instrumental variable estimation (",
                inst.method.list[[inst.method]],"'s transformation)\n",sep=""))
    }
  }
  
  print.form(formula(x),"Model Formula            : ",length.line)

  if (!is.null(instruments)){
    if (!is.null(endog)){
      print.form(endog,"Endogenous Variables     : ",length.line)
    }
    print.form(instruments,"Instrumental Variables   : ",length.line)
  }

  if (model.name=="ht"){
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
  if (model.name=="random" || model.name=="ht"){
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
  cat(paste("Residual Sum of Squares    : ",signif(x$ssr,digits),"\n",sep=""))
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
  tab.random <- summary(object$random)$CoefTable[,1:2,drop=F]
  tab.dim <- dim(tab.random)
  tab.within <- summary(object$within)$CoefTable[,1:2,drop=F]
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
  
  cat(paste(effect.plm.list[[effect]],"\n",sep=""))
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
