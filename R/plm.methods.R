summary.plm <- function(object,...){
  balanced <- attr(object,"pdim")$balanced
  model.name <- attr(object,"pmodel")$model
  effect <- attr(object,"pmodel")$effect
  std.err <- sqrt(diag(vcov(object)))
  b <- coefficients(object)
  z <- b/std.err
  p <- 2*(1-pnorm(abs(z)))
  CoefTable <- cbind(b,std.err,z,p)
  colnames(CoefTable) <- c("Estimate","Std. Error","t-value","Pr(>|t|)")
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

print.summary.plm <- function(x,digits= max(3, getOption("digits") - 2),width=getOption("width"),...){
  pmodel <- attr(x,"pmodel")
  pdim <- attr(x,"pdim")
  effect <- pmodel$effect
  endog <- pmodel$endog
  instruments <- pmodel$instruments
  model.name <- pmodel$model
#  cat("Model Description\n")
  cat(paste(effect.plm.list[effect]," ",sep=""))
  cat(paste(model.plm.list[model.name]," Model",sep=""))

  if (model.name=="random"){
    random.method <- attr(x,"pmodel")$random.method
    cat(paste(" (",random.method.list[random.method],"'s transformation)\n",sep=""))
  }
  else{
    cat("\n")
  }

  if (!is.null(instruments)){
    inst.method <- attr(x,"pmodel")$inst.method
    if (model.name!="ht"){
      cat(paste("Instrumental variable estimation (",
                inst.method.list[inst.method],"'s transformation)\n",sep=""))
    }
  }
  cat("\nCall:\n")
  print(x$call)

  if (!is.null(instruments)){
    if (!is.null(endog)){
      cat("Endogenous Variables:\n")
      print.form(endog,width)
    }
    cat("Instrumental Variables:\n")
    print.form(instruments,width)
  }

  if (model.name=="ht"){
    cat("\nTime--Varying Variables: ")
    names.x1 <- paste(x$varlist$x1,collapse=",")
    names.x2 <- paste(x$varlist$x2,collapse=",")
    names.z1 <- paste(x$varlist$z1,collapse=",")
    names.z2 <- paste(x$varlist$z2,collapse=",")
    cat(paste("exo (",names.x1,") ",sep=""))
    cat(paste("endo (",names.x2,")\n",sep=""))
    cat("Time--Invariant Variables: ")
    cat(paste("exo (",names.z1,") ",sep=""))
    cat(paste("endo (",names.z2,")\n",sep=""))

  }
  cat("\n")
  print(pdim)
  if (model.name=="random" || model.name=="ht"){
    cat("\nEffects:\n")
    printCoefmat(x$sigma2Table,digits)
    print.theta(x,digits)
  }
  cat("\nResiduals :\n")
  save.digits <- unlist(options(digits=digits))
  on.exit(options(digits=save.digits))
  print(sumres(x))
  
  cat("\nCoefficients :\n")
  printCoefmat(x$CoefTable,digits=digits)
  cat("\n")
  cat(paste("Total Sum of Squares: ",signif(x$tss,digits),"\n",sep=""))
  cat(paste("Residual Sum of Squares: ",signif(x$ssr,digits),"\n",sep=""))
  cat(paste("Multiple R-Squared: ",signif(x$rsq,digits),"\n",sep=""))
  fstat <- x$fstatistic
  cat(paste("F-statistic: ",signif(fstat$statistic),
            " on ",fstat$parameter["df1"]," and ",fstat$parameter["df2"],
            " DF, p-value: ",signif(fstat$p.value,digits),"\n",sep=""))
  invisible(x)
}

