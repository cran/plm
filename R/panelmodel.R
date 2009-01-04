
terms.panelmodel <- function(x,...){
  terms(formula(x))
}

vcov.panelmodel <- function(object,...){
  object$vcov
}

fitted.panelmodel <- function(object,...){
  object$fitted.values 
}

residuals.panelmodel <- function(object,...){
  object$residuals
}

df.residual.panelmodel <- function(object,...){
  object$df.residual
}

coef.panelmodel <- function(object,...){
  object$coefficients
}

plmformat <- function(object,names,rnames,
                      df.residual,model.name,pdim,pmodel,indexes,cl,...){
  interc <- pmodel$has.intercept
  coefficients <- coefficients(object)
  residuals <- residuals(object)
  fitted.values <- fitted.values(object)
  attr(residuals,"cm") <- attr(fitted.values,"cm") <- NULL
  if(model.name=="within"){
    vcov <- vcov(object)*object$df.residual/df.residual
  }
  else{
    vcov <- vcov(object)
  }
  cl$model.name <- model.name
  response.name <- deparse(pmodel$formula[[2]])
  mymodel <- object$model
  mymodel <- object$model[[2]]
  mymodel <- data.frame(object$model[[1]],mymodel)
  if (model.name == "random" && interc == TRUE) names(mymodel) <- c(response.name,c("(intercept)",names))
  else names(mymodel) <- c(response.name,names)
  attr(mymodel,"terms") <- terms(pmodel$formula)
  rownames(mymodel) <- names(residuals) <- names(fitted.values) <- rnames
  if (interc && model.name != "within") coefnames <- c("(intercept)",names) else coefnames <- names
  names(coefficients) <- rownames(vcov) <- colnames(vcov) <- coefnames
  object <- list(coefficients=coefficients,residuals=residuals,
                 fitted.values=fitted.values,vcov=vcov,df.residual=df.residual,model=mymodel,
                 call=cl)
  object <- structure(object,pdim=pdim,pmodel=pmodel,indexes=indexes,class=c("plm","panelmodel"))
}

print.panelmodel <- function(x,digits=max(3, getOption("digits") - 2), width = getOption("width"),...){
  cat("\nModel Formula: ")
  print(formula(x))
  cat("\nCoefficients:\n")
  print(coef(x),digits=digits)
  cat("\n")
  invisible(x)
}

print.form <- function(x,length.line){
  x <- deparse(x,width.cutoff=length.line)
  n <- length(x)
  cat(paste(x[1],"\n",sep=""))
  if (n>1){
    for (i in 2:n){
      cat(paste(x[i],"\n",sep=""))
    }
  }
}

print.theta <- function(x,digits){
  effect <- attr(x,"pmodel")$effect
  pdim <- attr(x,"pdim")
  if (effect!="twoways"){
    if (pdim$balanced){
      cat(paste("theta: ",signif(x$theta,digits)," \n"))
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

## model.matrix.panelmodel <- function(object, ...){
##   data <- model.frame(object, ...)
##   form <- attr(data,"terms")
##   m <- model.matrix(form,data)
##   if (attr(object,"pmodel")$model.name == "random" && attr(object,"pmodel")$has.intercept == TRUE){
##     m[,1] <- data[["(intercept)"]]
##   }
##   if (attr(object,"pmodel")$model.name == "within" && attr(object,"pmodel")$has.intercept == TRUE){
##     m <- m[, -1, drop = FALSE]
##   }
##   m
## }

model.matrix.panelmodel <- function(object, ...){
  m <- as.matrix(model.frame(object, ...)[, -1, drop = FALSE])
#  if (attr(object,"pmodel")$model.name == "random" && attr(object,"pmodel")$has.intercept == TRUE){
#    m[,1] <- data[["(intercept)"]]
#  }
  if (attr(object,"pmodel")$model.name %in% c("between","pooling") && attr(object,"pmodel")$has.intercept == TRUE){
    m <- cbind('(intercept)'=1,m)
  }
  m
}


