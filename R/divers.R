twosls <- function(y, X, W){
  Xhat <- lm(X ~ W)$fit
  if(!is.matrix(Xhat)){
    Xhat <- matrix(Xhat, ncol = 1)
    colnames(Xhat) <- colnames(X)
  }
  model <- lm(y~Xhat-1)
  yhat <- as.vector(crossprod(t(X),coef(model)))
  model$residuals <- y - yhat
  model
}


sumres <- function(x){
  sr <- summary(residuals(x))
  srm <- sr["Mean"]
  if (abs(srm)<1e-10){
    sr <- sr[c(1:3,5:6)]
  }
  sr
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

rbindl <- function(x){
  n <- length(x)
  d <- dim(x[[1]])
  s <- c()
  for (i in 1:n){
    s <- rbind(s,x[[i]])
  }
}

expand.formula <- function(x){
  if (!any(class(x) == "Formula")) stop("not a Formula object")
  extra <- attr(x, "extra")
  if (length(x) != 2) stop("not a two part formula")
  xs <- structure(x, class = "formula")
  has.response <- attr(terms(xs),"response") == 1
  if (has.response){
    y <- x[[2]]
    rhs <- x[[3]]
  }
  else{
    y <- NULL
    rhs <- x[[2]]
  }
  firstpart <- rhs[[2]]
  secondpart <- rhs[[3]]
  
  if (has.response){
    one <- do.call("~",list(y,firstpart))
    two <- do.call("~",list(y,secondpart))
  }
  else{
    one <- do.call("~",list(firstpart))
    two <- do.call("~",list(secondpart))
  }    
  two <- update(one,two)
  result <- as.Formula(one,two)
  class(result) <- c("pFormula", class(result))
  attr(result,"extra") <- extra
  result
}


lev2var <- function(x, ...){
  # take a data.frame and returns a vector of variable names, the
  # names of the vector being the names of the effect
  
  is.fact <- sapply(x, is.factor)
  not.fact <- names(x)[!is.fact]
  names(not.fact) <- not.fact
  x <- x[is.fact]
  wl <- lapply(x,levels)
  # nl is the number of levels for each factor
  nl <- sapply(wl,length)
  # nf is a vector of length equal to the total number of levels
  # containing the name of the factor
  nf <- rep(names(nl),nl)
  result <- unlist(wl)
  names(result) <- nf
  result <- paste(names(result),result,sep="")
  names(nf) <- result
  c(nf,not.fact)
}
