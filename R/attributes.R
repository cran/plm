pvar <- function(x, ...){
  UseMethod("pvar")
}

pvar.pdata.frame <- function(x, ...){
  attr(x,"pvar")
}

pvar.matrix <- function(x,id,time, ...){
  x <- as.data.frame(x)
  pvar.default(x,id,time)
}

pvar.data.frame <- function(x,id,time, ...){
  id <- x[[id]]
  time <- x[[time]]
  pvar.default(x,id,time)
}

pvar.default <- function(x,id,time, ...){
  name.var <- names(x)
  time.variation=rep(TRUE,length(x))
  id.variation=rep(TRUE,length(x))
  for (i in 1:length(x)){
    tv <- tapply(x[[i]],id,myvar)
    ti <- tapply(x[[i]],time,myvar)
    vari <- myvar(x[[i]])
    if (mean(tv[!is.na(tv)])/vari<1E-06) time.variation[i]=FALSE
    if (mean(ti[!is.na(ti)])/vari<1E-06) id.variation[i]=FALSE
  }
  names(id.variation) <- names(time.variation) <- name.var
  dim.var <- list(id.variation=id.variation,time.variation=time.variation)
  class(dim.var) <- "pvar"
  dim.var
}

print.pvar <- function(x,y=NULL, ...){
  varnames <- names(x$time.variation)
  if(any(!x$time.variation)){
    var <- varnames[x$time.variation==FALSE]
    if (!is.null(y)) var <- var[-which(var==y$id)]
    if (length(var)!=0) cat(paste("no time variation   : ",paste(var,collapse=" "),"\n"))
  }
  if(any(!x$id.variation)){
    var <- varnames[x$id.variation==FALSE]
    if (!is.null(y)) var <- var[-which(var==y$time)]
    if(length(var)!=0) cat(paste("no individual variation : ",paste(var,collapse=" "),"\n"))
  }
}

pdim <- function(x, ...){
  UseMethod("pdim")
}

pdim.pdata.frame <- function(x, ...){
  attr(x,"pdim")
}

pdim.data.frame <- function(x,id,time, ...){
  id <- x[[id]]
  time <- x[[time]]
  pdim(id,time)
}

pdim.default <- function(x,y, ...){
  if (length(x) != length(y)) stop("The length of the two vectors differs\n")
  z <- table(x,y)
  Ti <- apply(z,1,sum)
  nt <- apply(z,2,sum)
  n <- nrow(z)
  T <- ncol(z)
  N <- length(x)
  nT <- list(n=n,T=T,N=N)
  id.names <- rownames(z)
  time.names <- colnames(z)
  panel.names <- list(id.names=id.names,time.names=time.names)
  if (any(as.vector(z)==0)){
    balanced <- FALSE
  }
  else balanced <- TRUE
  if (any(as.vector(z)>1)) stop(cat("duplicate couples (time-id)\n"))
  Tint <- list(Ti=Ti,nt=nt)
  z <- list(nT=nT,Tint=Tint,balanced=balanced,panel.names=panel.names)
  class(z) <- "pdim"
  z
}  

print.pdim <- function(x, ...){
  if (x$balanced){
    cat("Balanced Panel\n")
    cat(paste("Number of Individuals        :  ",x$nT$n,"\n",sep=""))
    cat(paste("Number of Time Obserbations  :  ",x$nT$T,"\n",sep=""))
    cat(paste("Total Number of Observations :  ",x$nT$N,"\n",sep=""))
  }
  else{
    cat("Unbalanced Panel\n")
    cat(paste("Number of Individuals        :  ",x$nT$n,"\n",sep=""))
    cat(paste("Number of Time Obserbations  :  from ",min(x$Tint$Ti)," to ",max(x$Tint$Ti),"\n",sep=""))
    cat(paste("Total Number of Observations :  ",x$nT$N,"\n",sep=""))
  }
}

indexes <- function(x){
  if (class(x)[1]!="pdata.frame"){
    stop("indexes function only for pdata.frame\n")
  }
  attr(x,"indexes")
}

print.indexes <- function(x, ...){
  cat(paste("Individual index : ",x$id,"\nTime index       : ",x$time,"\n",sep=""))
}
