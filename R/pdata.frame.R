pdata.frame <- function(x,id,time=NULL){
  data.name <- paste(deparse(substitute(x)))
  time.index.name <- paste(deparse(substitute(time)))
  id.index.name <- paste(deparse(substitute(id)))

  options(warn=-1)
  id.integer <- !is.na(as.numeric(id.index.name))
  options(warn=0)
  if(id.integer && length(id)==1){
    if(time.index.name!="NULL"){warning("The time argument will be ignored\n")}
    N <- nrow(x)
    if( (N%%id)!=0){
      stop("unbalanced panel, the id variable should be indicated\n")
    }
    else{
      T <- N%/%id
      n <- N%/%T
      time <- rep(1:T,n)
      id <- rep(seq(1:n),rep(T,n))
      id.index.name <- "id"
      time.index.name <- "time"
      x[[id.index.name]] <- id <- as.factor(id)
      x[[time.index.name]] <- time <- as.factor(time)
    }
  }
  else{

    id.index.name <- paste(deparse(substitute(id)))
    id <- x[[id.index.name]] <- as.factor(x[[id.index.name]])
    if (time.index.name=="NULL"){
      Ti <- table(id)
      n <- length(Ti)
      time <- c()
      for (i in 1:n){
        time <- c(time,1:Ti[i])
      }
      time.index.name <- "time"
      x[[time.index.name]] <- time <- as.factor(time)
    }
    else{
      x[[time.index.name]] <- time <- as.factor(x[[time.index.name]])
    }
  }
  indexes <- list(id=id.index.name,time=time.index.name)
  var.names <- names(x)
  pdim <- pdimcheck(id,time)
  
  for (i in 1:length(x)){
    attr(x[[i]],"data") <- data.name
  }
  pvar <- pvarcheck(id,time,x)
  structure(x,class=c("pdata.frame","data.frame"),pvar=pvar,pdim=pdim,indexes=indexes)
}

pvarcheck <- function(id,time,data){
  name.var <- names(data)
  id.index.name <- paste(substitute(id))
  time.index.name <- paste(substitute(time))
  exist.time <- exists(time.index.name,sys.frame(sys.parent()))
  exist.id <- exists(id.index.name,sys.frame(sys.parent()))
  if (!exist.id) id <- data[[id.index.name]]
  if (!exist.time) time <- data[[time.index.name]]
  data <- as.data.frame(data)
  time.variation=rep(TRUE,length(data))
  id.variation=rep(TRUE,length(data))
  for (i in 1:length(data)){
    tv <- tapply(data[[i]],id,myvar)
    ti <- tapply(data[[i]],time,myvar)
    vari <- myvar(data[[i]])
    if (mean(tv[!is.na(tv)])/vari<1E-06) time.variation[i]=FALSE
    if (mean(ti[!is.na(ti)])/vari<1E-06) id.variation[i]=FALSE
  }
  names(id.variation) <- names(time.variation) <- name.var
  dim.var <- list(id.variation=id.variation,time.variation=time.variation)
  dim.var
}


pdimcheck <- function(id,time,data=NULL){
  id.index.name <- paste(substitute(id))
  time.index.name <- paste(substitute(time))
  if(!exists(id.index.name)){
    if(is.null(data)){
      stop(paste("variable ",id.index.name," not found\n",sep=""))
    }
    else{
      id <- data[[id.index.name]]
    }
  }
  if(!exists(time.index.name)){
    if(is.null(data)){
      stop(paste("variable ",time.index.name," not found\n",sep=""))
    }
    else{
      time <- data[[time.index.name]]
    }
  }
  z <- table(id,time)
  Ti <- apply(z,1,sum)
  nt <- apply(z,2,sum)
  n <- nrow(z)
  T <- ncol(z)
  N <- length(id)
  nT <- list(n=n,T=T,N=N)
  id.names <- rownames(z)
  time.names <- colnames(z)
  panel.names <- list(id.names=id.names,time.names=time.names)
  if (any(as.vector(z)==0)){
    balanced <- FALSE
    cat("unbalanced panel\n")
  }
  else balanced <- TRUE
  if (any(as.vector(z)>1)) stop(cat("duplicate couples (time-id)\n"))
  Tint <- list(Ti=Ti,nt=nt)
  list(nT=nT,Tint=Tint,balanced=balanced,panel.names=panel.names)
}

summary.pdata.frame <- function(object,...){
  zz <- summary.data.frame(object)
  attr(zz,"pdim") <- attr(object,"pdim")
  attr(zz,"pvar") <- attr(object,"pvar")
  attr(zz,"indexes") <- attr(object,"indexes")
  attr(zz,"varnames") <- names(object)
  
  class(zz) <- c("summary.pdata.frame","table")
  return(zz)
}

print.summary.pdata.frame <- function(x,...){
  pdim <- attr(x,"pdim")
  indexes <- attr(x,"indexes")
  pvar <- attr(x,"pvar")
  varnames <- attr(x,"varnames")

  if(pdim$balanced) bal <- "Balanced" else bal <- "Unbalanced"
  cat(paste("\n",bal," panel with ",pdim$nT$n," individuals and ",pdim$nT$T," time observations\n\n",sep=""))
  cat(paste("Individual index : ",indexes$id,"\nTime index      : ",indexes$time,"\n\n",sep=""))
  
  if(any(!pvar$time.variation)){
    var <- varnames[pvar$time.variation==FALSE]
    var <- var[-which(var==indexes$id)]
    if (length(var)!=0) cat(paste("The following variables don't have time variation   : ",paste(var,collapse=" "),"\n\n"))
  }
  if(any(!pvar$id.variation)){
    var <- varnames[pvar$id.variation==FALSE]
    var <- var[-which(var==indexes$time)]
    if(length(var)!=0) cat(paste("The following variables don't have individual variation : ",paste(var,collapse=" "),"\n\n"))
  }
  attr(x,"pdim") <- attr(x,"pvar") <- attr(x,"indexes") <- attr(x,"varnames") <- NULL
  print.table(x)
}
 
  
