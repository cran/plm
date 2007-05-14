pread.table <- function(file,id,time=NULL,name=NULL,...){
  x <- read.table(file,...)
  if (is.null(name)){
#    name <- gsub("([^(\.]*)\.(.*)","\\1",file)
    vecname <- strsplit(file,"/")[[1]]
    vecname <- vecname[length(vecname)]
    vecname <- strsplit(vecname,"\\.")[[1]]
#    name <- vecname[length(vecname)]
    name <- vecname[1]
  }
  pdata.frame(x,id=id,time=time,name=name)
}

pdata.frame <- function(x,id,time=NULL,name=NULL){
  id.name <- id
  time.name <- time
  if (is.null(name)){
    data.name <- paste(deparse(substitute(x)))
  }
  else{
    data.name <- name
  }
  # coerce character vectors to factors
  x.char <- names(x)[sapply(x,is.character)]
    for (i in x.char){
    x[[i]] <- factor(x[[i]])
  }
  # check and remove complete NA series
  na.check <- sapply(x,function(x) sum(!is.na(x))==0)
  na.serie <- names(x)[na.check]
  if (length(na.serie)>0){
    if (length(na.serie)==1){
      cat(paste("serie",na.serie," is NA and has been removed\n"))
    }
    else{
      cat(paste("series ",paste(na.serie,collapse=",")," are NA and have been removed\n"))
    }
  }
  x <- x[,!na.check]
  # check and remove cst series
  cst.check <- sapply(x,function(x) myvar(as.numeric(x))==0)
  cst.serie <- names(x)[cst.check]
  if (length(cst.serie)>0){
    if (length(cst.serie)==1){
      cat(paste("serie",cst.serie," is constant and has been removed\n"))
    }
    else{
      cat(paste("series ",paste(na.serie,collapse=",")," are constants and have been removed\n"))
    }
  }
  x <- x[,!cst.check]
  
  if(is.numeric(id.name)){
    if(!is.null(time.name)){warning("The time argument will be ignored\n")}
    N <- nrow(x)
    if( (N%%id.name)!=0){
      stop("unbalanced panel, the id variable should be indicated\n")
    }
    else{
      T <- N%/%id.name
      n <- N%/%T
      time <- rep(1:T,n)
      id <- rep(seq(1:n),rep(T,n))
      id.name <- "id"
      time.name <- "time"
      x[[id.name]] <- id <- as.factor(id)
      x[[time.name]] <- time <- as.factor(time)
    }
  }
  else{
    id <- x[[id.name]] <- as.factor(x[[id.name]])
    if (is.null(time.name)){
      Ti <- table(id)
      n <- length(Ti)
      time <- c()
      for (i in 1:n){
        time <- c(time,1:Ti[i])
      }
      time.name <- "time"
      x[[time.name]] <- time <- as.factor(time)
    }
    else{
      x[[time.name]] <- time <- as.factor(x[[time.name]])
    }
  }
  x <- x[order(id,time),]
  indexes <- list(id=id.name,time=time.name)
  class(indexes) <- "indexes"
  var.names <- names(x)
  pdim <- pdim(x,id.name,time.name)

  for (i in names(x)){
    attr(x[[i]],"data") <- data.name
    oldclass <- attr(x[[i]],"class")
    class(x[[i]]) <- c("pserie",oldclass)
    if (length(unique(x[[i]])) < length(levels(x[[i]]))){
      x[[i]] <- x[[i]][,drop=TRUE]
    }
  }
  pvar <- pvar(x,id.name,time.name)
  x <- structure(x,class=c("pdata.frame","data.frame"),pvar=pvar,pdim=pdim,indexes=indexes)
  assign(data.name,x,pos=1)
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
  centre("Indexes")
  print(indexes)
  centre("Panel Dimensions")
  print(pdim)
  if (sum(!pvar$id.variation)!=1 | sum(!pvar$time.variation)!=1){
    centre("Time/Individual Variation")
    print(pvar,indexes)
  }
  centre("Descriptive Statistics")
  attr(x,"pdim") <- attr(x,"pvar") <- attr(x,"indexes") <- attr(x,"varnames") <- NULL
  print.table(x)
}
 
as.data.frame.pdata.frame <- function(x,...){
  rn <- attr(x,"row.names")
  attr(x,"pdim") <- attr(x,"pvar") <- attr(x,"indexes") <- NULL
  x <- lapply(x,function(x){
    if (length(class(x))==1){
      attr(x,"class") <- NULL
    }
    else{
      class(x) <- class(x)[-1]
    }
    attr(x,"data") <- NULL
    x
  }
              )
  attr(x,"row.names") <- rn
  class(x) <- "data.frame"
  x
}

pfix <- function(x,...){
  subx <- substitute(x)
  if (is.name(subx)){
    subx <- deparse(subx)
  }
  if (!is.character(subx) || length(subx) != 1)
    stop("'fix' requires a name")
  parent <- parent.frame()
  ox <- x
  x <- as.data.frame(x)
  assign(subx, x, env = .GlobalEnv)
  if (exists(subx, envir=parent, inherits = TRUE))
    x <- edit(get(subx, envir=parent), title = subx, ...)
  else {
    x <- edit(function(){}, title = subx, ...)
    environment(x) <- .GlobalEnv
  }
  assign(subx, ox, env = .GlobalEnv)
}
  

"fix" <-
    function (x, ...)
{
    subx <- substitute(x)
    if (is.name(subx)){
        subx <- deparse(subx)
        print(subx)
      }
    if (!is.character(subx) || length(subx) != 1)
        stop("'fix' requires a name")
    parent <- parent.frame()
    if (exists(subx, envir=parent, inherits = TRUE))
        x <- edit(get(subx, envir=parent), title = subx, ...)
    else {
        x <- edit(function(){}, title = subx, ...)
        environment(x) <- .GlobalEnv
    }
    assign(subx, x, env = .GlobalEnv)
}
