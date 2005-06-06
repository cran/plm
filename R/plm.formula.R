plm <- function(y,...){
  UseMethod("plm")
}

plm.formula <-  function(y,instruments=NULL,endog=NULL,data,effect="individual",theta="swar",trinst="baltagi",model=NULL,np=FALSE,...){
  formula <- y
  pvar <- attr(data,"pvar")
  pdim <- attr(data,"pdim")
  indexes <- attr(data,"indexes")
  id.index.name <- indexes$id
  time.index.name <- indexes$time
  name.data <- paste(deparse(substitute(data)))
  formod <- paste(deparse(formula),collapse="")

  if(!is.null(instruments)){
    if(!is.null(endog)){
      forinst <- paste(deparse(instruments[[2]]),collapse="")
      forendog <- paste(deparse(endog[[2]]),collapse="")
      forexog <- paste(deparse(formula[[3]]),collapse="")
      vecinst <- strsplit(forinst,"\\s*\\+\\s*")[[1]]
      vecendog <- strsplit(forendog,"\\s*\\+\\s*")[[1]]
      vecexog <- strsplit(forexog,"\\s*\\+\\s*")[[1]]
      vecinst <- c(vecexog[!vecexog%in%vecendog],vecinst)
      instruments <- as.formula(paste(" ~ ",paste(vecinst,collapse=" + ")))
    }
    forinst <- paste(deparse(instruments[[2]]),collapse="")
    fortot <-  as.formula(paste(formod,"+",forinst,"+",id.index.name,"+",time.index.name))
  }
  else fortot <- as.formula(paste(formod,"+",id.index.name,"+",time.index.name))
  mf <- model.frame(fortot,data=data)
  nr <- nrow(mf)
  id <- mf[[id.index.name]]
  time <- mf[[time.index.name]]

  if(nr!=nrow(data)){
    warning("Sample size has changed\n")
    pdim <- pdimcheck(id,time)
  }

  model.choice <- "oneway"
  if(effect=="double") model.choice <- "double"
  
  if(effect=="double" && (!is.null(instruments) | !is.null(endog) | !pdim$balanced)){
    stop("double effects models not implemented, use time dummies instead\n")
  }

  mtx <- terms(formula)
  X <- model.matrix(mtx,mf)[,-1,drop=F]
  pvar <- pvarcheck(id,time,as.data.frame(X))
  K <- ncol(X)
  y <- model.response(mf)
  pmodel=list(model=model,formula=formula,effect=effect,theta=theta,trinst=trinst)
  if(!is.null(instruments)){
    mtw <- terms(instruments)
    W <- model.matrix(mtw,mf)
  }
  else W <- NULL
  zz <- plm(y,X,W,id,time,pvar,pdim,pmodel)

  if(np){
    z <- nopool(formula,data=data,effect=effect)
    zz$nopool <- z
  }
  zz
}
