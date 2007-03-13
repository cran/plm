plm.extract <- function(formula,data,instruments=NULL,endog=NULL){
  pvar <- attr(data,"pvar")
  pdim <- attr(data,"pdim")
  indexes <- attr(data,"indexes")
  time.names <- pdim$panel.names$time.names
  id.names <- pdim$panel.names$id.names
  id.index.name <- indexes$id
  time.index.name <- indexes$time
  name.data <- paste(deparse(substitute(data)))
  formod <- paste(deparse(formula),collapse="")
  init.instruments <- instruments
  if(!is.null(instruments)){
    if(!is.one.side.formula(instruments)){
      stop("instruments must be a one side formula")
    }
    if(!is.null(endog)){
      if(!is.one.side.formula(endog)){
        stop("endog must be a one side formula")
      }

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
  else{
    fortot <- as.formula(paste(formod,"+",id.index.name,"+",time.index.name))
  }
  mf <- model.frame(fortot,data=data)
  nr <- nrow(mf)
  if(nr!=nrow(data)){
    for (i in names(mf)){
      if (is.factor(mf[[i]])){
        x <- mf[[i]]
        if (length(unique(x)) < length(levels(x))){
          mf[[i]] <- mf[[i]][,drop=TRUE]
        }
      }
    }
  }
  id <- mf[[id.index.name]]
  time <- mf[[time.index.name]]
  mtx <- terms(formula)
  X <- model.matrix(mtx,mf)[,-1,drop=F]
  y <- model.response(mf)
  pdim <- pdim(id,time)
  pvar <- pvar(X,id,time)
  endog <- list(endog=endog,instruments=init.instruments)
  if(!is.null(instruments)){
    mtw <- terms(instruments)
    W <- model.matrix(mtw,mf)
  }
  else W <- NULL
  extr <- list(y=y,X=X,W=W,id=id,time=time)
  extr <- structure(extr,pdim=pdim,pvar=pvar,indexes=indexes,data=name.data)
  extr
}

plm <-  function(formula,data,effect="individual",model=NULL,instruments=NULL,endog=NULL,
                         random.method="swar",inst.method="bvk",...){

  if (!any(class(data) %in% "pdata.frame")){
    stop("argument data should be a pdata.frame\n")
  }

  cl <- match.call()
  model.name <- model

  if(!(effect %in% names(effect.plm.list))){
    stop(paste("effect must be one of ",oneof(effect.plm.list)))
  }
  if (!is.null(model.name)){
    if(!(model.name %in% names(model.plm.list))){
      stop("model must be one of ",oneof(model.plm.list))
    }
  }
  if(!(random.method %in% names(random.method.list))){
    stop(paste("random.method must be one of ",oneof(random.method.list)))
  }
  if(!(inst.method %in% names(inst.method.list))){
    stop(paste("inst.method  must be one of ",oneof(inst.method.list)))
  }
  init.instruments <- instruments

  if(effect=="twoways") model.choice <- "twoways" else model.choice <- "oneway"
  if(effect=="twoways" && (!is.null(instruments) | !is.null(endog) )){
    stop("twoways effects models not implemented, use time dummies instead\n")
  }
#  endog <- list(endog=endog,instruments=init.instruments)
  pmodel <- list(model.name=model.name,formula=formula,effect=effect,
                 random.method=random.method,inst.method=inst.method,
                 instruments=instruments,endog=endog)

  extr <- plm.extract(formula,data,instruments,endog)
  X <- extr$X
  y <- extr$y
  W <- extr$W
  id <- extr$id
  time <- extr$time
  pdim <- attr(extr,"pdim")
  pvar <- attr(extr,"pvar")
  indexes <- attr(extr,"indexes")
  if (is.null(model.name)){
    within <- plm.within(y,X,W,id,time,pvar,pdim,pmodel,indexes,cl,...)
    pooling <- plm.pooling(y,X,W,id,time,pvar,pdim,pmodel,indexes,cl,...)
    between <- plm.between(y,X,W,id,time,pvar,pdim,pmodel,indexes,cl,...)
    random <- plm.random(y,X,W,id,time,pvar,pdim,pmodel,indexes,cl,...)

    attr(within,"pmodel")$model <- "within"
    attr(pooling,"pmodel")$model <- "pooling"
    attr(random,"pmodel")$model <- "random"
    attr(between,"pmodel")$model <- "between"
    
    result <- list(pooling=pooling,within=within,between=between,random=random)
    result <- structure(result,pdim=pdim,pmodel=pmodel,pvar=pvar)
    class(result) <- "plms"
  }
  else{
    result <- switch(model.name,
                     "within"=plm.within(y,X,W,id,time,pvar,pdim,pmodel,indexes,cl,...),
                     "between"=plm.between(y,X,W,id,time,pvar,pdim,pmodel,indexes,cl,...),
                     "pooling"=plm.pooling(y,X,W,id,time,pvar,pdim,pmodel,indexes,cl,...),
                     "random"=plm.random(y,X,W,id,time,pvar,pdim,pmodel,indexes,cl,...),
                     "ht"=plm.ht(y,X,W,id,time,pvar,pdim,pmodel,indexes,cl,...),
                     )
  }
  result$call <- cl
  result
}
