plm <-  function(formula,instruments=NULL,endog=NULL,data,effect="individual",
                         random.method="swar",inst.method="bvk",model=NULL,np=FALSE,...){

  cl <- match.call()
  model.name <- model
  if(!(effect %in% names(effect.list))){
    stop(paste("effect must be one of",oneof(effect.list)))
  }

  if (!is.null(model.name)){
    if(!(model.name %in% names(model.list))){
      stop("model must be one of",oneof(effect.list))
    }
  }

  if(!(random.method %in% names(random.method.list))){
    stop(paste("random.method must be one of",oneof(random.method.list)))
  }

  if(!(inst.method %in% names(inst.method.list))){
    stop(paste("inst.method  must be one of",oneof(inst.method.list)))
  }

  pvar <- attr(data,"pvar")
  pdim <- attr(data,"pdim")
  indexes <- attr(data,"indexes")
  id.index.name <- indexes$id
  time.index.name <- indexes$time
  name.data <- paste(deparse(substitute(data)))
  formod <- paste(deparse(formula),collapse="")
  init.instruments <- instruments
  
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
    for (i in names(data)){
      x <- data[[i]]
      if (length(unique(x)) < length(levels(x))){
        data[[i]] <- data[[i]][,drop=TRUE]
      }
      pdim <- pdim(id,time)
    }
  }
    
  if(effect=="twoways") model.choice <- "twoways" else model.choice <- "oneway"
  if(effect=="twoways" && (!is.null(instruments) | !is.null(endog) | !pdim$balanced)){
    stop("twoways effects models not implemented, use time dummies instead\n")
  }

  mtx <- terms(formula)
  X <- model.matrix(mtx,mf)[,-1,drop=F]
  pvar <- pvar(X,id,time)
  K <- ncol(X)
  y <- model.response(mf)

  endog <- list(endog=endog,instruments=init.instruments)
  pmodel <- list(model.name=model.name,formula=formula,effect=effect,
                 random.method=random.method,inst.method=inst.method,
                 instruments=instruments,endog=endog)
  if(!is.null(instruments)){
    mtw <- terms(instruments)
    W <- model.matrix(mtw,mf)
  }
  else W <- NULL
  zz <- plm.fit(y,X,W,id,time,pvar,pdim,pmodel,indexes,cl)
  zz$call <- cl
    if(np){
    z <- nopool(formula,data=data,effect=effect)
    zz$nopool <- z
  }
    zz
}
