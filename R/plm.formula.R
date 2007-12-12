plm <-  function(formula, data, subset, na.action, effect="individual",
                 model = "within", instruments = NULL,
                 random.method = "swar", inst.method = "bvk", index = NULL, pvar = TRUE, ...){
  
  model.name <- model
  if (model.name=="ht") pvar <- TRUE
  data.name <- paste(deparse(substitute(data)))
  new.data.name <- "mydata"
  data2 <- data2plm.data(data,index)
  data <- data2$data
  id.name <- data2$id.name
  time.name <- data2$time.name
  for (i in 1:length(data)){
    attr(data[[i]],"data") <- new.data.name
    attr(data[[i]],"class") <- c("pserie",attr(data[[i]],"class"))
  }
  indexes <- list(id=id.name,time=time.name)
  class(indexes) <- "indexes"
  attr(data,"indexes") <- indexes
  nframe <- length(sys.calls())
  assign(new.data.name,data,env=sys.frame(which=nframe))
   
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
  if(effect=="twoways" && !is.null(instruments)){
    stop("twoways effects models not implemented, use time dummies instead\n")
  }

  osftemp <- ~a
  x <- osftemp
  rhs <- formula[[3]]
  if(length(rhs)>1 && rhs[[1]]=="|"){
    instruments <- osftemp
    instruments[[2]] <- rhs[[3]]
    x[[2]] <- rhs[[2]]
    formula[[3]] <- rhs[[2]]
  }
  else{
    x[[2]] <- rhs
  }
  if(!is.null(instruments)){
    instruments.init <- instruments
    instruments <- update(x,instruments.init)
  }
  cl <- match.call()
  mf <- match.call(expand.dots=FALSE)
  m <- match(c("formula","data","subset","na.action"),names(mf),0)
  mf <- mf[c(1,m)]
  mf$drop.unused.levels <- TRUE
  mf$data <- as.name(new.data.name)
  mf[[1]] <- as.name("model.frame")
  mindexes <- mf
  if(!is.null(instruments)){
    minst <- mf
    minst[["formula"]] <- instruments
  }
  mindexes[["formula"]] <- formula(paste("~",id.name,"+",time.name,sep="",collapse=""))
  mf$formula <- formula
  mf <- eval(mf,sys.frame(which=nframe))
  if(!is.null(instruments)){
    minst <- eval(minst,sys.frame(which=nframe))   
  }
  mindexes <- eval(mindexes,sys.frame(which=nframe))
  y <- model.response(mf,"numeric")
  if(!is.null(instruments)){
    int.row.names <- intersect(attr(mf,"row.names"),
                               intersect(attr(minst,"row.names"),
                                         attr(mindexes,"row.names")))
  }
  else{
    int.row.names <- intersect(attr(mf,"row.names"),
                               attr(mindexes,"row.names"))
  }    
  mf <- mf[int.row.names,]
  mindexes <- mindexes[int.row.names,]
  y <- y[int.row.names]
  attr(mf,"row.names") <- attr(mindexes,"row.names") <- int.row.names
  if(!is.null(instruments)){
    attr(minst,"row.names") <- attr(mindexes,"row.names") <- int.row.names
    minst <- minst[int.row.names,,drop=F]
    W <- model.matrix(instruments,minst)
  }
  else{
    W <- NULL
  }
  cl$formula <- formula
  if (!is.null(instruments)) cl$instruments <- instruments.init
  X <- model.matrix(formula,mf)[,-1,drop=FALSE]
  id <- mindexes[[id.name]][drop=T]
  time <- mindexes[[time.name]][drop=T]
  pmodel <- list(model.name=model.name,formula=formula,effect=effect,
                 random.method=random.method,inst.method=inst.method,
                 instruments=instruments)
  pdim <- pdim(id,time)
  if (pvar){
    pvar <- pvar(X,id,time)
  }
  else{
    id.variation <- time.variation <- rep(TRUE,ncol(X))
    names(id.variation) <- names(time.variation) <- colnames(X)
    pvar <- structure(list(id.variation=id.variation,time.variation=time.variation),class="pvar")
  }
  result <- switch(model.name,
                   "within"=plm.within(y,X,W,id,time,pvar,pdim,pmodel,indexes,cl,...),
                   "between"=plm.between(y,X,W,id,time,pvar,pdim,pmodel,indexes,cl,...),
                   "pooling"=plm.pooling(y,X,W,id,time,pvar,pdim,pmodel,indexes,cl,...),
                   "random"=plm.random(y,X,W,id,time,pvar,pdim,pmodel,indexes,cl,...),
                   "ht"=plm.ht(y,X,W,id,time,pvar,pdim,pmodel,indexes,cl,...),
                   "fd"=plm.fd(y,X,W,id,time,pvar,pdim,pmodel,indexes,cl,...)
                   )
  result$call <- cl
  result$indexes <- list(id=id,time=time)
  result
}
