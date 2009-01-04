plm <-  function(formula, data, subset, na.action,
                 effect=c('individual','time','twoways'),
                 model = c('within','random','ht','between','pooling','fd'),
                 instruments = NULL,
                 random.method = c('swar','walhus','amemiya','nerlove'),
                 inst.method = c('bvk','baltagi'),
                 index = NULL, pvar = TRUE, ...){

  # check and match the arguments
  effect <- match.arg(effect)
  model.name <- match.arg(model)
  random.method <- match.arg(random.method)
  inst.method <- match.arg(inst.method)

  if (!is.null(instruments)) formula <- as.Formula(formula,instruments)
  else formula <- Formula(formula)
  
  if (length(formula) == 2){
    iv.estimation <- TRUE
    instruments <- formula(formula,part="second",response=FALSE)
    complete.formula <- formula(formula,part="both")
    formula <- formula(formula)
    instruments <- update(formula,instruments)
  }
  else{
    iv.estimation <- FALSE
    complete.formula <- formula
  }
  
  # force pvar to be computed if hausman taylor model is chosen
  if (model.name=="ht") pvar <- TRUE
  if(effect=="twoways") model.choice <- "twoways" else model.choice <- "oneway"

  # eval the pdata.frame
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
  
  # eval the model.frame
  cl <- match.call()
  mf <- match.call(expand.dots=FALSE)
  m <- match(c("formula","data","subset","na.action"),names(mf),0)
  mf <- mf[c(1,m)]
  mf$drop.unused.levels <- TRUE
  mf$data <- as.name(new.data.name)
  mf[[1]] <- as.name("model.frame")
  mindexes <- mf
  mf$formula <- complete.formula
  mf <- eval(mf,sys.frame(which=nframe))
  mindexes[["formula"]] <- formula(paste("~",id.name,"+",time.name,sep="",collapse=""))
  mindexes <- eval(mindexes,sys.frame(which=nframe))
  int.row.names <- intersect(attr(mf,"row.names"),
                             attr(mindexes,"row.names"))
  mf <- mf[as.character(int.row.names),,drop=F]
  for (i in names(mf)) if(is.factor(mf[[i]])) mf[[i]] <- mf[[i]][drop=TRUE]

  mindexes <- mindexes[as.character(int.row.names),,drop=F]
  attr(mf,"row.names") <- attr(mindexes,"row.names") <- int.row.names

  y <- model.response(mf,"numeric")
  X <- model.matrix(formula,mf)
  if (attr(terms(formula),'intercept')==1){
    X <- X[,-1,drop=FALSE]
    has.intercept <- TRUE
  }
  else{
    has.intercept <- FALSE
  }
  if (iv.estimation) W <- model.matrix(instruments,mf) else W <- NULL
  
  if (model.name == "ht"){
    l <- lapply(mf,function(x) if (is.factor(x)) levels(x))
    l <- l[!sapply(l,is.null)]
    v <- c()
    for (i in 1:length(l)){
      nli <- names(l)[[i]]
      r <- paste(nli,l[[i]],sep="")
      ov <- rep(nli,length(r))
      names(ov) <- r
      v <- c(v,ov)
    }
    attr(X,"var.effects") <- v
  }
  id <- mindexes[[id.name]][drop=T]
  time <- mindexes[[time.name]][drop=T]
  pmodel <- list(model.name = model.name,
                 formula = formula,
                 effect = effect,
                 random.method = random.method,
                 inst.method = inst.method,
                 instruments = instruments,
                 has.intercept = has.intercept)
  
  pdim <- pdim(id,time)
  if (pvar) pvar <- pvar(X,id,time)
  else{
    id.variation <- time.variation <- rep(TRUE,ncol(X))
    names(id.variation) <- names(time.variation) <- colnames(X)
    pvar <- structure(list(id.variation=id.variation,time.variation=time.variation),class="pvar")
  }

  result <- switch(model.name,
                   "within" =plm.within (y,X,W,id,time,pvar,pdim,pmodel,indexes,cl,...),
                   "between"=plm.between(y,X,W,id,time,pvar,pdim,pmodel,indexes,cl,...),
                   "pooling"=plm.pooling(y,X,W,id,time,pvar,pdim,pmodel,indexes,cl,...),
                   "random" =plm.random (y,X,W,id,time,pvar,pdim,pmodel,indexes,cl,...),
                   "ht"     =plm.ht     (y,X,W,id,time,pvar,pdim,pmodel,indexes,cl,...),
                   "fd"     =plm.fd     (y,X,W,id,time,pvar,pdim,pmodel,indexes,cl,...)
                   )
  result$call <- cl
  result$indexes <- list(id=id,time=time)
  result
}

