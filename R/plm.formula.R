plm <-  function(formula, data, subset, na.action,
                 effect=c('individual','time','twoways'),
                 model = c('within','random','ht','between','pooling','fd'),
                 random.method = c('swar','walhus','amemiya','nerlove'),
                 inst.method = c('bvk','baltagi'),
                 index = NULL,
                 ...){

  dots <- list(...)
  if (!is.null(dots$instruments)){
    formula <- as.formula(paste(deparse(formula),"|",deparse(dots$instruments[[2]])))
    deprec.instruments <- paste("the use of the instruments argument is deprecated,",
                                "use two-part formulas instead")
    warning(deprec.instruments)
  }

  # check and match the arguments
  effect <- match.arg(effect)
  if (!any(is.na(model))) model <- match.arg(model)
  random.method <- match.arg(random.method)
  inst.method <- match.arg(inst.method)
  
  # gestion des index (a revoir en relation avec pdata.frame)
  data <- plm.data(data,index)
  data[[1]] <- factor(data[[1]])
  data[[2]] <- factor(data[[2]])
  index <- names(data)[1:2]
  names(index) <- c("id","time")
  new.data.name <- "mydata"
  attr(data,"indexes") <- as.list(index)
  nframe <- length(sys.calls())
  assign(new.data.name,data, env = sys.frame(which = nframe))
  for (i in 1:length(mydata)){
    attr(mydata[[i]], "data") <- new.data.name
    attr(mydata[[i]], "class") <- c("pserie", attr(mydata[[i]], "class"))
  }
  
  ##    the logic of that is that from now we know the name and the data
  ##    and will later on attach it to every variable it contains so that
  ##    the transformation functions (diff, within ...) can be used
  
  # Create a Formula object, with index as extra attributes
  formula <- pFormula(formula, extra = index)
  # in case of 2part formula, check whether the second part should be updated
  if (length(formula) == 2) formula <- expand.formula(formula)
  # eval the model.frame
  cl <- match.call()
  mf <- match.call(expand.dots=FALSE)
  m <- match(c("formula", "data", "subset", "na.action"),names(mf),0)
  mf <- mf[c(1,m)]
  mf$drop.unused.levels <- TRUE
  mf$data <- as.name(new.data.name)
  mf[[1]] <- as.name("model.frame")
  mf$formula <- formula
#  mf$formula <- do.call("pFormula", list(mf$formula, extra = index))
  mf$include.extra <- TRUE
  data <- eval(mf,sys.frame(which = nframe))
  class(data) <- c("pdata.frame", "data.frame")
  if (is.na(model)) attr(data, "formula") <- formula
  if (!is.na(model)){
    result <- switch(model,
                     "within"  = plm.within (formula, data, effect),
                     "between" = plm.between(formula, data, effect),
                     "pooling" = plm.pooling(formula, data),
                     "random"  = plm.random (formula, data, effect, random.method, inst.method),
                     "ht"      = plm.ht     (formula, data),
                     "fd"      = plm.fd     (formula, data)
                     )
    result$call <- cl
  }
  else{
    result <- data
  }
  result
}



