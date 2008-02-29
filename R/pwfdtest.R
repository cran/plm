
pwfdtest<-function(x, data, h0=c("fd","fe"), ...) {
  ## first-difference-based serial correlation test for panel models
  ## ref.: Wooldridge (2003), par. 10.6 
  if(!require(car)) stop("Library 'car' is needed")

  ## fetch fd residuals
  fdmod<-plm(x,data,model="fd")
  FDres<-fdmod$residuals

  ## indices (full length! must reduce by 1st time period)
   ## this is an ad-hoc solution for the fact that the 'fd' model
   ## carries on the full indices while losing the first time period
   time <- as.numeric(fdmod$indexes$time)
   id <- as.numeric(fdmod$indexes$id)


   ## fetch dimensions and adapt to those of indices
   n<-attr(fdmod,"pdim")$nT$n
 
   ## (re)create groupwise-separated index from 1 to nT 
   ## - dropping first time period
   ## - correcting Ti=Ti+1
   Ti<-attr(fdmod,"pdim")$Tint$Ti
   redind<-vector("list",n)
   tfirst<-0
   for(i in 1:n) {
     redind[[i]]<-(tfirst+2):(tfirst+Ti[i]+1)
     tfirst<-max(redind[[i]])
    }

   ## reduce indices by 1st time period
   redind<-unlist(redind)
   time<-time[redind]
   id<-id[redind]
    
  N <- length(FDres)
  FDres.1 <- c(NA,FDres[1:(N-1)])

  lagid <- id-c(NA,id[1:(N-1)])
  FDres.1[lagid!=0] <- NA

  ## make (panel) dataframe for auxiliary regression
  auxdata <- as.data.frame(cbind(id,time))
  auxdata$FDres<-FDres
  auxdata$FDres.1<-FDres.1

  ## pooling model FDres vs. lag(FDres),
  ## with intercept (might as well do it w.o.)
  auxmod<-plm(FDres~FDres.1,auxdata,model="pooling")

  switch(match.arg(h0), 
             fd = {h0des<-"differenced"
                   ## theoretical rho under H0: no serial 
                   ## corr. in differenced errors is 0
                   rho.H0 <- 0},
             fe = {h0des<-"original"
                   ## theoretical rho under H0: no serial 
                   ## corr. in original errors is -0.5
                   rho.H0 <- -0.5})

  ## test H0: rho=rho.H0 with HAC t-test (HC0-3 parm may be passed)
  myvcov<-function(x) pvcovHC(x, method="arellano", ...)

  myH0<-paste("FDres.1 = ", as.character(rho.H0), sep="")

  lhtest<-linear.hypothesis(model=auxmod, myH0, vcov.=myvcov, ...)
  
  ##(insert usual htest features)  
  FDARstat <- lhtest[2,3]
  names(FDARstat) <- dimnames(lhtest)[[2]][3] 
  ## this is either 'F' or 'Chisq' and is the name of 3rd
  ## column because we are supplying a vcov matrix
  pFDAR<-lhtest[2,4]

  dname <- paste(deparse(substitute(x)))
  RVAL <- list(statistic = FDARstat, parameter = NULL,
               method = "Wooldridge's first-difference test for serial correlation in panels",
               alternative = paste("serial correlation in", h0des, "errors"),
               p.value = pFDAR,
               data.name =   dname)
  class(RVAL) <- "htest"
  return(RVAL)

}