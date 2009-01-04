#  pbltest pgtest

#### pbgtest

pbgtest <- function (x, ...) 
{
    UseMethod("pbgtest")
}

pbgtest.formula<-function(x, data, model="random",order=NULL, index=NULL, ...) {
  ## formula method for pbgtest;
  ## defaults to a RE model
  mymod <- plm(formula=x, data=data, model=model, index=index,...)
  pbgtest(mymod,order)
}

pbgtest.panelmodel<-function(x, order = NULL, ...) {
  ## residual serial correlation test based on the residuals of the demeaned
  ## model (see Wooldridge p.288) and the regular bgtest() in {lmtest}


  ## structure:
  ## 1: take demeaned data from 'plm' object
  ## 2: est. auxiliary model by OLS on demeaned data
  ## 3: apply bgtest() to auxiliary model and return the result

  ## retrieve demeaned data
  demX <- x$model[[2]]
  demy <- x$model[[1]]
  
  ## ...and group numerosities
  Ti <- attr(x, "pdim")$Tint$Ti
  
  ## set lag order to minimum group numerosity if not specified by user
  ## (check whether this is sensible)
  
  if(is.null(order)) order <- min(Ti)
  
  ## bg test on the demeaned model:
  
  ## check package availability and load if necessary
  lm.ok <- require("lmtest")
  if(!lm.ok) stop("package lmtest is needed but not available")
  
  ## bgtest is the bgtest, exception made for the method attribute
  bgtest <- bgtest(lm(demy~demX-1),order=order)
  bgtest$method <- "Breusch-Godfrey/Wooldridge test for serial correlation in panel models"
  bgtest$alternative <- "serial correlation in idiosyncratic errors"
  bgtest$data.name <- paste(deparse(x$call$formula))
  names(bgtest$statistic) <- "chisq"
  # is it really working ? I'm wondering how the order of the test is taken into account ???
  return(bgtest)
}

### pwtest

pwtest <- function(x, data, effect = c("individual","time"), index=NULL, ...) {

  ## "RE" test à la Wooldridge, see 10.4.4
  ## (basically the scaled and standardized estimator for sigma from REmod)
  ## does not rely on normality or homoskedasticity; 
  ## H0: composite errors uncorrelated

  ## ref. Wooldridge, p.264

  ######### from here generic testing interface from
  ######### plm to my code

  ## tind is actually not needed here

  data <- plm.data(data,index=index)

  ## extract indices
  indices <- list(uno=as.numeric(data[,1]),
                due=as.numeric(data[,2])
                )  

  ## if effect="individual" std., else swap
  myinds <- switch(match.arg(effect),
          individual=indices,
          time=list(indices[[2]],indices[[1]])
          )

  index <- myinds[[1]]
  tindex <- myinds[[2]]

  ## reduce X,y to model matrix values (no NAs)
  X <- model.matrix(x,data=data)
  y <- model.response(model.frame(x,data=data))
  ## reduce index accordingly
  names(index) <- row.names(data)
  ind <- index[which(names(index)%in%row.names(X))]

  ######### till here. 
  ######### Add ordering here if needed.

  ## this doesn't need ordering of obs.

  ## det. number of groups and df
  n <- length(unique(ind))
  k <- dim(X)[[2]]
  ## det. max. group numerosity
  t <- max(tapply(X[,1],ind,length))
  ## det. total number of obs. (robust vs. unbalanced panels)
  nT <- length(ind)

  ## ref. Wooldridge, p.264
    
    ## estimate pooled model and extract resids
    poolmod <- lm.fit(X,y)
    u <- resid(poolmod)

    ## est. random effect variance
    ## "pre-allocate" an empty list of length n
    tres <- vector("list", n)

    ## list of n "empirical omega-blocks"
    ## with averages of xproducts of t(i) residuals
    ## for each group 1..n 
    ## (possibly different sizes if unbal., thus a list
    ## and thus, unlike Wooldridge (eq.10.37), ve divide 
    ## every block by *his* t(t-1)/2)
    unind <- unique(ind)
    for(i in 1:n) {
      ut <- u[ind==unind[i]]
      tres[[i]] <- ut%o%ut
      }

    ## sum over all upper triangles of emp. omega blocks:
    ## define aux. function
    uptrisum <- function(x) {
              uts <- sum(x[upper.tri(x,diag=FALSE)])
              return(uts)}

    ## det. # of upper triangle members (n*t(t-1)/2 if balanced)
    ti <- sapply(tres, function(x) dim(x)[[1]])
    uptrinum <- sum(ti*(ti-1)/2)  # don't need this!!

    ## ...apply to list and sum over resulting vector (df corrected)
    W <- sum(sapply(tres,uptrisum)) # /sqrt(n) simplifies out

    ## calculate se(Wstat) as in 10.40
    seW <- sqrt( sum( sapply(tres,uptrisum)^2 ) )

    ## NB should we apply a df correction here, maybe that of the standard
    ## RE estimator? (see page 261) 

    Wstat <- W/seW
    names(Wstat) <- "z"
    pW <- 2*pnorm(abs(Wstat),lower.tail=FALSE) # unlike LM, test is two-tailed!

  ##(insert usual htest features)
  dname <- paste(deparse(substitute(formula)))
  RVAL <- list(statistic = Wstat, parameter = NULL,
               method = paste("Wooldridge's test for unobserved ",
                               match.arg(effect),"effects "),
               alternative = "unobserved effect",
               p.value = pW,
               data.name =   dname)
  class(RVAL) <- "htest"
  return(RVAL)

  }


### pwartest

pwartest <- function(x, data, index=NULL, ...) {
  ## small-sample serial correlation test for FE models
  ## ref.: Wooldridge (2003) 10.5.4 
  if(!require(car)) stop("Library 'car' is needed")
  data <- plm.data(data,index=index)
  ## fetch within residuals
  femod <- plm(x,data,model="within")
  FEres <- resid(femod)

  ## this is a bug fix for incorrect naming of the "data" attr.
  ## for the pseries in pdata.frame()
  attr(FEres, "data") <- paste(deparse(substitute(data)))
#  FEres.1 <- lag(FEres,k=1)
  N <- length(FEres)
#  FEres.1 <- c(NA,FEres[2:N]-FEres[1:(N-1)])
  FEres.1 <- c(NA,FEres[1:(N-1)])
  id <- as.numeric(femod$indexes$id)
  lagid <- id-c(NA,id[1:(N-1)])
  FEres.1[lagid!=0] <- NA
  data$FEres <- FEres
  data$FEres.1 <- FEres.1
  ## pooling model FEres vs. lag(FEres)
  auxmod <- plm(FEres~FEres.1,data,model="pooling")

  ## calc. theoretical rho under H0: no serial corr. in errors
#  t.<-attr(data, "pdim")$nT$T
  t. <- pdim(data)$nT$T
  rho.H0 <- -1/(t.-1)

  ## test H0: rho=rho.H0 with HAC t-test (HC0-3 parm may be passed)
  myvcov <- function(x) pvcovHC(x, method="arellano", ...)
  myH0 <- paste("FEres.1 = ", as.character(rho.H0), sep="")
  linear.hypothesis(model=auxmod,"FEres.1=-0.125")
  lhtest <- linear.hypothesis(model=auxmod, myH0, vcov.=myvcov, ...)
#  print(lhtest)
#  print(lhtest$Df)
#  print(names(lhtest))
#  print.default(lhtest)
  ##(insert usual htest features)  
  FEARstat <- lhtest[2,3]
  names(FEARstat) <- dimnames(lhtest)[[2]][3]
  if (names(FEARstat)=="Chisq") names(FEARstat) <- "chisq"
  ## this is either 'F' or 'Chisq' and is the name of 3rd
  ## column because we are supplying a vcov matrix
  pFEAR <- lhtest[2,4]

  dname <- paste(deparse(substitute(x)))
  RVAL <- list(statistic = FEARstat,
               parameter = NULL,
               method = "Wooldridge's test for serial correlation in FE panels",
               alternative = "serial correlation",
               p.value = pFEAR,
               data.name =   dname)
  class(RVAL) <- "htest"
  return(RVAL)

}

### pbsytest

pbsytest <- function (x, ...){
  UseMethod("pbsytest")
}


pbsytest.formula <- function(x,data,test=c("AR","RE","J"), index=NULL, ...) {

  ## Bera., Sosa-Escudero and Yoon type LM test for random effects
  ## under serial correlation (H0: no random effects) or the inverse;
  ## ref. Baltagi 2005, pages 96-97;
  ## original ref. Bera, Sosa-Escudero and Yoon, JE 101 (2001)
  ## test="AR" you get the serial corr. test robust vs. RE
  ## test="RE" you get the RE test robust vs. serial corr.

  ######### from here generic testing interface from
  ######### plm to my code
  formula <- x

  data <- plm.data(data,index=index)

  ## extract indices
  index <- as.numeric(data[,1])
  tindex <- as.numeric(data[,2])

  ## reduce X,y to model matrix values (no NAs)
  X <- model.matrix(formula,data=data)
  y <- model.response(model.frame(formula,data=data))
  ## reduce index accordingly
  names(index) <- row.names(data)
  ind <- index[which(names(index)%in%row.names(X))]
  names(tindex) <- row.names(data)  
  tind <- index[which(names(tindex)%in%row.names(X))]

######### till here. 
  ######### Add ordering here if needed.

  ## this needs ordering of obs. on time, regardless 
  ## whether before that on groups or after

  ######### Ordering and numerosity check

  ## order by group, then time
  oo <- order(ind,tind)
  X <- X[oo,]
  y <- y[oo]
  ind <- ind[oo]
  tind <- tind[oo]

  ## det. number of groups and df
  n <- length(unique(ind))
  k <- dim(X)[[2]]
  ## det. max. group numerosity
  t <- max(tapply(X[,1],ind,length))
  ## det. total number of obs. (robust vs. unbalanced panels)
  nT <- length(ind)

  ######### End ordering block

  ## begin test
    
            ## estimate pooled model
  poolmod <- lm.fit(X,y)

            ## extract pooled res. for BP statistic
  poolres <- resid(poolmod)

            ## calc. A and B:
  S1 <- sum( tapply(poolres,ind,sum)^2 )
  S2 <- sum( poolres^2 )
            
  A <- S1/S2-1
		
  unind <- unique(ind)
  uu <- rep(NA,length(unind))
  uu1 <- rep(NA,length(unind))
  for(i in 1:length(unind)) {
    u.t <- poolres[ind==unind[i]]
    u.t.1 <- u.t[-length(u.t)]
    u.t <- u.t[-1]
    uu[i] <- crossprod(u.t)
    uu1[i] <- crossprod(u.t,u.t.1)
  }
  
  B <- sum(uu1)/sum(uu)
  
  switch(match.arg(test),
         AR={LM <- (n * t^2 * (B - (A/t))^2) / ((t-1)*(1-(2/t)))
             df <- c(df=1)
#                  names(LM) <- paste("Chisq(", df, ")", sep="")
             names(LM) <- "chisq"
             pLM <- pchisq(LM,df=1,lower.tail=FALSE)
             tname <- "Bera, Sosa-Escudero and Yoon locally robust test"
             myH0 <- "AR(1) errors sub random effects"
           },
         RE={LM <- (A - 2*B) * sqrt( (n * t) / (2*(t-1)*(1-(2/t))) )
#                  names(LM) <- "Z"
             names(LM) <- "z"
             df <- NULL
             pLM <- pnorm(LM,lower.tail=FALSE)
             tname <- "Bera, Sosa-Escudero and Yoon locally robust test"
             myH0 <- "random effects sub AR(1) errors"
           },              
         J={LM <- (n * t^2) / (2*(t-1)*(t-2)) * (A^2 - 4*A*B + 2*t*B^2) 
            df <- c(df=2)
#                 names(LM) <- paste("Chisq(", df, ")", sep="")
            names(LM) <- "chisq"
            pLM <- pchisq(LM,df=1,lower.tail=FALSE)
            tname <- "Baltagi and Li AR-RE joint test"
            myH0 <- "AR(1) errors or random effects"
          }
         )

  dname <- paste(deparse(substitute(formula)))
  RVAL <- list(statistic = LM,
               parameter = df,
               method = tname,
               alternative = myH0,
               p.value = pLM,
               data.name =   dname)
  class(RVAL) <- "htest"
  return(RVAL)

  }

### pdwtest

pdwtest <- function (x, ...) 
{
    UseMethod("pdwtest")
}

pdwtest.formula <- function(x, data, model="random", index=NULL, ...) {
  ## formula method for pdwtest;
  ## defaults to a RE model
  mymod <- plm(formula=x, data=data, model=model, index=index, ...)
  pdwtest(mymod, ...)
}

pdwtest.panelmodel <- function(x,...) {
  ## residual serial correlation test based on the residuals of the demeaned
  ## model and the regular dwtest() in {lmtest}
  ## reference Baltagi (page 98) for FE application, Wooldridge page 288 for
  ## the general idea.


  ## structure:
  ## 1: take demeaned data from 'plm' object
  ## 2: est. auxiliary model by OLS on demeaned data
  ## 3: apply bgtest() to auxiliary model and return the result

  ## retrieve demeaned data
  demX <- x$model[[2]]
  demy <- as.numeric(x$model[[1]])  ## pseries' give problems

 ## dw test on the demeaned model:

 ## check package availability and load if necessary
  lm.ok <- require("lmtest")
  if(!lm.ok) stop("package lmtest is needed but not available")
  
 ## ARtest is the bgtest, exception made for the method attribute
  ARtest <- dwtest(lm(demy~demX-1),...)
  ARtest$method <- "Durbin-Watson test for serial correlation in panel models"
  ARtest$alternative <- "serial correlation in idiosyncratic errors"
  ARtest$data.name <- paste(deparse(x$call$formula))
#  names(ARtest$statistic) <- "chisq"
  return(ARtest)
}


### pbltest


######### Baltagi and Li's LM_rho|mu ########
## ex Baltagi and Li (1995) Testing AR(1) against MA(1)...,
## JE 68, 133-151, test statistic (one-sided) is LM_4;
## see also idem (1997), Monte carlo results..., 
## Annales d'Econometrie et Statistique 48, formula (8)

## from version 2: disposes of Kronecker products,
## thus much faster and feasible on large NT (original 
## is already infeasible for NT>3000, this takes 10'' 
## on N=3000, T=10 and even 20000x10 (55'') is no problem;
## lme() hits the memory limit at ca. 20000x20)

pbltest <- function(x,data,alternative=c("twosided","onesided"), index=NULL, ...) {
 ## TODO: lme fails if there are any NAs in the data: reduce!

 ## usage: pbltest(modd.reg,data=dati0)

 ## this version (pbltest0) based on a "formula, pdataframe" interface

  data <- plm.data(data,index=index)
  gindex <- names(data)[1]
  tindex <- names(data)[2]

  ## sort according to indices
  eval(parse(text=paste("data <- data[order(data$",gindex,",data$",tindex,"),]",sep="")))

  ## make random effects formula
  rformula <- NULL
  eval(parse(text=paste("rformula <- ~1|",gindex,sep="")))

  require(nlme)

  ## est. MLE model
  mymod <- lme(x,data=data,random=rformula,method="ML")

  nt. <- mymod$dims$N
  n. <- as.numeric(mymod$dims$ngrps[1])
  t. <- nt./n.
  Jt <- matrix(1,ncol=t.,nrow=t.)/t.
  Et <- diag(1,t.)-Jt
  ## make 'bidiagonal' matrix (see BL, p.136)
  G <- matrix(0,ncol=t.,nrow=t.)
  for(i in 2:t.) {
    G[i-1,i] <- 1
    G[i,i-1] <- 1
    }
  
  ## retrieve composite (=lowest level) residuals
  uhat <- residuals(mymod,level=0)

  ## sigma2.e and sigma2.1 as in BL
  ## break up residuals by group to get rid of Kronecker prod.
  ## data have to be balanced and sorted by group/time, so this works
  uhat.i <- vector("list",n.)
  for(i in 1:n.) {
    uhat.i[[i]] <- uhat[t.*(i-1)+1:t.]
    } 
  s2e <- rep(NA,n.)
  s21 <- rep(NA,n.)
  for(i in 1:n.) {
    u.i <- uhat.i[[i]]
    s2e[i] <- as.numeric(crossprod(u.i,Et) %*% u.i)
    s21[i] <- as.numeric(crossprod(u.i,Jt) %*% u.i)
    }
  sigma2.e <- sum(s2e) / (n.*(t.-1)) 
  sigma2.1 <- sum(s21) / n. 

  ## calc. score under the null:
  star1 <- (Jt/sigma2.1 + Et/sigma2.e) %*% G %*% (Jt/sigma2.1 + Et/sigma2.e)
  star2 <- rep(NA,n.)
  ## again, do this group by group to avoid Kronecker prod.
  for(i in 1:n.) {
    star2[i] <- as.numeric(crossprod(uhat.i[[i]],star1) %*% uhat.i[[i]])
    }
  star2 <- sum(star2)
  Drho <- (n.*(t.-1)/t.) * (sigma2.1-sigma2.e)/sigma2.1 + sigma2.e/2 * star2
  ## star2 is (crossprod(uhat, kronecker(In, star1)) %*% uhat)
  
  ## components for the information matrix
  a <- (sigma2.e-sigma2.1)/(t.*sigma2.1)
  j.rr <- n. * (2 * a^2 * (t.-1)^2 + 2*a*(2*t.-3) + (t.-1))
  j.12 <- n.*(t.-1)*sigma2.e / sigma2.1^2
  j.13 <- n.*(t.-1)/t. * sigma2.e * (1/sigma2.1^2 - 1/sigma2.e^2)
  j.22 <- (n. * t.^2) / (2 * sigma2.1^2)
  j.23 <- (n. * t.) / (2 * sigma2.1^2)
  j.33 <- (n./2) * (1/sigma2.1^2 + (t.-1)/sigma2.e^2)
  
  ## build up information matrix
  Jmat <- matrix(nrow=3,ncol=3)
  Jmat[1,] <- c(j.rr,j.12,j.13)
  Jmat[2,] <- c(j.12,j.22,j.23)
  Jmat[3,] <- c(j.13,j.23,j.33)
  
  J11 <- n.^2 * t.^2 * (t.-1) / (det(Jmat) * 4*sigma2.1^2 * sigma2.e^2)
  ## this is the same as J11 <- solve(Jmat)[1,1], see BL page 73

  switch(match.arg(alternative),
         onesided = {
           LMr.m <- Drho * sqrt(J11)
           pval <- pnorm(LMr.m,lower.tail=FALSE)
                                        #    names(LMr.m) <- "Z"
           names(LMr.m) <- "z"
           method1 <- "one-sided"
           method2 <- "H0: rho = 0, HA: rho > 0"
           parameter <- NULL
         },
         twosided = {
           LMr.m <- Drho^2 * J11
           pval <- pchisq(LMr.m,1,lower.tail=FALSE)
#    names(LMr.m) <- "Chisq(1)"
           names(LMr.m) <- "chisq"
           parameter <- c(df=1)
           method1 <- "two-sided"
           method2 <- "H0: rho = 0, HA: rho != 0" 
         }
         )
  dname <- paste(deparse(substitute(x)))
  method <- paste("Baltagi and Li", method1,"LM test")
  alternative <- "AR(1)/MA(1) errors in RE panel models"

  res <- list(statistic = LMr.m,
              p.value = pval,
              method = method,
              alternative = alternative,
              parameter = parameter,
              data.name = dname)

  class(res) <- "htest"
  res
}

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
