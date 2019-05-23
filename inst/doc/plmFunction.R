## ----setup, echo=FALSE---------------------------------------------------
library("knitr")
opts_chunk$set(message = FALSE, warning = FALSE)

## ----texreg, echo = FALSE, results = "hide"------------------------------
library("texreg")
extract.plm <- function(model, include.rsquared = TRUE, include.adjrs = TRUE, 
    include.nobs = TRUE, include.ercomp = TRUE, ...) {
  s <- summary(model, ...)
  coefficient.names <- rownames(coef(s))
  coefficients <- coef(s)[, 1]
  standard.errors <- coef(s)[, 2]
  significance <- coef(s)[, 4]
  
  rs <- s$r.squared[1]
  adj <- s$r.squared[2]
  n <- length(model$residuals)
  
  gof <- numeric()
  gof.names <- character()
  gof.decimal <- logical()

  if (include.ercomp == TRUE){
      if (model$args$model == "random"){
          se <- sqrt(ercomp(model)$sigma)
          gof <- c(gof, se)
          gof.names <- c(gof.names, paste("s_", names(se), sep = ""))
          gof.decimal <- c(gof.decimal, rep(TRUE, length(se)))
      }
  }  
  if (include.rsquared == TRUE) {
    gof <- c(gof, rs)
    gof.names <- c(gof.names, "R$^2$")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  if (include.adjrs == TRUE) {
    gof <- c(gof, adj)
    gof.names <- c(gof.names, "Adj.\ R$^2$")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  if (include.nobs == TRUE) {
    gof <- c(gof, n)
    gof.names <- c(gof.names, "Num.\ obs.")
    gof.decimal <- c(gof.decimal, FALSE)
  }  
  tr <- createTexreg(
      coef.names = coefficient.names, 
      coef = coefficients, 
      se = standard.errors, 
      pvalues = significance, 
      gof.names = gof.names, 
      gof = gof, 
      gof.decimal = gof.decimal
  )
  return(tr)
}

setMethod("extract", signature = className("plm", "plm"), 
    definition = extract.plm)


## ----grunfeld------------------------------------------------------------
library("plm")
data("Grunfeld", package = "plm")
ols <- plm(inv ~ value + capital, Grunfeld, model = "pooling")
between <- update(ols, model = "between")
within <- update(ols, model = "within")
walhus <- update(ols, model = "random", random.method = "walhus", random.dfcor = 3)
amemiya <- update(walhus, random.method = "amemiya")
swar <- update(amemiya, random.method = "swar")

## ----grunfeldresults, echo = TRUE----------------------------------------
library("texreg")
screenreg(list(ols = ols, between = between, within = within, 
            walhus = walhus, amemiya = amemiya, swar = swar),
        digits = 5, omit.coef = "(Intercept)")

## ----ercompamemiya-------------------------------------------------------
ercomp(amemiya)

## ----produc--------------------------------------------------------------
data("Produc", package = "plm")
PrSwar <- plm(log(gsp) ~ log(pcap) + log(pc) + log(emp) + unemp, Produc, 
           model = "random", random.method = "swar", random.dfcor = 3)
summary(PrSwar)

## ----grunfeld2ways-------------------------------------------------------
Grw <- plm(inv ~ value + capital, Grunfeld, model = "random", effect = "twoways", 
           random.method = "walhus", random.dfcor = 3)
Grs <- update(Grw, random.method = "swar")
Gra <- update(Grw, random.method = "amemiya")
screenreg(list("Wallace-Hussain" = Grw, "Swamy-Arora" = Grs, "Amemiya" = Gra), digits = 5)

## ----produc2ways---------------------------------------------------------
data("Produc", package = "plm")
Prw <- plm(log(gsp) ~ log(pcap) + log(pc) + log(emp) + unemp, Produc, 
           model = "random", random.method = "walhus", 
           effect = "twoways", random.dfcor = 3)
Prs <- update(Prw, random.method = "swar")
Pra <- update(Prw, random.method = "amemiya")
screenreg(list("Wallace-Hussain" = Prw, "Swamy-Arora" = Prs, "Amemiya" = Pra), digits = 5)

## ----hedonic-------------------------------------------------------------
data("Hedonic", package = "plm")
form <- mv ~ crim + zn + indus + chas + nox + rm + 
    age + dis + rad + tax + ptratio + blacks + lstat
HedStata <- plm(form, Hedonic, model = "random", index = "townid", 
                random.models = c("within", "between"))
HedEviews <- plm(form, Hedonic, model = "random", index = "townid", 
                 random.models = c("within", "Between"))
HedEviewsWH <- update(HedEviews, random.models = "pooling")
screenreg(list(EViews = HedEviews, Stata = HedStata, "Wallace-Hussain" = HedEviewsWH), 
          digits = 5, single.row = TRUE)

## ------------------------------------------------------------------------
data("Crime", package = "plm")
crbalt <- plm(lcrmrte ~ lprbarr + lpolpc + lprbconv + lprbpris + lavgsen +
              ldensity + lwcon + lwtuc + lwtrd + lwfir + lwser + lwmfg + lwfed +
              lwsta + lwloc + lpctymle + lpctmin + region + smsa + factor(year)
              | . - lprbarr - lpolpc + ltaxpc + lmix,
              data = Crime, model = "random", inst.method = "baltagi")
crbvk <- update(crbalt, inst.method = "bvk")
crwth <- update(crbalt, model = "within")
screenreg(list(FE2SLS = crwth, EC2SLS = crbalt, G2SLS = crbvk), 
          single.row = TRUE, digits = 5, omit.coef = "(region)|(year)",
          reorder.coef = c(1:16, 19, 18, 17))

## ------------------------------------------------------------------------
data("Wages", package = "plm")
ht <- plm(lwage ~ wks + south + smsa + married + exp + I(exp^2) + 
            bluecol + ind + union + sex + black + ed | 
            bluecol + south + smsa + ind + sex + black | 
            wks + married + exp + I(exp^2) + union, 
          data = Wages, index = 595, 
          inst.method = "baltagi", model = "random", 
          random.method = "ht")

am <- update(ht, inst.method = "am")
  
screenreg(list("Hausman-Taylor" = ht, "Amemiya-MaCurdy" = am), 
          digits = 5, single.row = TRUE)

## ------------------------------------------------------------------------
data("Produc", package = "plm")
swar <- plm(form <- log(gsp) ~ log(pc) + log(emp) + log(hwy) + log(water) + log(util) + unemp, 
            Produc, index = c("state", "year", "region"), effect = "nested", random.method = "swar")
walhus <- update(swar, random.method = "walhus")
amem <- update(swar, random.method = "amemiya")
screenreg(list("Swamy-Arora" = swar, "Wallace-Hussain" = walhus, "Amemiya" = amem), digits = 5)

