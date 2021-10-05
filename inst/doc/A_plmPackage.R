## ----echo=FALSE,results='hide'------------------------------------------------
options(prompt= "R> ", useFancyQuotes = FALSE, scipen = 999)
library("knitr")
opts_chunk$set(message = FALSE, warning = FALSE)

## ----echo=TRUE, results='hide'------------------------------------------------
library("plm")

## -----------------------------------------------------------------------------
data("EmplUK", package="plm")
data("Produc", package="plm")
data("Grunfeld", package="plm")
data("Wages", package="plm")

## ----setdata1-----------------------------------------------------------------
head(Grunfeld)
E <- pdata.frame(EmplUK, index=c("firm","year"), drop.index=TRUE, row.names=TRUE)
head(E)
head(attr(E, "index"))

## -----------------------------------------------------------------------------
summary(E$emp)
head(as.matrix(E$emp))

## -----------------------------------------------------------------------------
head(lag(E$emp, 0:2))

## -----------------------------------------------------------------------------
head(diff(E$emp), 10)
head(lag(E$emp, 2), 10)
head(Within(E$emp))
head(between(E$emp), 4)
head(Between(E$emp), 10)

## ----results='hide'-----------------------------------------------------------
emp ~ wage + capital | lag(wage, 1) + capital
emp ~ wage + capital | . -wage + lag(wage, 1)

## ----fe_re--------------------------------------------------------------------
grun.fe <- plm(inv~value+capital, data = Grunfeld, model = "within")
grun.re <- plm(inv~value+capital, data = Grunfeld, model = "random")

## ----summary_re---------------------------------------------------------------
summary(grun.re)
ranef(grun.re)

## -----------------------------------------------------------------------------
fixef(grun.fe, type = "dmean")

## -----------------------------------------------------------------------------
summary(fixef(grun.fe, type = "dmean"))

## -----------------------------------------------------------------------------
grun.twfe <- plm(inv~value+capital, data=Grunfeld, model="within", effect="twoways")
fixef(grun.twfe, effect = "time")

## -----------------------------------------------------------------------------
grun.amem <- plm(inv~value+capital, data=Grunfeld,
                 model="random", random.method="amemiya")

## -----------------------------------------------------------------------------
ercomp(inv~value+capital, data=Grunfeld, method = "amemiya", effect = "twoways")

## ----2RE-amemiya--------------------------------------------------------------
grun.tways <- plm(inv~value+capital, data = Grunfeld, effect = "twoways",
                  model = "random", random.method = "amemiya")
summary(grun.tways)

## ----hedonic------------------------------------------------------------------
data("Hedonic", package = "plm")
Hed <- plm(mv~crim+zn+indus+chas+nox+rm+age+dis+rad+tax+ptratio+blacks+lstat,
           data = Hedonic, model = "random", index = "townid")
summary(Hed)

## ----hedonic-punbal-----------------------------------------------------------
punbalancedness(Hed)

## ----G2SLS--------------------------------------------------------------------
data("Crime", package = "plm")
cr <- plm(lcrmrte ~ lprbarr + lpolpc + lprbconv + lprbpris + lavgsen +
          ldensity + lwcon + lwtuc + lwtrd + lwfir + lwser + lwmfg + lwfed +
          lwsta + lwloc + lpctymle + lpctmin + region + smsa + factor(year)
          | . - lprbarr - lpolpc + ltaxpc + lmix,
          data = Crime, model = "random")
summary(cr)

## ----hausman-taylor-----------------------------------------------------------
ht <- plm(lwage ~ wks + south + smsa + married + exp + I(exp ^ 2) + 
              bluecol + ind + union + sex + black + ed |
              bluecol + south + smsa + ind + sex + black |
              wks + married + union + exp + I(exp ^ 2), 
          data = Wages, index = 595,
          model = "random", random.method = "ht", inst.method = "baltagi")
summary(ht)

## ----grunfeld.within----------------------------------------------------------
grun.varw <- pvcm(inv~value+capital, data=Grunfeld, model="within")
grun.varr <- pvcm(inv~value+capital, data=Grunfeld, model="random")
summary(grun.varr)

## ----gmm----------------------------------------------------------------------
emp.gmm <- pgmm(log(emp)~lag(log(emp), 1:2)+lag(log(wage), 0:1)+log(capital)+
                lag(log(output), 0:1) | lag(log(emp), 2:99),
                data = EmplUK, effect = "twoways", model = "twosteps")
summary(emp.gmm)

## ----gmm2---------------------------------------------------------------------
z2 <- pgmm(log(emp) ~ lag(log(emp), 1)+ lag(log(wage), 0:1) +
           lag(log(capital), 0:1) | lag(log(emp), 2:99) +
           lag(log(wage), 2:99) + lag(log(capital), 2:99),        
           data = EmplUK, effect = "twoways", model = "onestep", 
           transformation = "ld")
summary(z2, robust = TRUE)

## ----pggls--------------------------------------------------------------------
zz <- pggls(log(emp)~log(wage)+log(capital), data=EmplUK, model="pooling")
summary(zz)

## -----------------------------------------------------------------------------
zz <- pggls(log(emp)~log(wage)+log(capital), data=EmplUK, model="within")

## -----------------------------------------------------------------------------
znp <- pvcm(inv~value+capital, data=Grunfeld, model="within")
zplm <- plm(inv~value+capital, data=Grunfeld, model="within")
pooltest(zplm, znp)

## ----results='hide'-----------------------------------------------------------
pooltest(inv~value+capital, data=Grunfeld, model="within")

## -----------------------------------------------------------------------------
g <- plm(inv ~ value + capital, data=Grunfeld, model="pooling")
plmtest(g, effect="twoways", type="ghm")

## ----results='hide'-----------------------------------------------------------
plmtest(inv~value+capital, data=Grunfeld, effect="twoways", type="ghm")

## -----------------------------------------------------------------------------
gw <- plm(inv ~ value + capital, data=Grunfeld, effect="twoways", model="within")
gp <- plm(inv ~ value + capital, data=Grunfeld, model="pooling")
pFtest(gw, gp)

## ----results='hide'-----------------------------------------------------------
pFtest(inv~value+capital, data=Grunfeld, effect="twoways")

## -----------------------------------------------------------------------------
gw <- plm(inv~value+capital, data=Grunfeld, model="within")
gr <- plm(inv~value+capital, data=Grunfeld, model="random")
phtest(gw, gr)

## ----wtest--------------------------------------------------------------------
pwtest(log(gsp)~log(pcap)+log(pc)+log(emp)+unemp, data=Produc)

## ----pbsytestJoint------------------------------------------------------------
pbsytest(log(gsp)~log(pcap)+log(pc)+log(emp)+unemp, data=Produc, test="j")

## ----pbsytestAR---------------------------------------------------------------
pbsytest(log(gsp)~log(pcap)+log(pc)+log(emp)+unemp, data=Produc)

## ----pbsytestRE---------------------------------------------------------------
pbsytest(log(gsp)~log(pcap)+log(pc)+log(emp)+unemp, data=Produc, test="re")

## ----pbltest------------------------------------------------------------------
pbltest(log(gsp)~log(pcap)+log(pc)+log(emp)+unemp, 
        data=Produc, alternative="onesided")

## ----generalAR----------------------------------------------------------------
pbgtest(grun.fe, order = 2)

## ----pwartest-----------------------------------------------------------------
pwartest(log(emp) ~ log(wage) + log(capital), data=EmplUK)

## ----pwfdtest1----------------------------------------------------------------
pwfdtest(log(emp) ~ log(wage) + log(capital), data=EmplUK)

## ----pwfdtest2----------------------------------------------------------------
pwfdtest(log(emp) ~ log(wage) + log(capital), data=EmplUK, h0="fe")

## ----pcdtest1-----------------------------------------------------------------
pcdtest(inv~value+capital, data=Grunfeld)

## ----pcdtest2-----------------------------------------------------------------
pcdtest(inv~value+capital, data=Grunfeld, model="within")

## ----levinlin-----------------------------------------------------------------
data("HousePricesUS", package = "pder")
lprice <- log(pdata.frame(HousePricesUS)$price)
(lev <- purtest(lprice, test = "levinlin", lags = 2, exo = "trend"))
summary(lev) ### gives details

## ----ips----------------------------------------------------------------------
purtest(lprice, test = "ips", lags = 2, exo = "trend")

## ----phansi1------------------------------------------------------------------
### input is numeric (p-values), replicates Hanck (2013), Table 11 (left side)
pvals <- c(0.0001,0.0001,0.0001,0.0001,0.0001,0.0001,0.0050,0.0050,0.0050,
           0.0050,0.0175,0.0175,0.0200,0.0250,0.0400,0.0500,0.0575,0.2375,0.2475)
countries <- c("Argentina","Sweden","Norway","Mexico","Italy","Finland","France",
              "Germany","Belgium","U.K.","Brazil","Australia","Netherlands",
              "Portugal","Canada", "Spain","Denmark","Switzerland","Japan")
names(pvals) <- countries
h <- phansi(pvals)
print(h)
h$rejected # logical indicating the individuals with rejected individual H0

## ----phansi2, results='hide'--------------------------------------------------
### input is a (suitable) purtest object / different example
y <- data.frame(split(Grunfeld$inv, Grunfeld$firm))
obj <- purtest(y, pmax = 4, exo = "intercept", test = "madwu")
phansi(obj, alpha = 0.06) # test with significance level set to 6 %

## ----vcovHC1------------------------------------------------------------------
re <- plm(inv~value+capital, data = Grunfeld, model = "random")
summary(re, vcov = vcovHC) # gives usual summary output but with robust test statistics

library("lmtest")
coeftest(re, vcovHC, df = Inf)

## ----vcovHC2, results='hide'--------------------------------------------------
summary(re, vcov = vcovHC(re, method="white2", type="HC3"), df = Inf)
coeftest(re, vcovHC(re, method="white2", type="HC3"), df = Inf)

## ----waldtest-vcovHC----------------------------------------------------------
waldtest(re, update(re, . ~ . -capital),
         vcov=function(x) vcovHC(x, method="white2", type="HC3"))

## ----car-vcovHC---------------------------------------------------------------
library("car")
linearHypothesis(re, "2*value=capital", vcov. = vcovHC)

## ----re2----------------------------------------------------------------------
library(nlme)
reGLS <- plm(inv~value+capital, data=Grunfeld, model="random")

reML <- lme(inv~value+capital, data=Grunfeld, random=~1|firm)

coef(reGLS)

summary(reML)$coefficients$fixed

## ----vcmrand------------------------------------------------------------------
vcm <- pvcm(inv~value+capital, data=Grunfeld, model="random", effect="time")

vcmML <- lme(inv~value+capital, data=Grunfeld, random=~value+capital|year)

coef(vcm)

summary(vcmML)$coefficients$fixed

## ----vcmfixed-----------------------------------------------------------------
vcmf <- pvcm(inv~value+capital, data=Grunfeld, model="within", effect="time")

vcmfML <- lmList(inv~value+capital|year, data=Grunfeld)

## ----gglsre-------------------------------------------------------------------
sGrunfeld <- Grunfeld[Grunfeld$firm %in% 4:6, ]

ggls <- pggls(inv~value+capital, data=sGrunfeld, model="pooling")

gglsML <- gls(inv~value+capital, data=sGrunfeld,
              correlation=corSymm(form=~1|year))

coef(ggls)

summary(gglsML)$coefficients

## ----lmAR1--------------------------------------------------------------------
Grunfeld$year <- as.numeric(as.character(Grunfeld$year))
lmAR1ML <- gls(inv~value+capital,data=Grunfeld,
               correlation=corAR1(0,form=~year|firm))

## ----reAR1--------------------------------------------------------------------
reAR1ML <- lme(inv~value+capital, data=Grunfeld,random=~1|firm,
               correlation=corAR1(0,form=~year|firm))

## ----fetchcoefs---------------------------------------------------------------
summary(reAR1ML)$coefficients$fixed
coef(reAR1ML$modelStruct$corStruct, unconstrained=FALSE)

## ----LRar---------------------------------------------------------------------
lmML <- gls(inv~value+capital, data=Grunfeld)
anova(lmML, lmAR1ML)

## ----LRarsubRE----------------------------------------------------------------
anova(reML, reAR1ML)

## ----LRre---------------------------------------------------------------------
anova(lmML, reML)

## ----LRresubAR----------------------------------------------------------------
anova(lmAR1ML, reAR1ML)

