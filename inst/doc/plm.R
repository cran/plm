###################################################
### chunk number 1: 
###################################################
options(prompt= "R> ", useFancyQuotes = FALSE)


###################################################
### chunk number 2: 
###################################################
library("plm")
#library(nlme)
#library(MASS)
#library(kinship)
#source("~/R/plm/chargement.R")


###################################################
### chunk number 3: 
###################################################
data("EmplUK", package="plm")
data("Produc", package="plm")
data("Grunfeld", package="plm")
data("Wages",package="plm")



###################################################
### chunk number 4: setdata1
###################################################
Wages <- plm.data(Wages, index = 595)



###################################################
### chunk number 5: 
###################################################
log(emp)~lag(log(emp),1)+lag(log(emp),2)+lag(log(wage),2)+lag(log(wage),3)+diff(capital,2)+diff(capital,3)


###################################################
### chunk number 6: 
###################################################
dynformula(emp~wage+capital,log=list(capital=FALSE,TRUE),lag=list(emp=2,c(2,3)),diff=list(FALSE,capital=TRUE))


###################################################
### chunk number 7: 
###################################################
grun.fe <- plm(inv~value+capital,data=Grunfeld,model="within")
grun.re <- plm(inv~value+capital,data=Grunfeld,model="random")


###################################################
### chunk number 8: 
###################################################
summary(grun.re)


###################################################
### chunk number 9: 
###################################################
fixef(grun.fe)


###################################################
### chunk number 10: 
###################################################
summary(fixef(grun.fe))


###################################################
### chunk number 11: 
###################################################
grun.twfe <- plm(inv~value+capital,data=Grunfeld,model="within",effect="twoways")
fixef(grun.twfe,effect="time")


###################################################
### chunk number 12: 
###################################################
grun.amem <- plm(inv~value+capital,data=Grunfeld,model="random",random.method="amemiya")



###################################################
### chunk number 13: 
###################################################
grun.tways <- plm(inv~value+capital,data=Grunfeld,effect="twoways",model="random",random.method="amemiya")
summary(grun.tways)


###################################################
### chunk number 14: pbron
###################################################
emp.iv <- plm(emp~wage+capital|lag(wage,1)+capital,data=EmplUK,model="random")
emp.iv <- plm(emp~wage+capital|.-wage+lag(wage,1),data=EmplUK,model="random")


###################################################
### chunk number 15: 
###################################################
data("Crime", package = "plm")
cr <- plm(log(crmrte) ~ log(prbarr) + log(polpc) + log(prbconv) +
         log(prbpris) + log(avgsen) + log(density) + log(wcon) + 
         log(wtuc) + log(wtrd) + log(wfir) + log(wser) + log(wmfg) + 
         log(wfed) + log(wsta) + log(wloc) + log(pctymle) + log(pctmin) + 
         region + smsa + factor(year) | . - log(prbarr) -log(polpc) + 
         log(taxpc) + log(mix), data = Crime,
         model = "random")


###################################################
### chunk number 16: 
###################################################
ht <- plm(lwage~wks+south+smsa+married+exp+I(exp^2)+
          bluecol+ind+union+sex+black+ed | 
          sex+black+bluecol+south+smsa+ind,
          data=Wages,model="ht",index=595)
summary(ht)


###################################################
### chunk number 17: grunfeld.within
###################################################
grun.varw <- pvcm(inv~value+capital,data=Grunfeld,model="within")
grun.varr <- pvcm(inv~value+capital,data=Grunfeld,model="random")
summary(grun.varr)


###################################################
### chunk number 18: gmm
###################################################
emp.gmm <- pgmm(dynformula(log(emp)~log(wage)+log(capital)+log(output),lag=list(2,1,0,1)),EmplUK,effect="twoways",model="twosteps",gmm.inst=~log(emp),lag.gmm=list(c(2,99)))
summary(emp.gmm)


###################################################
### chunk number 19: pggls
###################################################
zz <- pggls(log(emp)~log(wage)+log(capital),data=EmplUK,model="random")
summary(zz)


###################################################
### chunk number 20: 
###################################################
zz <- pggls(log(emp)~log(wage)+log(capital),data=EmplUK,model="within")


###################################################
### chunk number 21: 
###################################################
znp <- pvcm(inv~value+capital,data=Grunfeld,model="within")
zplm <- plm(inv~value+capital,data=Grunfeld)
pooltest(zplm,znp)


###################################################
### chunk number 22: 
###################################################
pooltest(inv~value+capital,data=Grunfeld,model="within")


###################################################
### chunk number 23: 
###################################################
g <- plm(inv ~ value + capital,data=Grunfeld,model="pooling")
plmtest(g,effect="twoways",type="ghm")


###################################################
### chunk number 24: 
###################################################
plmtest(inv~value+capital,data=Grunfeld,effect="twoways",type="ghm")


###################################################
### chunk number 25: 
###################################################
gw <- plm(inv ~ value + capital,data=Grunfeld,effect="twoways",model="within")
gp <- plm(inv ~ value + capital,data=Grunfeld,model="pooling")
pFtest(gw,gp)


###################################################
### chunk number 26: 
###################################################
pFtest(inv~value+capital,data=Grunfeld,effect="twoways")


###################################################
### chunk number 27: 
###################################################
gw <- plm(inv~value+capital,data=Grunfeld,model="within")
gr <- plm(inv~value+capital,data=Grunfeld,model="random")
phtest(gw, gr)


###################################################
### chunk number 28: wtest
###################################################
pwtest(log(gsp)~log(pcap)+log(pc)+log(emp)+unemp, data=Produc)


###################################################
### chunk number 29: pbsytestJoint
###################################################
pbsytest(log(gsp)~log(pcap)+log(pc)+log(emp)+unemp,data=Produc,test="j")


###################################################
### chunk number 30: pbsytestAR
###################################################
pbsytest(log(gsp)~log(pcap)+log(pc)+log(emp)+unemp,data=Produc)


###################################################
### chunk number 31: pbsytestRE
###################################################
pbsytest(log(gsp)~log(pcap)+log(pc)+log(emp)+unemp,data=Produc,test="re")


###################################################
### chunk number 32: pbltest
###################################################
pbltest(log(gsp)~log(pcap)+log(pc)+log(emp)+unemp,data=Produc,alternative="onesided")


###################################################
### chunk number 33: generalAR
###################################################
## this can be taken away as soon as attached to plm.rnw
grun.fe <- plm(inv ~ value + capital, data = Grunfeld, model = "within")
pbgtest(grun.fe, order=2)


###################################################
### chunk number 34: pwartest
###################################################
pwartest(log(emp) ~ log(wage) + log(capital), data=EmplUK)


###################################################
### chunk number 35: pwfdtest1
###################################################
pwfdtest(log(emp) ~ log(wage) + log(capital), data=EmplUK)


###################################################
### chunk number 36: pwfdtest2
###################################################
pwfdtest(log(emp) ~ log(wage) + log(capital), data=EmplUK, h0="fe")


###################################################
### chunk number 37: pcdtest1
###################################################
pcdtest(inv~value+capital, data=Grunfeld)


###################################################
### chunk number 38: pcdtest2
###################################################
pcdtest(inv~value+capital, data=Grunfeld, model="within")


###################################################
### chunk number 39: 
###################################################
library("lmtest")
re <- plm(inv~value+capital,data=Grunfeld,model="random")
coeftest(re,vcovHC)


###################################################
### chunk number 40: 
###################################################
coeftest(re,vcovHC(re,method="white2",type="HC3"))


###################################################
### chunk number 41: 
###################################################
waldtest(re,update(re,.~.-capital),vcov=function(x) vcovHC(x,method="white2",type="HC3"))


###################################################
### chunk number 42: 
###################################################
library("car")
linear.hypothesis(re, "2*value=capital", vcov.=vcovHC)


###################################################
### chunk number 43: re
###################################################

reGLS<-plm(inv~value+capital,data=Grunfeld,model="random")

reML<-lme(inv~value+capital,data=Grunfeld,random=~1|firm)

coef(reGLS)

summary(reML)$coef$fixed



###################################################
### chunk number 44: vcmrand
###################################################
vcm<-pvcm(inv~value+capital,data=Grunfeld,model="random",effect="time")

vcmML<-lme(inv~value+capital,data=Grunfeld,random=~value+capital|year)

coef(vcm)

summary(vcmML)$coef$fixed



###################################################
### chunk number 45: vcmfixed
###################################################
vcmf<-pvcm(inv~value+capital,data=Grunfeld,model="within",effect="time")

vcmfML<-lmList(inv~value+capital|year,data=Grunfeld)



###################################################
### chunk number 46: gglsre
###################################################
sGrunfeld <- Grunfeld[Grunfeld$firm%in%4:6,]

ggls<-pggls(inv~value+capital,data=sGrunfeld,model="random")

gglsML<-gls(inv~value+capital,data=sGrunfeld,
            correlation=corSymm(form=~1|year))

coef(ggls)

summary(gglsML)$coef


###################################################
### chunk number 47: lmAR1
###################################################
Grunfeld$year <- as.numeric(as.character(Grunfeld$year))
lmAR1ML<-gls(inv~value+capital,data=Grunfeld,
             correlation=corAR1(0,form=~year|firm))


###################################################
### chunk number 48: reAR1
###################################################
reAR1ML<-lme(inv~value+capital,data=Grunfeld,random=~1|firm,
             correlation=corAR1(0,form=~year|firm))


###################################################
### chunk number 49: fetchcoefs
###################################################
summary(reAR1ML)$coef$fixed
coef(reAR1ML$modelStruct$corStruct,unconstrained=FALSE)


###################################################
### chunk number 50: LRar
###################################################
lmML<-gls(inv~value+capital,data=Grunfeld)
anova(lmML,lmAR1ML)


###################################################
### chunk number 51: LRarsubRE
###################################################
anova(reML,reAR1ML)


###################################################
### chunk number 52: LRre
###################################################
anova(lmML,reML)


###################################################
### chunk number 53: LRresubAR
###################################################
anova(lmAR1ML,reAR1ML)


