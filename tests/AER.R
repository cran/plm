### systemfit commente pour l'instant


#### Fatalities
library("AER")
## data from Stock and Watson (2007)
data("Fatalities")
## add fatality rate (number of traffic deaths
## per 10,000 people living in that state in that year)
Fatalities$frate <- with(Fatalities, fatal/pop * 10000)
## add discretized version of minimum legal drinking age
Fatalities$drinkagec <- cut(Fatalities$drinkage,
  breaks = 18:22, include.lowest = TRUE, right = FALSE)
Fatalities$drinkagec <- relevel(Fatalities$drinkagec, ref = 4)
## any punishment?
Fatalities$punish <- with(Fatalities,
  factor(jail == "yes" | service == "yes", labels = c("no", "yes")))
## plm package
library("plm")

## for comparability with Stata we use HC1 below
## p. 351, Eq. (10.2)
f1982 <- subset(Fatalities, year == "1982")
fm_1982 <- lm(frate ~ beertax, data = f1982)
coeftest(fm_1982, vcov = vcovHC(fm_1982, type = "HC1"))

## p. 353, Eq. (10.3)
f1988 <- subset(Fatalities, year == "1988")
fm_1988 <- lm(frate ~ beertax, data = f1988)
coeftest(fm_1988, vcov = vcovHC(fm_1988, type = "HC1"))

## pp. 355, Eq. (10.8)
fm_diff <- lm(I(f1988$frate - f1982$frate) ~ I(f1988$beertax - f1982$beertax))
coeftest(fm_diff, vcov = vcovHC(fm_diff, type = "HC1"))

## pp. 360, Eq. (10.15)
##   (1) via formula
fm_sfe <- lm(frate ~ beertax + state - 1, data = Fatalities)
##   (2) by hand
fat <- with(Fatalities,
  data.frame(frates = frate - ave(frate, state),
  beertaxs = beertax - ave(beertax, state)))
fm_sfe2 <- lm(frates ~ beertaxs - 1, data = fat)
##   (3) via plm()
fm_sfe3 <- plm(frate ~ beertax, data = Fatalities,
  index = c("state", "year"), model = "within")

coeftest(fm_sfe, vcov = vcovHC(fm_sfe, type = "HC1"))[1,]
coeftest(fm_sfe2, vcov = vcovHC(fm_sfe2, type = "HC1"))[1,] ## uses different df in sd and pval
coeftest(fm_sfe3, vcov = vcovHC(fm_sfe3, type = "HC1", method = "white1"))[,1] ## uses different df in pval


## pp. 363, Eq. (10.21)
## via lm()
fm_stfe <- lm(frate ~ beertax + state + year - 1, data = Fatalities)
coeftest(fm_stfe, vcov = vcovHC(fm_stfe, type = "HC1"))[1,]
## via plm()
fm_stfe2 <- plm(frate ~ beertax, data = Fatalities,
  index = c("state", "year"), model = "within", effect = "twoways")
coeftest(fm_stfe2, vcov = vcovHC) ## different


## p. 368, Table 10.1, numbers refer to cols.
fm1 <- plm(frate ~ beertax, data = Fatalities, index = c("state", "year"), model = "pooling")
fm2 <- plm(frate ~ beertax, data = Fatalities, index = c("state", "year"), model = "within")
fm3 <- plm(frate ~ beertax, data = Fatalities, index = c("state", "year"), model = "within",
  effect = "twoways")
fm4 <- plm(frate ~ beertax + drinkagec + jail + service + miles + unemp + log(income),
  data = Fatalities, index = c("state", "year"), model = "within", effect = "twoways")
fm5 <- plm(frate ~ beertax + drinkagec + jail + service + miles,
  data = Fatalities, index = c("state", "year"), model = "within", effect = "twoways")
fm6 <- plm(frate ~ beertax + drinkage + punish + miles + unemp + log(income),
  data = Fatalities, index = c("state", "year"), model = "within", effect = "twoways")
fm7 <- plm(frate ~ beertax + drinkagec + jail + service + miles + unemp + log(income),
  data = Fatalities, index = c("state", "year"), model = "within", effect = "twoways")
## summaries not too close, s.e.s generally too small
coeftest(fm1, vcov = vcovHC)
coeftest(fm2, vcov = vcovHC)
coeftest(fm3, vcov = vcovHC)
coeftest(fm4, vcov = vcovHC)
coeftest(fm5, vcov = vcovHC)
coeftest(fm6, vcov = vcovHC)
coeftest(fm7, vcov = vcovHC)


## Grunfeld

## Panel models
library("plm")
data("Grunfeld", package = "AER")

pg <- plm.data(subset(Grunfeld, firm != "American Steel"), c("firm", "year"))

fm_fe <- plm(invest ~ value + capital, model = "within", data = pg)
summary(fm_fe)
coeftest(fm_fe, vcov = vcovHC)

fm_reswar <- plm(invest ~ value + capital, data = pg,
  model = "random", random.method = "swar")
summary(fm_reswar)

## testing for random effects
fm_ols <- plm(invest ~ value + capital, data = pg, model = "pooling")
plmtest(fm_ols, type = "bp")
plmtest(fm_ols, type = "honda")

## Random effects models
fm_ream <- plm(invest ~ value + capital, data = pg, model = "random", random.method = "amemiya")
fm_rewh <- plm(invest ~ value + capital, data = pg, model = "random", random.method = "walhus")
fm_rener <- plm(invest ~ value + capital, data = pg, model = "random", random.method = "nerlove")

## Baltagi (2005), Tab. 2.1
rbind(
  "OLS(pooled)" = coef(fm_ols),
  "FE" = c(NA, coef(fm_fe)),
  "RE-SwAr" = coef(fm_reswar),
  "RE-Amemiya" = coef(fm_ream),
  "RE-WalHus" = coef(fm_rewh),
  "RE-Nerlove" = coef(fm_rener))

## Hausman test
phtest(fm_fe, fm_reswar)

library("plm")
ggr <- subset(Grunfeld, firm %in% c("General Motors", "US Steel",
  "General Electric", "Chrysler", "Westinghouse"))
ggr[c(26, 38), 1] <- c(261.6, 645.2)
ggr[32, 3] <- 232.6

pggr <- plm.data(ggr, c("firm", "year"))
## library("systemfit")
## fm_ols <- systemfit(invest ~ value + capital, data = pggr, method = "OLS")
## fm_pols <- systemfit(invest ~ value + capital, data = pggr, method = "OLS", pooled = TRUE)

## PSID7682

data("PSID7682")

library("plm")
psid <- plm.data(PSID7682, c("id", "year"))

## Baltagi & Khanti-Akom, Table I, column "HT"
## original Cornwell & Rupert choice of exogenous variables
psid_ht1 <- plm(log(wage) ~ weeks + south + smsa + married +
  experience + I(experience^2) + occupation + industry + union + gender + ethnicity + education |
  weeks + south + smsa + married + gender + ethnicity,
  data = psid, model = "ht")

## Baltagi & Khanti-Akom, Table II, column "HT"
## alternative choice of exogenous variables
psid_ht2 <- plm(log(wage) ~ occupation + south + smsa + industry +
  experience + I(experience^2) + weeks + married + union + gender + ethnicity + education |
  occupation + south + smsa + industry + gender + ethnicity,
  data = psid, model = "ht")

## Baltagi & Khanti-Akom, Table III, column "HT"
## original choice of exogenous variables + time dummies
## (see also Baltagi, 2001, Table 7.1)
psid$time <- psid$year
psid_ht3 <- plm(log(wage) ~ weeks + south + smsa + married +
  experience + I(experience^2) + occupation + industry + union + gender + ethnicity + education + time |
  weeks + south + smsa + married + gender + ethnicity + time,
  data = psid, model = "ht")


## USAirlines

## alternatively, use plm()

data("USAirlines")

library("plm")
usair <- plm.data(USAirlines, c("firm", "year"))
fm_full2 <- plm(log(cost) ~ log(output) + I(log(output)^2) + log(price) + load,
  data = usair, model = "within", effect = "twoways")
fm_time2 <- plm(log(cost) ~ log(output) + I(log(output)^2) + log(price) + load,
  data = usair, model = "within", effect = "time")
fm_firm2 <- plm(log(cost) ~ log(output) + I(log(output)^2) + log(price) + load,
  data = usair, model = "within", effect = "individual")
fm_no2 <- plm(log(cost) ~ log(output) + I(log(output)^2) + log(price) + load,
  data = usair, model = "pooling")
pFtest(fm_full2, fm_time2)
pFtest(fm_full2, fm_firm2)
pFtest(fm_full2, fm_no2)
