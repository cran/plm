## ----setup, echo=FALSE--------------------------------------------------------
library("knitr")
opts_chunk$set(message = FALSE, warning = FALSE)

## -----------------------------------------------------------------------------
library("plm")
data("SeatBelt", package = "pder")
SeatBelt$occfat <- with(SeatBelt, log(farsocc / (vmtrural + vmturban)))
pSB <- pdata.frame(SeatBelt)

## -----------------------------------------------------------------------------
formols <- occfat ~ log(usage) + log(percapin)
mfols <- model.frame(pSB, formols)
Xols <- model.matrix(mfols)
y <- pmodel.response(mfols)
coef(lm.fit(Xols, y))

## -----------------------------------------------------------------------------
coef(plm(formols, SeatBelt, model = "pooling"))

## -----------------------------------------------------------------------------
formiv1 <- occfat ~ log(usage) + log(percapin) | log(percapin) + ds + dp + dsp
formiv2 <- occfat ~ log(usage) + log(percapin) | . - log(usage) + ds + dp + dsp

## -----------------------------------------------------------------------------
mfSB1 <- model.frame(pSB, formiv1)
X1 <- model.matrix(mfSB1, rhs = 1)
W1 <- model.matrix(mfSB1, rhs = 2)
head(X1, 3) ; head(W1, 3)

## -----------------------------------------------------------------------------
library("Formula")
head(model.frame(Formula(formiv2), SeatBelt), 3)
head(model.frame(Formula(formiv2), SeatBelt, dot = "previous"), 3)

## -----------------------------------------------------------------------------
mfSB2 <- model.frame(pSB, formiv2)
X2 <- model.matrix(mfSB2, rhs = 1)
W2 <- model.matrix(mfSB2, rhs = 2)
head(X2, 3) ; head(W2, 3)

## -----------------------------------------------------------------------------
HX1 <- lm.fit(W1, X1)$fitted.values
head(HX1, 3)

## -----------------------------------------------------------------------------
coef(lm.fit(HX1, y))

## -----------------------------------------------------------------------------
coef(plm(formiv1, SeatBelt, model = "pooling"))

## -----------------------------------------------------------------------------
coef(AER::ivreg(formiv1, data = SeatBelt))

## ----eval = FALSE, include = FALSE--------------------------------------------
# X2 <- model.matrix(Formula(form1), mfSB, rhs = 2, dot = "previous")
# 
# formols <- occfat ~ log(usage) + log(percapin)  | . - log(usage) +  ds + dp + dsp
# 
# form1 <- occfat ~ log(usage) + log(percapin) + log(unemp) + log(meanage) +
#     log(precentb) + log(precenth) + log(densrur) + log(densurb) +
#     log(viopcap) + log(proppcap) + log(vmtrural) + log(vmturban) +
#     log(fueltax) + lim65 + lim70p + mlda21 + bac08
# form2 <- . ~ . |  . - log(usage) + ds + dp +dsp
# 
# jorm1 <- occfat ~ log(usage) + log(percapin) + log(unemp) + log(meanage) +
#     log(precentb) + log(precenth) + log(densrur) + log(densurb) +
#     log(viopcap) + log(proppcap) + log(vmtrural) + log(vmturban) +
#     log(fueltax) + lim65 + lim70p + mlda21 + bac08 | . - log(usage) +
#     ds + dp + dsp
# jorm2 <- noccfat ~ . | .

