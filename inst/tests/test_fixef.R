### Test of fixef
###
### (1): general tests
### (2): consistency with summary.plm
###
### see also:
###    * test file test_within_intercept.R for consistency checks
###      between functions fixef and within_intercept
###    * test file test_fixef_comp_lm_plm.R for a comparison of the fixed effects to LSDV models via lm()

############# (1): general run tests #############
library(plm)
data("Grunfeld", package = "plm")

# balanced models
gi  <- plm(inv ~ value + capital, data = Grunfeld, model = "within")
gt  <- plm(inv ~ value + capital, data = Grunfeld, model = "within", effect = "time")
gtw <- plm(inv ~ value + capital, data = Grunfeld, model = "within", effect = "twoways")


f_level             <- fixef(gi, type = "level")
f_level_robust_mat  <- fixef(gi, type = "level", vcov = vcovHC(gi)) # vcov is matrix
f_level_robust_func <- fixef(gi, type = "level", vcov = vcovHC)     # vcov is function

print(attr(f_level,             "se"))
print(attr(f_level_robust_func, "se"))

print(summary(f_level),             digits = 8)
print(summary(f_level_robust_func), digits = 8)

f_level_t             <- fixef(gt, type = "level")
f_level_t_robust_func <- fixef(gt, type = "level", vcov = vcovHC)     # vcov is function

print(attr(f_level_t,             "se"))
print(attr(f_level_t_robust_func, "se"))

print(summary(f_level_t),             digits = 8)
print(summary(f_level_t_robust_func), digits = 8)

f_level_d             <- fixef(gtw, type = "level")
f_level_d_robust_func <- fixef(gtw, type = "level", vcov = vcovHC)     # vcov is function

print(attr(f_level_d,             "se"))
print(attr(f_level_d_robust_func, "se"))

print(summary(f_level_d),             digits = 8)
print(summary(f_level_d_robust_func), digits = 8)

# just run tests for type = "dmean" and type = "dfirst"
fixef(gi, type = "dmean")
fixef(gt, type = "dmean")
fixef(gtw, effect = "individual", type = "dmean")
fixef(gtw, effect = "time",       type = "dmean")
fixef(gtw, effect = "twoways",    type = "dmean")


fixef(gi, type = "dfirst")
fixef(gt, type = "dfirst")
fixef(gtw, effect = "individual", type = "dfirst")
fixef(gtw, effect = "time",       type = "dfirst")
fixef(gtw, effect = "twoways",    type = "dfirst")

fixef(gtw, effect = "twoways",    type = "level")


# this errored until 2022-04-22:
mod2 <- plm(inv ~ value + capital, data = Grunfeld[Grunfeld$firm %in% c(1,2), ])
fixef(mod2, type = "dfirst")

# until 2022-04-22: this gave non-sane result c(NA, 0), now numeric(0)
mod1 <- plm(inv ~ value + capital, data = Grunfeld[Grunfeld$firm %in% c(1), ])
fixef(mod1, type = "dfirst")


############# (2): consistency with summary.plm #############
# compare summary.plm to summary.fixef( , type = "dfirst")
mod_pool <- plm(inv ~ value + capital + factor(firm), data = Grunfeld, model = "pooling")
sum_mod_pool <- summary(mod_pool)
f_dfirst <- fixef(gi, type = "dfirst")
sum_f_dfirst <- summary(f_dfirst)

if(!isTRUE(all.equal(sum_mod_pool[["coefficients"]][-c(1:3) , "Estimate"], sum_f_dfirst[ , "Estimate"], check.attributes = FALSE)))
  stop("estimates diverge: summary.plm vs. summary.fixef(..., type = \"dfirst\")")

if(!isTRUE(all.equal(sum_mod_pool[["coefficients"]][-c(1:3) , "Std. Error"], sum_f_dfirst[ , "Std. Error"], check.attributes = FALSE)))
  stop("standard errors diverge: summary.plm vs. summary.fixef(..., type = \"dfirst\")")

if(!isTRUE(all.equal(sum_mod_pool[["coefficients"]][-c(1:3) , "t-value"], sum_f_dfirst[ , "t-value"], check.attributes = FALSE)))
  stop("t-values diverge: summary.plm vs. summary.fixef(..., type = \"dfirst\")")

if(!isTRUE(all.equal(sum_mod_pool[["coefficients"]][-c(1:3) , "Pr(>|t|)"], sum_f_dfirst[ , "Pr(>|t|)"], check.attributes = FALSE)))
  stop("p-values diverge: summary.plm vs. summary.fixef(..., type = \"dfirst\")")


###### compare to package lfe:
## Standard errors are bootstrapped in lfe
##  -> different SE results compared to plm
##  -> different SE results for every new call
#
# library(lfe)
# data("Grunfeld", package = "plm")
# mod_felm <- felm(inv ~ value + capital | firm, data = Grunfeld)
# summary(mod_felm)
# 
# fe_lfe <- getfe(mod_felm, se = TRUE, bN = 50)
# print(fe_lfe)

# sum_f_level <- summary(f_level)
# print(sum_f_level)
