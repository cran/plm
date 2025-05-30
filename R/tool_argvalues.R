## This file contain named vectors of the acceptable values for different
## arguments used in plm functions.


random.method.list <- c(swar    = "Swamy-Arora",
                        walhus  = "Wallace-Hussain",
                        amemiya = "Amemiya",
                        nerlove = "Nerlove",
                        ht      = "Hausman-Taylor")

effect.plm.list <- c(individual = "Oneway (individual) effect",
                     time       = "Oneway (time) effect",
                     twoways    = "Twoways effects",
                     nested     = "Nested effects")

effect.pvcm.list <- c(individual  = "Oneway (individual) effect",
                      time        = "Oneway (time) effect")

effect.pggls.list <- c(individual = "Oneway (individual) effect",
                       time       = "Oneway (time) effect")

effect.pgmm.list <- c(individual = "Oneway (individual) effect",
                      twoways    = "Twoways effects")

model.plm.list <- c(pooling = "Pooling",
                    within  = "Within",
                    between = "Between",
                    random  = "Random Effect",
                    ht      = "Hausman-Taylor",
                    fd      = "First-Difference")

ht.method.list <- c(ht = "Hausman-Taylor estimator",
                    am = "Amemiya-MaCurdy estimator",
                    bms = "Breusch-Mizon-Schmidt estimator")

model.pvcm.list <- c(within = "No-pooling model",
                     random = "Random coefficients model")

model.pggls.list <- c(within  = "Within FGLS model",
                      random  = "General FGLS model",
                      pooling = "General FGLS model",
                      fd      = "First-Difference FGLS model")

model.pgmm.list <- c(onestep  = "One-step model",
                     twosteps = "Two-steps model")

model.pgmm.transformation.list <- c(d  = "Difference GMM",
                                    ld = "System GMM")

model.pcce.list <- c(ccemg = "Mean Groups model",
                     ccep  = "Pooled model")

model.pmg.list <- c(mg  = "Mean Groups model",
                    dmg = "Demeaned Mean Groups model",
                    cmg = "Common Correlated Effects Mean Groups model")

inst.method.list <- c(bvk     = "Balestra-Varadharajan-Krishnakumar",
                      baltagi = "Baltagi",
                      am      = "Amemiya-MaCurdy",
                      bms     = "Breusch-Mizon-Schmidt")

# robust.list and weights.list are not used anywhere...
robust.list <- c(white1   = "White 1",
                 white2   = "White 2",
                 arellano = "Arellano")

weights.list <- c(HC0 = "HC0",
                  HC1 = "HC1",
                  HC2 = "HC2",
                  HC3 = "HC3",
                  HC4 = "HC4")

oneof <- function(x){
  x <- names(x)
  last <- x[length(x)]
  x <- x[-length(x)]
  x <- paste(x,collapse=", ")
  x <- paste(x,last,sep=" and ")
  x
}

plm.arg <- c("formula", "data", "subset", "weights", "na.action", "effect", "model",
             "instruments", "random.method", "inst.method", "index")

pgmm.fsm.list <- c(I    = "I",
                   G    = "G",
                   GI   = "GI",
                   full = "full")
                              