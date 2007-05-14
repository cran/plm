length.line <- 80

random.method.list <- list(swar="Swamy-Arora",
                          walhus="Wallace-Hussein",
                          amemiya="Amemiya",
                          nerlove="Nerlove")

effect.plm.list <- list(individual="Oneway (individual) effect",
                        time="Oneway (time) effect",
                        twoways="Twoways effects")
effect.pvcm.list <- list(individual="Oneway (individual) effect",
                         time="Oneway (time) effect")
effect.pggls.list <- list(individual="Oneway (individual) effect",
                          time="Oneway (time) effect")
effect.pgmm.list <- list(individual="Oneway (individual) effect",
                         twoways="Twoways effects")

model.plm.list <- list(pooling="Pooling",
                       within="Within",
                       between="Between",
                       random="Random Effect",
                       ht="Hausman-Taylor",
                       fd="First--Difference"
                       )

model.pvcm.list <- list(within="No-pooling model",
                        random="Random coefficients model")
model.pggls.list <- list(within="Within model",
                         random="Random effects model")
model.pgmm.list <- list(onestep="One step model",
                        twosteps="Two steps model")

inst.method.list <- list(baltagi="Baltagi",
                    bvk="Balestra-Varadharajan-Krishnakumar")


robust.list <- list(white1="White 1",
                    white2="White 2",
                    arellano="Arellano")

weights.list <- list(HC0="HC0",
                    HC1="HC1",
                    HC2="HC2",
                    HC3="HC3",
                    HC4="HC4")

oneof <- function(x){
  x <- names(x)
  last <- x[length(x)]
  x <- x[-length(x)]
  x <- paste(x,collapse=", ")
  x <- paste(x,last,sep=" and ")
  x
}
