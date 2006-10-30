random.method.list <- list(swar="Swamy-Arora",
                          walhus="Wallace-Hussein",
                          amemiya="Amemiya",
                          nerlove="Nerlove")

effect.list <- list(individual="Oneway (individual) effect",
                    time="Oneway (time) effect",
                    twoways="Twoways effects")


inst.method.list <- list(baltagi="Baltagi",
                    bvk="Balestra-Varadharajan-Krishnakumar")

model.list <- list(pooling="Pooling",
                   within="Within",
                   between="Between",
                   random="Random Effect",
                   ht="Hausman-Taylor")

oneof <- function(x){
  x <- names(x)
  last <- x[length(x)]
  x <- x[-length(x)]
  x <- paste(x,collapse=", ")
  x <- paste(x,last,sep=" and ")
  x
}
