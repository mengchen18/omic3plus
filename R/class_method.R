setClass(
  "evaluatedGeneSet", 
  contains = "GeneSet",
  slots = c(pvalue="numeric", 
            score="numeric", 
            cluster = "integer", 
            regulation = "character"),
  
)



setClass(
  "evaluatedGeneSetCollection", 
  contains = "list",
  validity = function(object) {
    msg <- NULL
    if (!all(sapply(object, is, "evaluatedGeneSet")))
      msg <- c(msg, "members must all be 'evaluatedGeneSet' classes")
    tryCatch({
      if (anyDuplicated(names(object)))
        msg <- c(msg, "each setName must be distinct")
    }, error=function(err) {
      msg <<- c(msg, conditionMessage(err))
    })
    if (!is.null(msg))
      msg
    else
      TRUE
  })

setClass("annotatedVector", 
         slots = c(name = "character", 
                   value = "numeric", 
                   annotated = "logical"
         ))

setClass("gsAnnotatedVector", 
         slots = c(stats = "annotatedVector", 
                   annot = "evaluatedGeneSet")
)


gs1 <- new("evaluatedGeneSet", setName = "t1")
gs2 <- new("evaluatedGeneSet", setName = "t2")

gsc <- evaluatedGeneSetCollection(gs1, gs2)
gsc[[1]]@pvalue
class(gsc)

names(gsc)



a <- new("evaluatedGeneSet")


a
setClass()


library(GSEABase)




EnzymeIdentifier()


GeneSet()
## Gene set from ExpressionSet
data(sample.ExpressionSet)
gs1 <- GeneSet(sample.ExpressionSet[100:109])
gs1


gs1@organism
