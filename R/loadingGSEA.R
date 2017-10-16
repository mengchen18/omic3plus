setwd("~/GitHub/omic3plus/")

library(fgsea)
library(mogsa)
library(jsonlite)
library(GSEABase)

data("NCI60_4arrays")
db <- readRDS(file = "data/msigdb_v6_symbol.RDS")

x <- NCI60_4arrays$agilent
pc <- prcomp(t(x))
plot(pc)
plot(pc$x[, 1:2])

x1 <- pc$rotation[, 1]
ord <- order(x1)
x1 <- x1[ord]
barplot(x1)

s <- compGSEA(x = x1, genesets = db, fdr.thresh = 0.5, nperm = 10000, 
              nproc = 20, minSize = 5, maxSize = 500, trim.gs = TRUE)

a <- s$upreg
b <- s$downreg
c <- s$gs




function(sideRes, side = "upregulated") {
  # gsl <- db[sideRes$pathway]
  gsl <- gsi[sideRes$pathway]
  ds <- minjaccard(gsl)
  hcl <- hclust(ds, method = "complete")
  plot(hcl)
  
  cls <- cutree(hcl, h = 0.25)
  cls <- formatC(cls, width=nchar(max(cls)), flag="0")
  cls <- paste("gsCluster", cls, sep = "_")
  
  ##
  lsl <- lapply(1:length(gsl), function(i) {
    gs1 <- gsl[[i]]
    list(
      name = setName(gs1), 
      ids = geneIds(gs1), 
      des = description(gs1), 
      longDes = longDescription(gs1), 
      organism = organism(gs1), 
      pubmedIds = pubMedIds(gs1), 
      idType = geneIdType(gs1),
      listUrl = urls(gs1)[["GENESET_LISTING_URL"]],
      extUrl = urls(gs1)[["EXTERNAL_DETAILS_URL"]], 
      parent = cls[i],
      pvalue = fgseaRes$pval[i] 
    )
  })
  clsLabel <- sapply(lsl, "[[", "parent")
  
  clsl <- lapply(unique(cls), function(i) {
    list(
      children = lsl[i == clsLabel],
      name = i,
      parent = "upregulated"
    )
  })
  
  upreg <- list(
    children = clsl,
    name = "upregulated",
    parent = "total"
  )
}




#######################################
source("R/minjaccard.R")

gss <- gs[fgseaRes$pathway]
barplot(sapply(gss, length))
a <- minjaccard(gss)
am <- as.matrix(a)
which(am < 0.25 & am >0 , arr.ind = TRUE)




plotGseaTable(gs[fgseaRes$pathway[1:20]], x1, fgseaRes,
              gseaParam = 0.5)
fgseaRes$leadingEdge


