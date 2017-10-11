setwd("~/GitHub/omic3plus/")

library(fgsea)
library(mogsa)
library(jsonlite)

db <- readRDS(file = "data/msigdb_v6_symbol.RDS")
gs <- lapply(db, slot, "geneIds")
names(gs) <- sapply(db, slot, "setName")


data("NCI60_4arrays")
x <- NCI60_4arrays$agilent
pc <- prcomp(t(x))
plot(pc)
plot(pc$x[, 1:2])

x1 <- pc$rotation[, 1]
ord <- order(x1)
x1 <- x1[ord]
barplot(x1)

fgseaRes <- fgsea(pathways = gs, 
                  stats = x1,
                  minSize=15,
                  maxSize=500,
                  nperm=10000)
fgseaRes <- fgseaRes[fgseaRes$pval < 0.05, ]
fgseaRes <- fgseaRes[order(fgseaRes$pval), ]
head(fgseaRes)

plotGseaTable(gs[fgseaRes$pathway[1:20]], x1, fgseaRes,
              gseaParam = 0.5)


fgseaRes$leadingEdge


###

ll <- minjaccard(gs)
proxy::dist(fgseaRes$leadingEdge, method = "jaccard")

