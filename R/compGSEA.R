#' @param x Named vector of gene-level stats. Names should be the same as in 'genesets'
#' @param genesets An object of class \code{GeneSetCollection} or a list of genesets
#' @param fdr.thres The FDR threshould, genesets with FDR lower than the threshold will be returned
#' @param trim.gs A logical value indicating if the geneset should be trimed according to the names of \code{x}
#' @param nperm	Number of permutations to do. Minimial possible nominal p-value is about 1/nperm
#' @param ... other arguments passed to \code{\link{fgsea}}
#' @return return a list with 3 components: upreg - the upregulated genesets, downreg - the down-regulated 
#' @export
#' @import fgsea
#' @import GSEABase
#' 
#' 
#' 
compGSEA <- function(x, genesets, fdr.thresh = 0.05, trim.gs = FALSE, nperm = 1000, ...) {
  
  if (!inherits(genesets, "list"))
    stop("genesets should be either a \"GeneSetCollection\" or \"list\".")
  
  if (trim.gs)
    genesets <- lapply(genesets, intersect, names(x))
  
  # run gsea
  fgseaRes <- fgsea(pathways = genesets, stats = x, nperm = nperm, ...)
  fgseaRes <- fgseaRes[fgseaRes$padj < fdr.thresh, ]
  fgseaRes <- fgseaRes[order(fgseaRes$pval), ]
  
  posRes <- fgseaRes[fgseaRes$NES > 0]
  negRes <- fgseaRes[fgseaRes$NES < 0]
  
  list(upreg = posRes, downreg = negRes, gs = genesets[fgseaRes$pathway])
}

setGeneric("compGSEA")

setMethod(compGSEA, signature = c(genesets = "GeneSetCollection"), 
          definition = function(x, genesets, fdr.thresh = 0.05, trim.gs = FALSE, ...) {
            
            genesetslist <- geneIds(genesets)
            s <- compGSEA(x, genesetslist, fdr.thresh = fdr.thresh, trim.gs = FALSE, ...)
            gsn <- c(s$upreg$pathway, s$downreg$pathway)
            gss <- lapply(genesets[gsn], function(x) {
              geneIds(x) <- s$gs[[setName(x)]]
              x
            })
            s$gs <- GeneSetCollection(gss)
            s
          })