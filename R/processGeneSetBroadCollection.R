setwd("~/GitHub/omic3plus/")
library(GSEABase)
library(XML)

dd <- XML::xmlTreeParse("data/msigdb_v6.0.xml")
gs <- dd$doc$children
gs <- gs$MSIGDB

# entrez gene id
ab <- xmlApply(gs, function(x) {
  x <- xmlToList(x, simplify = FALSE)

  if (x[["CATEGORY_CODE"]] == "ARCHIVED") 
    return (NULL)
  gss <- GeneSet(type = "BroadCollection", 
                 setName = x[["STANDARD_NAME"]] ,
                 geneIdType = EntrezIdentifier(), 
                 shortDescription = x[["DESCRIPTION_BRIEF"]], 
                 longDescription = x[["DESCRIPTION_FULL"]], 
                 organism = x[["ORGANISM"]],
                 pubMedIds = x[["PMID"]],
                 urls = x[c("GENESET_LISTING_URL", "EXTERNAL_DETAILS_URL")],
                 contributor = x[c("CONTRIBUTOR_ORG", "CONTRIBUTOR")],
                 creationDate = as.character(Sys.time()),
                 collectionType = BroadCollection(
                   category = tolower(x[["CATEGORY_CODE"]]), 
                   subCategory = x[["SUB_CATEGORY_CODE"]])
  )
  geneIds(gss) <- strsplit(x["MEMBERS_EZID"], ",")[[1]]
  gss
})
ab <- ab[!sapply(ab, is.null)]
gsc <- GeneSetCollection(ab)


# symbol

# entrez gene id
abs <- xmlApply(gs, function(x) {
  x <- xmlToList(x, simplify = FALSE)
  
  if (x[["CATEGORY_CODE"]] == "ARCHIVED") 
    return (NULL)
  gss <- GeneSet(type = "BroadCollection", 
                 setName = x[["STANDARD_NAME"]] ,
                 geneIdType = SymbolIdentifier(), 
                 shortDescription = x[["DESCRIPTION_BRIEF"]], 
                 longDescription = x[["DESCRIPTION_FULL"]], 
                 organism = x[["ORGANISM"]],
                 pubMedIds = x[["PMID"]],
                 urls = x[c("GENESET_LISTING_URL", "EXTERNAL_DETAILS_URL")],
                 contributor = x[c("CONTRIBUTOR_ORG", "CONTRIBUTOR")],
                 creationDate = as.character(Sys.time()),
                 collectionType = BroadCollection(
                   category = tolower(x[["CATEGORY_CODE"]]), 
                   subCategory = x[["SUB_CATEGORY_CODE"]])
  )
  geneIds(gss) <- strsplit(x["MEMBERS_SYMBOLIZED"], ",")[[1]]
  gss
})
abs <- abs[!sapply(abs, is.null)]
gscs <- GeneSetCollection(ab)

saveRDS(gsc, file = "data/msigdb_v6_entrez.RDS")
saveRDS(gscs, file = "data/msigdb_v6_symbol.RDS")

