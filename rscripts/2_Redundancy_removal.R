### Redundancy removal after GSEA

#### Apply the last Rscript----
library(rstudioapi)
library(fgsea)
library(org.Hs.eg.db)
setwd(dirname(getActiveDocumentContext()$path))

source("./1_GSEA_based_on_rank.R")
#----

## An example can be found here: https://bioconductor.org/packages/devel/bioc/vignettes/fgsea/inst/doc/fgsea-tutorial.html

### GSEA required data----
# ranked lists are already annotated
head(rnk_g1)
head(rnk_g2)
head(rnk_g3)
head(rnk_g4)

# (.gmt) file
pathways <- gmtPathways(gmt.file = "../raw_data/Human_GOBP_AllPathways_no_GO_iea_July_03_2023_symbol.gmt")

# GSEA results extract tables
gsea1 <- gsea_results1@result 
gsea2 <- gsea_results1@result
gsea3 <- gsea_results1@result
gsea4 <- gsea_results1@result

# Convert to fgsea
rownames(gsea1) <- NULL
gsea1 <- gsea1[,-2]
colnames(gsea1)[1] <- "pathway"
#----

### Collapse GSEA results----
collapsedPathways <- collapsePathways(gsea1[order("pvalue")]["p.adjust" < 0.01], 
                                      pathways, rnk_g1)
#----





