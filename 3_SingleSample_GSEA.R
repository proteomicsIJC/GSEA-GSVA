########################################
#### ssGSEA function and workflow   ####
########################################
library(Biobase)
library(GSVA)
library(GSVAdata)
library(tidyverse)
library(dplyr)

### source the function
source("./functions/ssGSEA.R")

data(c2BroadSets)
class(c2BroadSets)
data(gbm_VerhaakEtAl)
expression_matrix <- exprs(gbm_eset)
expression_matrix <- expression_matrix[,c(1:2), drop = F]

data(brainTxDbSets)
pathways <- as.list(brainTxDbSets)
pathways <- fgsea::gmtPathways("./raw_data/Human_GOBP_AllPathways_no_GO_iea_July_03_2023_symbol.gmt")
pathways <- pathways[1:400]

remove(c2BroadSets)
remove(gbm_eset)
remove(brainTxDbSets)

library(org.Hs.eg.db) #### In order that the enrichment works the library has to be loaded onto to environment <3
tictoc::tic()
result <- ssGSEA(expression_matrix = expression_matrix, pathways = pathways, recycle = F,
                 max_length = 100, min_length = 10, min_Ratio = 0,
                 organism = "org.Hs.eg.db", use_enrichr = F, 
                 collapse = F, 
                 alpha = 0.25, scale_size = F, scale_ratio = F, curated_results = F)
tictoc::toc()


