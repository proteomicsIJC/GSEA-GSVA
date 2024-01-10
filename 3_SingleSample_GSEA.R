########################################
#### ssGSEA function and workflow   ####
########################################
library(Biobase)
library(GSVA)
library(GSVAdata)
library(tidyverse)
library(dplyr)

### source the function

data(c2BroadSets)
class(c2BroadSets)
data(gbm_VerhaakEtAl)
expression_matrix <- exprs(gbm_eset)
expression_matrix <- expression_matrix[,c(1:2)]

data(brainTxDbSets)
pathways <- as.list(brainTxDbSets)
pathways <- fgsea::gmtPathways("./raw_data/Human_GOBP_AllPathways_no_GO_iea_July_03_2023_symbol.gmt")
pathways <- pathways[1:4000]

remove(c2BroadSets)
remove(gbm_eset)
remove(brainTxDbSets)

tictoc::tic()
result <- ssGSEA(expression_matrix = expression_matrix, pathways = pathways,  
                 max_length = 500, min_length = 30, min_Ratio = 0,
                 organism = "org.Hs.eg.db", use_enrichr = T, 
                 collapse = F,
                 alpha = 0.25, scale_size = F, scale_ratio = F, curated_results = T)
tictoc::toc()

