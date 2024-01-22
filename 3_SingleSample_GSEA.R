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

### Import some data expression and pathway data expression
data(c2BroadSets)
class(c2BroadSets)
data(gbm_VerhaakEtAl)
expression_matrix <- exprs(gbm_eset)
expression_matrix <- expression_matrix[,c(1:2), drop = F]

data(brainTxDbSets)
pathways <- as.list(brainTxDbSets)
pathways <- fgsea::gmtPathways("./raw_data/Human_GOBP_AllPathways_no_GO_iea_July_03_2023_symbol.gmt")
pathways <- pathways[1:4000]

remove(c2BroadSets)
remove(gbm_eset)
remove(brainTxDbSets)

### Execute the function
library(org.Hs.eg.db) #### In order that the enrichment works the library has to be loaded onto to environment <3
tictoc::tic()
result_1 <- ssGSEA(expression_matrix = expression_matrix, pathways = pathways, recycle = F,
                 max_length = 100, min_length = 10, min_Ratio = 0,
                 organism = "org.Hs.eg.db", use_enrichr = F, 
                 collapse = F, 
                 alpha = 0.25, scale_size = F, scale_ratio = F, curated_results = F)
tictoc::toc()


### Save the enriched path list
saveRDS(pathways_enriched, file = "./results/path_list_max100_min10_0ratio.rds")

### Load enriched path list
pathways_onto_the_analysis <- readRDS("./results/path_list_max100_min10_0ratio.rds")

tictoc::tic()
result_2 <- ssGSEA(expression_matrix = expression_matrix, pathways = pathways, recycle = T,
                 max_length = 100, min_length = 10, min_Ratio = 0,
                 organism = "org.Hs.eg.db", use_enrichr = F, 
                 collapse = F, 
                 alpha = 0.25, scale_size = F, scale_ratio = F, curated_results = F)
tictoc::toc()


