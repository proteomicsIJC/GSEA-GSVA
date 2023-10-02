### GSEA Analysis

### Libraries and WD----
library(rstudioapi)
setwd(dirname(getActiveDocumentContext()$path))

library(dplyr)
library(tidyr)
library(limma)
library(sva)
library(janitor)
library(seqinr)
library(ggplot2)
library(ggrepel)
library(ggfortify)
library(reshape2)
library(mice)
library(gplots)
library(pheatmap)
library(openxlsx)
library(stringr)
library(reshape2)
library(clusterProfiler)
library(msigdbr)
library(fgsea)
library(org.Hs.eg.db)

source("./functions/first_accession.R")
#----

# A Top-Table Type object should had been already created
### Get the data----
g1_vs_rest <- openxlsx::read.xlsx("./raw_data/TT_up_and_down_1_vs_rest.xlsx", sheet = 1)

# Create a basic file system
wd <- getwd()
dir.create(file.path(wd,"../raw_data"))
dir.create(file.path(wd,"../results"))
dir.create(file.path(wd,"../plots"))
#----


### Filter out unwanted columns and reshape the TT dataframe
rownames(g1_vs_rest) <- g1_vs_rest$Accession

##original_dataframe <- g1_vs_rest

# Filter the dataframe by adj.PValue
g1_vs_rest <- g1_vs_rest %>% 
  filter(adj.P.Val < 0.05)

# Get the wanted columns
# adj.PVAL
adj.pval_g1 <- grep(pattern = "^adj.P.Val", colnames(g1_vs_rest))

# PVAL
pval_g1 <- grep(pattern = "^P.Value", colnames(g1_vs_rest))

# logFC
logfc_g1 <- grep(pattern = "^logFC", colnames(g1_vs_rest))

g1_vs_rest <- g1_vs_rest[,c(c(1,2,3),
                            c(logfc_g1),
                            c(pval_g1,adj.pval_g1))]
#----

### Complete the annotation of the dataframe-----
# Obtain a list consisting of all the accession names in the dataset
g_acs_1 <- strsplit(g1_vs_rest$Accession, "\\;")

# Use the first_accession function
good_names_1 <- first_accession(g_acs_1)

# Introduce the data to the dataset
g1_vs_rest$Accession_1 <- good_names_1
g1_vs_rest <- g1_vs_rest %>%
  relocate(Accession_1, .after = Accession)

# Same for gene_names 
g_acs_1 <- strsplit(g1_vs_rest$Gene.names, "\\;")

good_names_1 <- first_accession(g_acs_1)

g1_vs_rest$Gene.names_1 <- good_names_1
g1_vs_rest <- g1_vs_rest %>%
  relocate(Gene.names_1, .after = Gene.names)

#### MAYBE ADD Uniprotws to complete the annotation in case it is necessayry !!!!

# Remove duplicated Gene_names and NAs
g1_vs_rest <- g1_vs_rest %>% 
  filter(!duplicated(Gene.names_1)) %>% 
  filter(!is.na(Gene.names_1))


#### GSEA analysis----
# For each dataset, create a ranked list, list will be ranked following the formula
# rank = (-log10(pvalue)*logFC)
g1_vs_rest$rank <- (-log10(g1_vs_rest$P.Value))*g1_vs_rest$logFC 

# Sort the dataset by rank
g1_vs_rest <- g1_vs_rest %>% 
  arrange(desc(rank))

# Create a named vector
rnk_g1 <- g1_vs_rest$rank
names(rnk_g1) <- g1_vs_rest$Gene.names_1
rnk_g1 <- sort(rnk_g1, decreasing = TRUE)

# Gene sets
## http://baderlab.org/GeneSets
pathways <- gmtPathways(gmt.file = "../raw_data/Human_GOBP_AllPathways_no_GO_iea_July_03_2023_symbol.gmt")

# Gesea analysis
fgseaRes1 <- fgsea(pathways = pathways, 
                   stats    = rnk_g1,
                   minSize  = 30,
                   maxSize  = 500, gseaParam = 0.5)
#----
