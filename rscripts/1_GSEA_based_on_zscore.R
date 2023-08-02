### TMT Standard Workflow

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

source("../functions/first_accession.R")
#----

# A Top-Table Type object should had been already created
# Also a meta-data object

### Get the data----
data <- openxlsx::read.xlsx("../raw_data/Clean_Top_Table_k6_ANOVA.xlsx", sheet = 1)

meta_data <- read.table("../raw_data/meta_data.tsv", sep = "\t", dec = ".", header = T)
meta_data$pacient_group <- as.character(meta_data$pacient_group)

  
# Create a basic file system
wd <- getwd()
dir.create(file.path(wd,"../raw_data"))
dir.create(file.path(wd,"../results"))
dir.create(file.path(wd,"../plots"))
#----


### Filter out unwanted columns and reshape the TT dataframe
rownames(data) <- data$Accession
original_dataframe <- data

# Filter the dataframe by adj.PValue
data <- data %>% 
  filter(adj.P.Val < 0.05)

# Get the wanted columns
# adj.PVAL
adj.pval <- grep(pattern = "^adj.P.Val", colnames(data))
# logFC
logfc <- grep(pattern = "vs", colnames(data))
# Other cols
clusters <- grep(pattern = "cluster", colnames(data))

data <- data[,c(c(1,2,3),
                        c(logfc),
                        c(clusters,adj.pval))]
colnames(data)

# Retrive the expression matrix
expression_matrix <- original_dataframe[,c(5:60)]
#----

### Complete the annotation of the dataframe-----
# Obtain a list consisting of all the accession names in the dataset
g_acs <- strsplit(data$Accession, "\\;")

# Use the first_accession function
good_names <- first_accession(g_acs)

# Introduce the data to the dataset
data$Accession_1 <- good_names
data <- data %>%
  relocate(Accession_1, .after = Accession)

# Same for gene_names 
g_acs <- strsplit(data$Gene.names, "\\;")
good_names <- first_accession(g_acs)
data$Gene.names_1 <- good_names
data <- data %>%
  relocate(Gene.names_1, .after = Gene.names)

# Remove all genes with some NA
data <- na.omit(data)

# Retrieve the annotation
annoation_original <- data[,c(1:5)]
annoation <- data[,c(1,4)] 
#----

### Get differentially expressed genes, for 1 vs 1 or more groups---- 
# Subset the expression matrix
expression_matrix <- expression_matrix[rownames(expression_matrix) %in% rownames(data),]
expression_matrix$Accession <- rownames(expression_matrix)

# Transform the expression matrix to long format
expression_matrix_long <- tidyr::pivot_longer(expression_matrix, cols = c(1:56), 
                                              names_to = "sample", values_to = "intensity")

expression_matrix_long <- merge(x = expression_matrix_long, y = meta_data, by = "sample")

# Calculate the mean intensity for each group
## for pacient groups
expression_matrix_zscore <- expression_matrix_long %>% 
  group_by(pacient_group, Accession) %>% 
  summarise(mean_intensity=mean(intensity))

## for protein clusters
##expression_matrix_zscore <- expression_matrix_long %>% 
##  group_by(Accession) %>% 
##  summarise(mean_intensity=mean(intensity))

# Calculate z-score for each gene
## for pacient groups
expression_matrix_zscore <- expression_matrix_zscore %>% 
  group_by(Accession) %>% 
  mutate(zscore = scale(mean_intensity, center = TRUE, scale = TRUE)[,1])

## for protein clusters
##expression_matrix_zscore <- expression_matrix_zscore %>% 
##  mutate(zscore = scale(mean_intensity, center = TRUE, scale = TRUE)[,1])
##expression_matrix_zscore <- merge(x = expression_matrix_zscore, y = prot_cluster, by = "Accession")

## check-point for pacient groups
ee <- expression_matrix_zscore$mean_intensity[expression_matrix_zscore$Accession == "A0FGR8" & expression_matrix_zscore$pacient_group == "3"]
esd <- sd(expression_matrix_zscore$mean_intensity[expression_matrix_zscore$Accession == "A0FGR8"])
eman <- mean(expression_matrix_zscore$mean_intensity[expression_matrix_zscore$Accession == "A0FGR8"])
(ee - eman)/esd

## check-point for protein clusters
##ee <- expression_matrix_zscore$mean_intensity[expression_matrix_zscore$Accession == "A0FGR8"]
##esd <- sd(expression_matrix_zscore$mean_intensity)
##eman <- mean(expression_matrix_zscore$mean_intensity)
##(ee - eman)/esd

## annotate if z-score is up or down
expression_matrix_zscore <- expression_matrix_zscore %>% 
  mutate(up_or_down = (ifelse(
    zscore > 0, "up","down"
  )))

# Get gene names
expression_matrix_zscore <- merge(expression_matrix_zscore, annoation, by.x = "Accession", by.y = "Accession")

# Split the data by clustering codification <- ADD IN LINE 141
zscore_pac1 <- expression_matrix_zscore[expression_matrix_zscore$pacient_group == "1",]
zscore_pac2 <- expression_matrix_zscore[expression_matrix_zscore$pacient_group == "2",]
zscore_pac3 <- expression_matrix_zscore[expression_matrix_zscore$pacient_group == "3",]
zscore_pac4 <- expression_matrix_zscore[expression_matrix_zscore$pacient_group == "4",]

#### ADD A WAY TO FILTER OUT GENES !!!!!! <- SPEAK THIS TOMORROW !!!!

# Create a named vector and order it by decreasing order
z_score_value_pac1 <- zscore_pac1$zscore
names(z_score_value_pac1) <- zscore_pac1$Gene.names_1
z_score_value_pac1 <- sort(z_score_value_pac1, decreasing = TRUE)

z_score_value_pac2 <- zscore_pac2$zscore
names(z_score_value_pac2) <- zscore_pac2$Gene.names_1
z_score_value_pac2 <- sort(z_score_value_pac2, decreasing = TRUE)

z_score_value_pac3 <- zscore_pac3$zscore
names(z_score_value_pac3) <- zscore_pac3$Gene.names_1
z_score_value_pac3 <- sort(z_score_value_pac3, decreasing = TRUE)

z_score_value_pac4 <- zscore_pac4$zscore
names(z_score_value_pac4) <- zscore_pac4$Gene.names_1
z_score_value_pac4 <- sort(z_score_value_pac4, decreasing = TRUE)
#----


#### GSEA analysis----
HS_hallmark_sets <- msigdbr(species = "Homo sapiens", category = "H")

#### ADD GO AND REACTOME !!!!! 

gsea_results <- GSEA(
  geneList = z_score_value_pac2, # Ordered ranked gene list
  minGSSize = 25, # Minimum gene set size
  maxGSSize = 500, # Maximum gene set set
  pvalueCutoff = 0.05, # p-value cutoff
  eps = 0, # Boundary for calculating the p value
  seed = TRUE, # Set seed to make results reproducible
  pAdjustMethod = "BH", # Benjamini-Hochberg correction
  TERM2GENE = dplyr::select(
    HS_hallmark_sets,
    gs_name,
    gene_symbol
  )
)
#----
