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

source("./functions/first_accession.R")
#----

# A Top-Table Type object should had been already created
# Also a meta-data object

### Get the data----
data <- openxlsx::read.xlsx("./raw_data/Clean_Top_Table_k6_ANOVA.xlsx", sheet = 1)
meta_data <- read.table("./raw_data/meta_data.tsv", sep = "\t", dec = ".", header = T)
meta_data$pacient_group <- as.character(meta_data$pacient_group)
prot_cluster <- read.table("./raw_data/clust_prot_k6_ANOVA.tsv",header = T, sep = "\t", dec = ".")

# Create a basic file system
wd <- getwd()
dir.create(file.path(wd,"./raw_data"))
dir.create(file.path(wd,"./results"))
dir.create(file.path(wd,"./plots"))
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

# Visualize the generation of the new protein groups
head(data[,c("Accession","Accession_1")])

# Remove all genes with some NA
data <- na.omit(data)
#----

### Get differentially expressed genes, for 1 vs 1 or more groups---- 
# Subset the expression matrix
expression_matrix <- expression_matrix[rownames(expression_matrix) %in% rownames(data),]
expression_matrix$Accession <- rownames(expression_matrix)

# Transform the expression matrix to long format
expression_matrix_long <- tidyr::pivot_longer(expression_matrix, cols = c(1:56), 
                                              names_to = "sample", values_to = "intensity")

expression_matrix_long <- merge(x = expression_matrix_long, y = meta_data, by = "sample")
expression_matrix_long <- merge(x = expression_matrix_long, y = prot_cluster, by = "Accession")

# Calculate the mean intensity for each group
## for pacient groups
##expression_matrix_long2 <- expression_matrix_long %>% 
##  group_by(pacient_group, Accession) %>% 
##  summarise(mean_intensity=mean(intensity))

## for protein clusters
expression_matrix_long2 <- expression_matrix_long %>% 
  group_by(Accession) %>% 
  summarise(mean_intensity=mean(intensity))

# Calculate z-score for each gene
## for pacient groups
##expression_matrix_long3 <- expression_matrix_long2 %>% 
##  group_by(Accession) %>% 
##  mutate(zscore = scale(mean_intensity, center = TRUE, scale = TRUE)[,1])

## for protein clusters
expression_matrix_long3 <- expression_matrix_long2 %>% 
  mutate(zscore = scale(mean_intensity, center = TRUE, scale = TRUE)[,1])
expression_matrix_long3 <- merge(x = expression_matrix_long3, y = prot_cluster, by = "Accession")

## check-point for pacient groups
##ee <- expression_matrix_long3$mean_intensity[expression_matrix_long3$Accession == "A0FGR8" & expression_matrix_long3$cluster == "3"]
##esd <- sd(expression_matrix_long3$mean_intensity[expression_matrix_long3$Accession == "A0FGR8"])
##eman <- mean(expression_matrix_long3$mean_intensity[expression_matrix_long3$Accession == "A0FGR8"])

## check-point for protein clusters
##ee <- expression_matrix_long3$mean_intensity[expression_matrix_long3$Accession == "A0FGR8"]
##esd <- sd(expression_matrix_long3$mean_intensity)
##eman <- mean(expression_matrix_long3$mean_intensity)
#----




