### TMT Standard Workflow

### Libraries and WD----
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
library(rstudioapi)
library(reshape2)
source("./functions/first_accession.R")
setwd(dirname(getActiveDocumentContext()$path))
#----

# A Top-Table Type object should had been already created

### Get the data----
data <- openxlsx::read.xlsx("./raw_data/Clean_Top_Table_k6_ANOVA.xlsx", sheet = 1)

# Create a basic file system
wd <- getwd()
dir.create(file.path(wd,"./raw_data"))
dir.create(file.path(wd,"./results"))
dir.create(file.path(wd,"./plots"))
#----


### Filter out unwanted columns and reshape the TT dataframe
rownames(data) <- data$Accession
data2 <- data

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
expression_matrix <- data2[,c(5:60)]
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




