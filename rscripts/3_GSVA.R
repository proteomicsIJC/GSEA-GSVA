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
setwd(dirname(getActiveDocumentContext()$path))
#----

### Function deffinition----
source("./functions/general")
source("./functions/TMT_MaxQuant")
#----

### Get the data----
peaks <- read.csv2("./raw_data/proteins.csv", header = T, check.names = F, sep = ",", dec = ".")
maxquant <- read.table("./raw_data/proteinGroups.txt", header = T, check.names = F, sep = "\t", dec = ".")

# Meta data
to_get_vect <- read.csv2("./raw_data/to_get_vector.csv", sep = ",", header = T,row.names = 1, check.names = F)
to_get_groups <- read.csv2("./raw_data/to_get_name.csv", sep = ",", header = T)

# Contaminants
cont <- readLines("./raw_data/contaminants.fasta")

# Check if folders exist, if NO create them
wd <- getwd()
dir.create(file.path(wd,"./raw_data"))
dir.create(file.path(wd,"./results"))
dir.create(file.path(wd,"./plots"))
file.remove(file.path(wd,"./results/used_parameters.txt"))
file.create(file.path(wd, "./results/used_parameters.txt"))
#----