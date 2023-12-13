#############################
### .gmt file formating #####
#############################


#### IMPORTANT ##########################################################################
### This scrtipt will only work if the dbplyr version is less or equal to the 2.3.4 #####
###  in case you update the version please run:
###  remove.packages("dbplyr")
### and then
### devtools::install_version("dbplyr", version = "2.3.4")
### in case you only want to intsall the package run:
### devtools::install_version("dbplyr", version = "2.3.4")
#########################################################################################

### Libraries and WD 
library(rstudioapi)
setwd(dirname(getActiveDocumentContext()$path))

library(msigdb)
library(BiocFileCache)
library(GSEABase)

### Download  the data using and msigdb
gsc <- getMsigdb(version = "7.5", org = "hs")

msigdb_human <- getMsigdb(org = "hs", id = "SYM", version = "7.5")

# add the kegg data
msigdb_human <- appendKEGG(msigdb_human)

### Enter the data
# which collection do we have ?
listCollections(msigdb_human)

# list subcollections
listSubCollections(msigdb_human)

# Get the wanted collections
hallmark <- subsetCollection(gsc = msigdb_human, collection = "h")

c2 <- subsetCollection(gsc = msigdb_human, collection = "c2")
c2_cp <- subsetCollection(gsc = c2, subcollection = c("CP:BIOCARTA", "CP:PID","CP:REACTOME",
                                                         "CP:WIKIPATHWAYS","CP:KEGG"))

c5 <- subsetCollection(gsc = msigdb_human, collection = "c5")
c5_gobp <- subsetCollection(gsc = c5, subcollection = "GO:BP")

### Time to make a list for each 
hallmark_list <- list()
for (i in 1:length(hallmark)){
  set_unique <- hallmark[i]
  genes_of_set_unique <- geneIds(set_unique)
  hallmark_list[[i]] <- genes_of_set_unique[[1]]
  names(hallmark_list)[i] <- names(genes_of_set_unique)
}

c2_cp_list <- list()
for (i in 1:length(c2_cp)){
  set_unique <- c2_cp[i]
  genes_of_set_unique <- geneIds(set_unique)
  c2_cp_list[[i]] <- genes_of_set_unique[[1]]
  names(c2_cp_list)[i] <- names(genes_of_set_unique)
}

c5_cp_list <- list()
for (i in 1:length(c5_gobp)){
  set_unique <- c5_gobp[i]
  genes_of_set_unique <- geneIds(set_unique)
  c5_cp_list[[i]] <- genes_of_set_unique[[1]]
  names(c5_cp_list)[i] <- names(genes_of_set_unique)
}

### curate the file c5_cp_list with the quickGO dataset



