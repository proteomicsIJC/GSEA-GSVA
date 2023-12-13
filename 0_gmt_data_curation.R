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
hallmark <- subsetCollection(gsc = gsc, collection = "h")
c2_cp <- subsetCollection(gsc = gsc, collection = "c2", subcollection = "CP")
c5_gobp <- subsetCollection(gsc = gsc, collection = "c5", subcollection = "GO:BP")




