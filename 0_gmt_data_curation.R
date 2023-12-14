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
library(ActivePathways)
library(tidyverse)

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
not_good_bp1 <- readr::read_tsv("./raw_data/QuickGO-annotations-1702502595807-20231213.tsv")
not_good_bp2 <- readr::read_tsv("./raw_data/QuickGO-annotations-1702502639380-20231213.tsv")

not_good_bp <- rbind(not_good_bp1,not_good_bp2)
remove(not_good_bp1);remove(not_good_bp2)

not_good_bp <- not_good_bp %>% 
  subset(select = c(`GO NAME`, SYMBOL)) %>% 
  dplyr::mutate(`GO NAME` = toupper(`GO NAME`)) %>% 
  dplyr::mutate(`GO NAME` = gsub(pattern = " ", replacement = "_", x = `GO NAME`)) %>% 
  dplyr::mutate(`GO NAME` = paste0("GOBP_", `GO NAME`))

not_good_bp_collapse <- not_good_bp %>% 
  dplyr::group_by(`GO NAME`) %>% 
  dplyr::summarise(SYMBOL = paste0(SYMBOL, collapse = ","))

bad_go_bp_list <- list()
for (i in 1:length(rownames(not_good_bp_collapse))){
  genes_of_path <- not_good_bp_collapse$SYMBOL[i]
  genes_of_path <- unlist(strsplit(genes_of_path, split = ","))
  bad_go_bp_list[[i]] <- genes_of_path
  names(bad_go_bp_list)[i] <- not_good_bp_collapse$`GO NAME`[i]
}

### subset c5_cp_list based on the bad_go_bp_list
# Identify shared elements
shared_names <- intersect(names(c5_cp_list), names(bad_go_bp_list))

# Keep non-shared elements or elements with different values
c5_cp_list_good <- c5_cp_list
c5_cp_list_good <- lapply(names(c5_cp_list_good), function(name){
  if (name %in% shared_names){
    shared_elements <- intersect(c5_cp_list_good[[name]], bad_go_bp_list[[name]])
    if (length(shared_elements) > 0){
      c5_cp_list_good[[name]] <- setdiff(c5_cp_list_good[[name]], shared_elements)
    }
  }
  c5_cp_list_good[[name]]
})
names(c5_cp_list_good) <- names(c5_cp_list)

list1 <- c5_cp_list
list2 <- c5_cp_list_good
shared_names <- intersect(names(list1), names(list2))

# Extract elements with different number of elements
different_elements <- lapply(shared_names, function(name) {
  if (length(list1[[name]]) != length(list2[[name]])) {
    list(name = name, list1 = list1[[name]], list2 = list2[[name]])
  }
})

# Remove NULL elements from the list
different_elements <- different_elements[!sapply(different_elements, is.null)]

# View the elements with different number of elements
head(print(different_elements))

### Save the list
paths_to_gmt <- c(hallmark_list, c2_cp_list, c5_cp_list_good)
paths_to_gmt_data <- paths_to_gmt
for (i in 1:length(names(paths_to_gmt_data))){
  paths_to_gmt_data[[i]] <- paste0(paths_to_gmt_data[[i]], collapse = ",")
}

paths_to_gmt_data <- as.data.frame(paths_to_gmt_data)
to_pivot_long <- colnames(paths_to_gmt_data)
paths_to_gmt_data <- paths_to_gmt_data %>% 
  pivot_longer(cols = all_of(to_pivot_long),
               names_to = "paths",
               values_to = "genes")


