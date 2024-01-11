#########################
## minestrone for GSEA ##
#########################

## collapsed_list = list of the collapsing process

minestrone_for_GSEA <- function(collapsed_list, clean.names = F){
  ### make the soup
  soup <- as.data.frame(collapsed_list$parentPathways)
  soup$child <- rownames(soup)
  colnames(soup)[1] <- "parent"
  for (i in 1:length(rownames(soup))){
    if (is.na(soup$parent[i])){
      soup$parent[i] <- soup$child[i] 
    }
  }
  if (clean.names){
  # clean parent names
  soup <- soup %>% 
    mutate(parent = gsub("%.*%","",parent)) %>% 
    mutate(parent = gsub("GO\\:","",parent)) %>% 
    mutate(parent = gsub("[0-9]{5,7}$","",parent)) %>% 
    mutate(parent = gsub("R\\-HSA\\-[0-9]{5,7}\\.[0-9]{1}$","",parent)) %>%
    mutate(parent = gsub("HALLMARK\\_","",parent)) %>% 
    mutate(parent = gsub("\\_"," ",parent))
  
  # clean child names
  soup <- soup %>% 
    mutate(child = gsub("%.*%","",child)) %>% 
    mutate(child = gsub("GO\\:","",child)) %>% 
    mutate(child = gsub("[0-9]{5,7}$","",child)) %>% 
    mutate(child = gsub("R\\-HSA\\-[0-9]{5,7}\\.[0-9]{1}$","",child)) %>%
    mutate(child = gsub("HALLMARK\\_","",child)) %>% 
    mutate(child = gsub("\\_"," ",child))
  }
  return(soup)
}
