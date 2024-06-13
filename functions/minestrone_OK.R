#########################
## minestrone for GSEA ##
#########################

## collapsed_list = list of the collapsing process

minestrone_OK <- function(collapsed_list){
  child <- names(collapsed_list$parentPathways)
  parent <- (collapsed_list$parentPathways)
  
  soup <- as.data.frame(tibble::tibble(
    child = child,
    parent = parent
  ))
  
  for (i in 1:length(rownames(soup))){
    if (is.na(soup$parent[i])){
      soup$parent[i] <- soup$child[i]
    }}
  return(soup)
  }