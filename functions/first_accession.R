first_accession <- function(charlist){
  empty <- c()
  final_list <- list()
  
  for (j in 1:length(charlist)){
    
    for (i in 1:length(charlist[j])){
      befla <- charlist[[j]][1]
      befla <- c(empty,befla)}
    
    final_list[[j]] <- befla
    the_only_good <- unlist(final_list)
  }
  
  return(the_only_good)
}
