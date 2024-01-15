############ 
## ssGSEA ##
############

## expression_matrix = expression matrix with colnames as samples and rownames as genes
## pathways = a pathway list
## recycle = if you want to recycle the last functional annotation
## max_lenght = max length of a path to be considered
## min_lenght = min length of a path to be considered
## min_Ratio = the minimum ratio of a path to be considered
## use_enrichr = functionaly enrich the paths present in the dataset
## organism = in case use use_enrichr is TRUE, the organims
## collapse = collapse enriched paths only if the collapsing if use_enrichr is TRUE
## alpha = exponent of the positive ranking genes
## scale_size = normalize by GeneSet size
## scale_ratio = normalize by GeneSet Ratio
## curated_results = create a curated version of the results
## max_pvalue_collapse = max Pvalue to take into account in order to collapse a path into another

ssGSEA <- function(expression_matrix = expression_matrix, pathways = pathways, recycle = F,
                   max_length = 250, min_length = 10, min_Ratio = 0,
                   use_enrichr = F, organism = "org.Hs.eg.db",
                   collapse = F, max_pvalue_collapse = 0.05,
                   alpha = 0.25, scale_size = F, scale_ratio = F, curated_results = F){
  ### GENSET WORK
  # Once the list is done, time to work with the data of the list
  if (exists("pathways_onto_the_analysis", envir = .GlobalEnv) && recycle){
    cat(paste0("Using existing pathsways_enriched object from the environment","\n"))
    pathways_onto_the_analysis <- get("pathways_onto_the_analysis", envir = .GlobalEnv)
    } else {
      cat(paste0("Creating an enriched version of the annotation","\n"))
      paths_nested <- list()
      for (i in 1:length(pathways)){
        ## Extract the genes of the path
        path_now_genes <- pathways[[i]]
        ## Calculate the lenght of the path and the ratio of it
        path_now_ratio <- length(intersect(path_now_genes, rownames(expression_matrix)))/length(path_now_genes)
        # Do the nested list
        paths_nested[[i]] <- list(genes = path_now_genes,
                                  long  = length(intersect(path_now_genes, rownames(expression_matrix))),
                                  ratio = path_now_ratio)
        names(paths_nested)[i] <- names(pathways)[i]
        
        ## Filter the nested list
        # Define the function to filter
        filter_path <- function(enriched_path_list, min_length, max_length, min_ratio){
          filtered_paths <- lapply(enriched_path_list, function(path) {
            if (path$long >= min_length && path$long <= max_length &&
                path$ratio >= min_ratio) {
              return(path)
              } else {
                return(NULL)
                }
            })
          filtered_paths <- filtered_paths[!sapply(filtered_paths, is.null)]
          return(filtered_paths)
          }
        
        # run the function
        pathways_enriched <- filter_path(enriched_path_list = paths_nested,
                                         min_length = min_length, max_length = max_length,
                                         min_ratio = min_Ratio)
        pathways_enriched <<- pathways_enriched
        
        # once the list is filtered, remove the extra data
        filter_the_info <- function(enriched_path_list){
          remove_extra_info <- lapply(enriched_path_list, function(path) {
            path$long <- NULL  # Remove length
            path$ratio <- NULL   # Remove ratio
            return(path)
            })
          return(remove_extra_info)
          }
        pathways_removed <- filter_the_info(enriched_path_list = pathways_enriched)
        cat(paste0("Enriched path: ",i,"/",length(pathways)),"\r")
        }
      cat(paste0("Enriched path annotation can be consulted in the 'pathway_enriched' object","\n"))
      pathways_onto_the_analysis <- list()
      cat("\n")
      cat("Filtering pathways by length and ratio", "\n")
      cat("\n")
      for (j in 1:length(pathways_removed)){
        pathways_onto_the_analysis[[j]] <- unlist(pathways_removed[[j]], use.names = F)
        names(pathways_onto_the_analysis)[j] <- names(pathways_removed)[j]}
      pathways_onto_the_analysis <<- pathways_onto_the_analysis}
  
  if ((use_enrichr)){
    cat(paste0("Functional annotation of the data","\n"))
    #### USE enrichr to curate the data more
    # annotation to data frame
    pathways_data <- do.call(rbind, lapply(names(pathways_onto_the_analysis), function(path) {
      data.frame(path = path, gene = pathways[[path]])
    }))
    colnames(pathways_data) <- c("ID","geneID")
    # Back genes
    backgGenes <- AnnotationDbi::keys(get(organism)) # EntrezID
    backgGenes <- AnnotationDbi::select(get(organism), keys = backgGenes, columns = "UNIPROT") # Converts entrezID to GeneSymbol
    backgGenes <- backgGenes$UNIPROT
    backgGenes_kk <- AnnotationDbi::select(get(organism), keys = backgGenes, columns = "SYMBOL", keytype = "UNIPROT") # Converts entrezID to Symbols
    backgGenes_kk <- backgGenes_kk$SYMBOL
    
    # Execute cluster profiler
    cluster_prof <- clusterProfiler::enricher(gene = c(rownames(expression_matrix)), TERM2GENE = pathways_data,
                                              universe = backgGenes_kk,
                                              pvalueCutoff = 0.05, qvalueCutoff = 0.01,
                                              minGSSize = min_length, maxGSSize = max_length)
    
    # Save clustser profiler data
    cat(paste0("Functional annotation result can be consulted in the 'cluster_prof_table' note that only significantly enriched are annotated","\n"))
    cluster_prof_table <- cluster_prof@result
    cluster_prof_table <- cluster_prof_table %>% 
      filter(p.adjust <= 0.05)
    cluster_prof_table <<- cluster_prof_table
    pathways_onto_the_analysis <- pathways_onto_the_analysis[names(pathways_onto_the_analysis) %in% cluster_prof_table$ID]
    
    if (collapse){
      cat("\n")
      cat(paste0("Collapse of the functional annotation data","\n"))
      # Set the information to use in the collapse
      # gene symbols
      genes <- cluster_prof@gene
      universe <- genes
      # pathways to do the collapse
      paths_to_collapse <- cluster_prof@geneSets
      # functional annotation
      functional_annot <- cluster_prof_table
      functional_annot <- functional_annot %>%
        dplyr::arrange(pvalue)
      
      paths_to_collapse <- paths_to_collapse[functional_annot$ID]
      paths_to_collapse <- lapply(paths_to_collapse, intersect, universe)
      
      # set a parent paths
      parentPaths <- setNames(rep(NA, length(paths_to_collapse)), names(paths_to_collapse))
      for (k in 1:length(paths_to_collapse)) {
        ## set the path to check if is parent
        p <- names(paths_to_collapse)[k]
        if (!is.na(parentPaths[p])) {
          next
        }
        ## set the paths to check if are childs
        paths_to_check <- setdiff(names(which(is.na(parentPaths))),p)
        
        ## Initialize the minPval vector
        minPval <- setNames(rep(1, length(paths_to_check)), paths_to_check)
        
        ## Our universe (u2)
        u2 <- paths_to_collapse[[p]]
        
        ## Do the individual enrichment
        enrich_u2 <- clusterProfiler::enricher(gene = u2, TERM2GENE = pathways_data,
                                               universe = backgGenes_kk,
                                               pvalueCutoff = 0.05, qvalueCutoff = 0.01,
                                               minGSSize = min_length, maxGSSize = max_length)
        
        enrich_u2 <- enrich_u2@result 
        
        ## subset the Pvalues
        enrich_u2_pval <- enrich_u2$pvalue
        names(enrich_u2_pval) <- enrich_u2$ID
        enrich_u2_pval <- enrich_u2_pval[names(enrich_u2_pval) %in% paths_to_check]
        
        ## minPval new
        # common names
        common_minPval <- intersect(names(enrich_u2_pval), names(minPval))
        # minimum of common names
        min_of_common <- pmin(enrich_u2_pval[common_minPval], minPval[common_minPval])
        # find uncommon names
        uncommon_names <- setdiff(names(minPval), names(enrich_u2_pval))
        # get the uncommon values from minPval
        uncommon_minPval <- minPval[uncommon_names]
        minPval <- c(uncommon_minPval,min_of_common)
        parentPaths[names(which(minPval < max_pvalue_collapse))] <- p
        cat(paste0("Now collapsing path:",k,"/",length(paths_to_collapse)),"\r")
        
      }
      cat("\n")
      cat("Collapsed process can be consulted in the 'collapsed_pathways' object")
      cat("\n")
      collapsed_pathways <<- list(mainPathways = names(which(is.na(parentPaths))),
                                  parentPathways = parentPaths)
      pathways_onto_the_analysis <- pathways_onto_the_analysis[collapsed_pathways$mainPathways]
    }}
  cat(paste0("Pathway data processing finished","\n"))
  cat("\n")
  cat(paste0("Working on the expression matrix data","\n"))
  ### EXPRESSION MATRIX WORK
  # Save the expression matrix shape
  expressed_genes <- rownames(expression_matrix)
  sample_names <- colnames(expression_matrix)
  num_genes <- nrow(expression_matrix)
  num_samples <- ncol(expression_matrix)
  
  # Create a ranked matrix
  ### Rank 1 is hte lowest expressed and rank N the highest!
  ranked_matrix <<- matrixStats::colRanks(x = expression_matrix, preserveShape = T, ties.method = "average")
  
  ## This list will have the results of the enrichment
  ssgsea_result <- list()
  ### For each same
  for (sample in 1:length(colnames(ranked_matrix))){
    ### This list will be appended, it will contain the enrichment score of each set of genes in the sample
    sample_list <- list()
    ### Transform the sample row of the ranked matrix to a ranked named vector, like in the classical GSEA
    sample_now <- ranked_matrix[,sample, drop = F]
    sample_now <- as.data.frame(sample_now)
    sample_name <- colnames(sample_now) ### This is quite practical, here the sample name is stored
    colnames(sample_now) <- "rank"
    sample_now$genes <- rownames(sample_now)
    sample_now <- sample_now %>% 
      arrange(rank)
    ranked_genes <- sample_now$rank
    names(ranked_genes) <- sample_now$genes
    ### for each set (pathways_onto_the_analysis) calculate the enrichment
    for (set in 1:length(pathways_onto_the_analysis)){
      set_now_name <- names(pathways_onto_the_analysis)[set]
      set_now_genes <- pathways_onto_the_analysis[[set]]
      
      ## indicator pos and negative
      # indicator pos: indicator that means that the gene is in the gene set and it will go to the positive part of the formula
      # indicator neg: indicator that means that the gene is not in the gene set and it will go to the negative part of the formula
      indicator_pos <<- names(ranked_genes) %in% set_now_genes 
      indicator_neg <<- !indicator_pos
      
      ## positive part of the formula
      rank_alpha <- (ranked_genes*indicator_pos)^alpha
      positive_part <- cumsum(rank_alpha)/sum(rank_alpha)
      
      ## negative part of the formula
      negative_part <- cumsum(indicator_neg)/sum(indicator_neg)
      
      ## difference calculation
      diff <- positive_part - negative_part
      
      ## now in case we want, normalize by gene set size or gene set ratio
      if (scale_size) {
        diff <- diff/length(set_now_genes)}
      if (scale_ratio) {
        ratio_of_set <- pathways_enriched[[set_now_name]]$ratio
        diff <- diff*ratio_of_set}
      
      ## final sum
      diff <- sum(diff)*(-1)
      ## append the sample enrichment score list
      sample_list[[set]] <- diff
      names(sample_list)[set] <- set_now_name
      cat(paste0("Enrichment for Set ",set,"/",length(pathways_onto_the_analysis),"\r"))
    }
    ssgsea_result[[sample]] <- sample_list
    names(ssgsea_result)[sample] <- sample_names[sample]
    cat(paste0("Enrichment for sample ",sample,"/",num_samples," named: ",sample_names[sample],"\n"))
  }
  ### MAKE A MATRIX FROM THE LIST
  # Save the data names
  sample_names <- names(ssgsea_result)
  pathway_names <- unique(unlist(lapply(ssgsea_result, names)))
  
  # Initialize the matrix
  result_matrix <- matrix(0, nrow = length(pathway_names), ncol = length(sample_names),
                          dimnames = list(pathway_names, sample_names))
  # Create the matrix
  for (paths in 1:length(pathway_names)){
    for (samples in 1:length(sample_names)){
      if (pathway_names[paths] %in% names(ssgsea_result[[sample_names[samples]]])){
        result_matrix[paths, samples] <- ssgsea_result[[sample_names[samples]]][[pathway_names[paths]]]
      }
    }
  }
  # Finishing the result matrix
  result_matrix <- as.data.frame(result_matrix)
  cols_ <- colnames(result_matrix)
  result_matrix$path <- rownames(result_matrix)
  
  result_matrix <- result_matrix %>% 
    dplyr::relocate(path, .before = all_of(cols_))
  rownames(result_matrix) <- NULL
  
  if (curated_results){
    cat(paste0("Curating dataset","\n"))
    cols_ <- colnames(result_matrix)
    
    pathways_onto_the_analysis <- pathways[result_matrix$path]
    
    pathways_data <- as.data.frame(tibble(
      path = rep("Not done", times = length(pathways_onto_the_analysis)),
      genes = rep("Not done", times = length(pathways_onto_the_analysis))))
    
    for (i in 1:length(pathways_onto_the_analysis)){
      pathways_data$path[i] <- names(pathways_onto_the_analysis)[i]
      pathways_data$Genes_of_Set[i] <- paste0(pathways_onto_the_analysis[[i]], collapse = ";")}
    
    pathways_data$N_Genes_of_Set <- 0
    for (i in 1:length(rownames(pathways_data))){
      pathways_data$N_Genes_of_Set[i] <- length(unlist(strsplit(pathways_data$Genes_of_Set[i], split = ";"),
                                                       use.names = F))}
    pathways_data$Present_Genes <- pathways_data$Genes_of_Set
    for (i in 1:length(rownames(pathways_data))){
      pathways_data$Present_Genes[i] <- paste0(intersect(unlist(strsplit(pathways_data$Genes_of_Set[i], split = ";")),rownames(expression_matrix)),
                                               collapse = ";")}
    pathways_data$N_Present_Genes <- 0
    for (i in 1:length(rownames(pathways_data))){
      pathways_data$N_Present_Genes[i] <- length(unlist(strsplit(pathways_data$Present_Genes[i], split = ";"),
                                                        use.names = F))}
    pathways_data$GeneRatio <- 0
    for (i in 1:length(rownames(pathways_data))){
      pathways_data$GeneRatio[i] <- pathways_data$N_Present_Genes[i]/pathways_data$N_Genes_of_Set[i]}
    
    pathways_data <- pathways_data %>% 
      subset(select = c(path,Genes_of_Set, Present_Genes,
                        N_Genes_of_Set, N_Present_Genes,
                        GeneRatio))
    
    result_matrix <- merge(result_matrix, pathways_data, by = "path")
  }
  
  colnames(result_matrix)[1] <- "GeneSet"
  return(result_matrix)
}

