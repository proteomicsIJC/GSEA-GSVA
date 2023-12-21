######################################################
#### ssGSEA function and workflow in development ####
#####################################################
library(Biobase)
library(GSVA)
library(GSVAdata)

data(c2BroadSets)
class(c2BroadSets)
data(gbm_VerhaakEtAl)
expression_matrix <- exprs(gbm_eset)
expression_matrix <- expression_matrix[,c(1:5)]

data(brainTxDbSets)
pathways <- as.list(brainTxDbSets)
pathways <- fgsea::gmtPathways("./raw_data/Human_GOBP_AllPathways_no_GO_iea_July_03_2023_symbol.gmt")
pathways <- pathways[1:4000]

remove(c2BroadSets)
remove(gbm_eset)
remove(brainTxDbSets)

##### HERE THE FUNCTION WOULD START THE WORK <3
ssGSEA <- function(expression_matrix = expression_matrix, pathways = pathways,
                   max_length, min_length, min_Ratio,
                   organism = "org.Hs.eg.db", use_enrichr = F, 
                   collapse = F,
                   alpha = 0.25){
  ### GENSET WORK
  # Once the list is done, time to work with the data of the list
  cat(paste0("Creating an enriched version of the annotation","\n"))
  paths_nested <- list()
  
  for (i in 1:length(pathways)){
  ## Extract the genes of the path
  path_now_genes <- pathways[[i]]
  ## Calculate the lenght of the path and the ratio of it
  path_now_length <- length(path_now_genes)
  path_now_ratio <- length(intersect(path_now_genes, rownames(expression_matrix)))/length(path_now_genes)
  
  # Do the nested list
  paths_nested[[i]] <- list(genes = path_now_genes,
                            long  = path_now_length,
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
  
  if ((use_enrichr)){
    cat(paste0("Functional annotation of the data","\n"))
  #### USE enrichr to curate the data more
  # annotation to data frame
    pathways_data <- do.call(rbind, lapply(names(pathways_onto_the_analysis), function(path) {
      data.frame(path = path, gene = pathways[[path]])
      }))
    colnames(pathways_data) <- c("ID","geneID")
    # Back genes
    organism <- get(organism)
    backgGenes <- keys(organism) # EntrezID
    backgGenes <- AnnotationDbi::select(organism, keys = backgGenes, columns = "UNIPROT") # Converts entrezID to GeneSymbol
    backgGenes <- backgGenes$UNIPROT
    backgGenes_kk <- AnnotationDbi::select(organism, keys = backgGenes, columns = "SYMBOL", keytype = "UNIPROT") # Converts entrezID to Symbols
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
        parentPaths[names(which(minPval < 0.05))] <- p
        cat(paste0("Now collapsing path:",k,"/",length(paths_to_collapse)),"\r")
        
      }
      cat("\n")
      cat("Collapsed process can be consulted in the 'collapsed_pathways' object")
      cat("\n")
      collapsed_pathways <<- list(mainPaths = names(which(is.na(parentPaths))),
                                 parent_paths = parentPaths)
      pathways_onto_the_analysis <- paths_to_collapse[names(pathways_onto_the_analysis) %in% collapsed_pathways$mainPaths]
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
  
  # Ranking of the genes
  cat(paste0("Ranked expression matrix can be consulted in the 'ranked_pathways' object","\n"))
  ranked_matrix <<- matrixStats::colRanks(expression_matrix, preserveShape = T)
  # Enrichment Score (ES) for each sample (column)
  ssgsea_result_list <- list()
  for (sample in 1:length(sample_names)){
    ssgsea_now <- list()
    sample_now <- sort(ranked_matrix[,sample_names[sample]], decreasing = T)
    
    for (set in 1:length(pathways_onto_the_analysis)){
      ## HERE APPEND ssgsea_now !!
      ## ALSO ADD NAMES will be the sets
      cat(paste0("Enrichment for Set ",set,"/",length(pathways_onto_the_analysis),"\r"))
    }
    ## HERE APPEND ssgsea_result_list !!
    ## ALSO ADD THE NAMES (will be the samples)
    cat(paste0("Enrichment for sample ",sample,"/",num_samples," named: ",sample_names[sample],"\n"))
  }

  return(sample_now)
  }

tictoc::tic()
result <- ssGSEA(expression_matrix = expression_matrix, pathways = pathways,  
                 max_length = 100, min_length = 1, min_Ratio = 0.2, 
                 organism = "org.Hs.eg.db", use_enrichr = F, collapse = F,
                 alpha = 0.25)
tictoc::toc()


