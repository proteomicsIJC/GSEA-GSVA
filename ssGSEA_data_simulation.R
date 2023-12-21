### Create an expression matrix of two samples
# Set a seed for reproducibility
set.seed(123)

# Define the number of genes and samples
num_genes <- 26  # A to Z
num_samples <- 3

# Create a matrix of random values using runif
expression_matrix <- matrix(runif(num_genes * num_samples), nrow = num_genes)

# Assign row and column names
row_names <- paste0("gene_", LETTERS)
col_names <- paste0("s", 1:num_samples)

rownames(expression_matrix) <- row_names
colnames(expression_matrix) <- col_names

# Print the expression matrix
print(expression_matrix)


### Create am small list of genes I can work with
# Set seed for reproducibility
set.seed(123)

# Function to generate random genes for a pathway
generate_random_genes <- function() {
  genes <- sample(paste0("gene_", LETTERS[1:26]), size = sample(5:10, 1))
  return(genes)
}

# Create a list of 6 random pathways with genes
pathways <- list()
for (i in 1:6) {
  pathway_name <- paste("Pathway", i, sep = "_")
  random_genes <- generate_random_genes()
  pathways[[pathway_name]] <- random_genes
}

### Pass the annotation to a dataframe
###pathways_data <- do.call(rbind, lapply(names(pathways), function(path) {
###  data.frame(path = path, gene = pathways[[path]])
###}))
###colnames(pathways_data) <- c("ID","geneID") 

###ei <- clusterProfiler::enricher(gene = c(rownames(expression_matrix)[1:3], "gene_F","gene_M"), TERM2GENE = pathways_data,
###                          pvalueCutoff = 0.05, qvalueCutoff = 0.1, 
###                          pAdjustMethod = "fdr", minGSSize = 1, maxGSSize = 100)
###ei@result
