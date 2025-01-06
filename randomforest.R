# Load necessary libraries
library(randomForest)
library(genefilter)

# Assuming your data is stored in a data frame named 'your_data'
# Make sure to replace 'your_data' and column names accordingly

colnames(Genelistvalues2) <- c( "gene","80","90", "100", "110","120","130","80STD","90STD", "100STD", "110STD","120STD","130STD","Slope","RSquare")

# Extract relevant columns (gene, TPM values, Slope, RSquare)
gene_data <- Genelistvalues2[, c("gene", "80", "90", "100", "110", "120", "130", "Slope", "RSquare")]


# Assuming 'filtered_data' is a data frame with the gene name in the 'gene' column
gene_names <- filtered_data$gene  # Extract gene names

# Create gene_expression matrix with gene names
gene_expression <- log2(filtered_data[, c("80", "90", "100", "110", "120", "130")])

# Add the 'gene' column to the gene_expression matrix
gene_expression <- cbind(gene = gene_names, gene_expression)

gene_expression <- cbind(gene = filtered_data$gene, gene_expression)

# Exclude the 'gene' column from the predictor matrix (if needed)
gene_expression <- gene_expression[, -1]

# Now 'gene_expression' contains the gene names and expression values



non_missing_counts <- apply(!is.na(gene_data[, 2:7]), 1, sum)



# Summary statistics of non-missing value counts
summary(non_missing_counts)

# Prepare data for modeling
slope <- filtered_data$Slope
r_square <- filtered_data$RSquare
# Define the response variable (G: Good, N: Not Good)
response <- ifelse(abs(slope) > 1 & r_square > 0.8, "G", "N")


# Ensure response is a vector, not a matrix
response <- as.vector(response)



# Assuming 'filtered_data' is a data frame with the gene name in the 'gene' column
gene_names <- filtered_data$gene  # Extract gene names

# Create gene_expression matrix with gene names
gene_expression <- log2(filtered_data[, c("80", "90", "100", "110", "120", "130")])

# Add the 'gene' column to the gene_expression matrix
gene_expression <- cbind(gene = gene_names, gene_expression)

# Exclude the 'gene' column from the predictor matrix (if needed)
gene_expression <- gene_expression[, -1]


set.seed(1234)
print(Sys.time())  # Print current time
rf_model <- randomForest(x = gene_expression, y = as.factor(response), ntree = 200)
print(Sys.time())  # Print time after model training

length(imp.temp) 



# Look at variable importance

# Assuming 'filtered_data' is a data frame with the gene name in the 'gene' column
gene_names <- filtered_data$gene  # Extract gene names

# Look at variable importance
imp.temp <- abs(rf_model$importance[, ])
t <- order(imp.temp, decreasing = TRUE)

# Plot variable importance
plot(c(1:length(imp.temp)), imp.temp[t], log = 'x', cex.main = 1.5,
     xlab = 'gene rank', ylab = 'variable importance', cex.lab = 1.5,
     pch = 16, main = 'Gene Importance', xlim = c(1, length(imp.temp) + 1))

# Add gene names as labels
text(c(1:length(imp.temp)), imp.temp[t], labels = gene_names[t], pos = 4, cex = 0.7)




# Get subset of expression values for 25 most 'important' genes
gn.imp <- rownames(gene_expression)[t]
gn.7 <- gn.imp[1:7]  # vector of top 25 genes, in order
t <- is.element(rownames(gene_expression), gn.7)
sig.gene_expression <- gene_expression[t, ]



# Get subset of expression values for the top 6 'important' genes
top_genes <- gene_names[t][1:6]  # Extract the top 6 gene names
t_top <- is.element(rownames(gene_expression), top_genes)
sig.gene_expression_top <- gene_expression[t_top, ]

response_top <- response[gene_names %in% top_genes]
