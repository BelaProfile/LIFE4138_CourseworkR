data_A_vs_F <- read.table(file = 'A_vs_F.deseq2.results.tsv', sep = '\t', header = TRUE)

library(ggplot2)
library(pheatmap)

#Set thresholds for filtering
padj_threshold <- 0.05
log2fc_up <- 1
log2fc_down <- -1


#Q1: Summary Statistics
#Analyze data_A_vs_F
data_A_vs_F$Significance <- "Non-significant"

upregulated <- subset(data_A_vs_F, padj < padj_threshold & log2FoldChange > log2fc_up)
cat("Upregulated genes:", nrow(upregulated), "\n")

downregulated <- subset(data_A_vs_F, padj < padj_threshold & log2FoldChange < log2fc_down)
cat("Downregulated genes:", nrow(downregulated), "\n")

cat("\nSummary Statistics for data_A_vs_F\n")
nonsignificant <- subset(data_A_vs_F, !(padj < 0.05 & (log2FoldChange > 1 | log2FoldChange < -1)))
cat("Number of non-significant genes:", nrow(nonsignificant), "\n")

missing_padj <- sum(is.na(data_A_vs_F$padj))
cat("Number of genes with missing padj:", missing_padj, "\n")

total_genes <- nrow(data_A_vs_F)
significant_genes <- nrow(upregulated) + nrow(downregulated)
nonsignificant_genes <- nrow(nonsignificant)
cat("Significant genes (Up + Down):", significant_genes, "\n")

#Print summary statistics
cat("Total genes:", nrow(data_A_vs_F), "\n")
cat("\nP-value Summary:\n")
print(summary(data_A_vs_F$pvalue))
cat("\nLog2 Fold Change Summary:\n")
print(summary(data_A_vs_F$log2FoldChange))




#Volcano Plot
volcano_plot <- ggplot(data_A_vs_F, aes(x = log2FoldChange, y = -log10(pvalue))) +
  geom_point(color = "yellow", alpha = 0.5) + #non-significant genes
  geom_point(data = upregulated, aes(x = log2FoldChange, y = -log10(pvalue)), color = "blue") +
  geom_point(data = downregulated, aes(x = log2FoldChange, y = -log10(pvalue)), color = "brown") +
  geom_hline(yintercept = -log10(padj_threshold), linetype = "dashed", color = "black") +
  geom_vline(xintercept = c(log2fc_up, log2fc_down), linetype = "dashed", color = "red") +
  labs(title = "Volcano Plot: A vs F", x = "Log2 Fold Change", y = "-Log10(p-value)")
print(volcano_plot)

#MA Plot
ma_plot <- ggplot(data_A_vs_F, aes(x = baseMean, y = log2FoldChange)) +
  geom_point(color = "brown", alpha = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  scale_x_log10() +
  labs(title = "MA Plot: A_vs_F", x = "Mean Expression (BaseMean)", y = "Log2 Fold Change")
print(ma_plot)

#Histogram of P-values
pvalue_histogram <- ggplot(data_A_vs_F, aes(x = pvalue)) +
  geom_histogram(bins = 50, fill = "dodgerblue", color = "black", alpha = 1) +
  geom_vline(xintercept = 0.05, linetype = "dashed", color = "red", linewidth = 1, 
             label = "Significance threshold (p = 0.05)") +
  labs(title = "Histogram of p-values: A_vs_F", x = "p-value", y = "Frequency")+
  annotate("text", x = 0.07, y = max(table(cut(data_A_vs_F$pvalue, breaks = 50))) * 0.8, 
    label = "Significance threshold", color = "red", angle = 90, size = 3)
print(pvalue_histogram)

#Heatmap
#Select top 20 upregulated and 20 downregulated genes
top_upregulated <- data_A_vs_F[data_A_vs_F$padj < 0.05 & data_A_vs_F$log2FoldChange > 1, ]
top_upregulated <- top_upregulated[order(top_upregulated$padj), ][1:20, ]

top_downregulated <- data_A_vs_F[data_A_vs_F$padj < 0.05 & data_A_vs_F$log2FoldChange < -1, ]
top_downregulated <- top_downregulated[order(top_downregulated$padj), ][1:20, ]

#Combine top genes
top_genes <- rbind(top_upregulated, top_downregulated)

#Simulate expression matrix (replace with actual expression data if available)
set.seed(123)
expression_matrix <- matrix(rnorm(40 * 6, mean = 10, sd = 2), nrow = 40, ncol = 6) # 40 genes, 6 samples
rownames(expression_matrix) <- top_genes$gene_id
colnames(expression_matrix) <- paste0("Sample_", 1:6)

#Z-score normalization
normalized_matrix <- t(scale(t(expression_matrix)))

pheatmap(
  normalized_matrix,
  color = colorRampPalette(c("blue", "white", "red"))(10),
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  fontsize_row = 8,      # Reduce font size for row names
  fontsize_col = 8,      # Adjust font size for column names
  show_rownames = TRUE,
  main = "Heatmap of Top Differentially Expressed Genes"
)



#Significant gene list
#Remove rows with any NA values in the required columns
filtered_data <- data_A_vs_F[complete.cases(data_A_vs_F[, c("log2FoldChange", "pvalue", "padj")]), ]

#Filter for significant genes
significant_genes <- filtered_data[filtered_data$padj < padj_threshold & abs(filtered_data$log2FoldChange) > 1, ]

#Separate upregulated and downregulated genes
upregulated_genes <- significant_genes[significant_genes$log2FoldChange > 0, ]
downregulated_genes <- significant_genes[significant_genes$log2FoldChange < 0, ]

#Create summary tables
upregulated_table <- upregulated_genes[, c("gene_id", "log2FoldChange", "pvalue", "padj")]
write.csv(upregulated_table, "upregulatedgenes_AvsF.csv", row.names = FALSE)
cat("Preview of Upregulated Genes:\n")
print(head(upregulated_table))

downregulated_table <- downregulated_genes[, c("gene_id", "log2FoldChange", "pvalue", "padj")]
write.csv(downregulated_table, "downregulatedgenes_AvsF.csv", row.names = FALSE)
cat("\nPreview of Downregulated Genes:\n")
print(head(downregulated_table))






### A_vs_D ###



data_A_vs_D <- read.table(file = 'A_vs_D.deseq2.results.tsv', sep = '\t', header = TRUE)


#Set thresholds for filtering same as above
#Q1: Summary Statistics
#Analyze data_A_vs_D
data_A_vs_D$Significance <- "Non-significant"

upregulated <- subset(data_A_vs_D, padj < padj_threshold & log2FoldChange > log2fc_up)
cat("Upregulated genes:", nrow(upregulated), "\n")

downregulated <- subset(data_A_vs_D, padj < padj_threshold & log2FoldChange < log2fc_down)
cat("Downregulated genes:", nrow(downregulated), "\n")

cat("\nSummary Statistics for data_A_vs_D\n")
nonsignificant <- subset(data_A_vs_D, !(padj < 0.05 & (log2FoldChange > 1 | log2FoldChange < -1)))
cat("Number of non-significant genes:", nrow(nonsignificant), "\n")

missing_padj <- sum(is.na(data_A_vs_D$padj))
cat("Number of genes with missing padj:", missing_padj, "\n")

total_genes <- nrow(data_A_vs_D)
significant_genes <- nrow(upregulated) + nrow(downregulated)
nonsignificant_genes <- nrow(nonsignificant)
cat("Significant genes (Up + Down):", significant_genes, "\n")

#Print summary statistics
cat("Total genes:", nrow(data_A_vs_D), "\n")
cat("\nP-value Summary:\n")
print(summary(data_A_vs_D$pvalue))
cat("\nLog2 Fold Change Summary:\n")
print(summary(data_A_vs_D$log2FoldChange))

#Volcano Plot
volcano_plot <- ggplot(data_A_vs_D, aes(x = log2FoldChange, y = -log10(pvalue))) +
  geom_point(color = "yellow", alpha = 0.5) + # non-significant genes
  geom_point(data = upregulated, aes(x = log2FoldChange, y = -log10(pvalue)), color = "blue") +
  geom_point(data = downregulated, aes(x = log2FoldChange, y = -log10(pvalue)), color = "brown") +
  geom_hline(yintercept = -log10(padj_threshold), linetype = "dashed", color = "black") +
  geom_vline(xintercept = c(log2fc_up, log2fc_down), linetype = "dashed", color = "red") +
  labs(title = "Volcano Plot: A vs D", x = "Log2 Fold Change", y = "-Log10(p-value)")
print(volcano_plot)

#MA Plot
ma_plot <- ggplot(data_A_vs_D, aes(x = baseMean, y = log2FoldChange)) +
  geom_point(color = "brown", alpha = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  scale_x_log10() +
  labs(title = "MA Plot: A_vs_D", x = "Mean Expression (BaseMean)", y = "Log2 Fold Change")
print(ma_plot)

#Histogram of P-values
pvalue_histogram <- ggplot(data_A_vs_D, aes(x = pvalue)) +
  geom_histogram(bins = 50, fill = "dodgerblue", color = "black", alpha = 1) +
  geom_vline(xintercept = 0.05, linetype = "dashed", color = "red", linewidth = 1, 
             label = "Significance threshold (p = 0.05)") +
  labs(title = "Histogram of p-values: A_vs_D", x = "p-value", y = "Frequency") +
  annotate("text", x = 0.07, y = max(table(cut(data_A_vs_D$pvalue, breaks = 50))) * 0.8, 
           label = "Significance threshold", color = "red", angle = 90, size = 3)
print(pvalue_histogram)

#Heatmap
#Select top 20 upregulated and 20 downregulated genes
top_upregulated <- data_A_vs_D[data_A_vs_D$padj < 0.05 & data_A_vs_D$log2FoldChange > 1, ]
top_upregulated <- top_upregulated[order(top_upregulated$padj), ][1:20, ]

top_downregulated <- data_A_vs_D[data_A_vs_D$padj < 0.05 & data_A_vs_D$log2FoldChange < -1, ]
top_downregulated <- top_downregulated[order(top_downregulated$padj), ][1:20, ]

#Combine top genes
top_genes <- rbind(top_upregulated, top_downregulated)

#Simulate expression matrix (replace with actual expression data if available)
set.seed(123)
expression_matrix <- matrix(rnorm(40 * 6, mean = 10, sd = 2), nrow = 40, ncol = 6) # 40 genes, 6 samples
rownames(expression_matrix) <- top_genes$gene_id
colnames(expression_matrix) <- paste0("Sample_", 1:6)

#Z-score normalization
normalized_matrix <- t(scale(t(expression_matrix)))

pheatmap(
  normalized_matrix,
  color = colorRampPalette(c("blue", "white", "red"))(10),
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  fontsize_row = 8,
  fontsize_col = 8,
  show_rownames = TRUE,
  main = "Heatmap of Top Differentially Expressed Genes"
)

#Significant Gene List
#Remove rows with any NA values in the required columns
filtered_data <- data_A_vs_D[complete.cases(data_A_vs_D[, c("log2FoldChange", "pvalue", "padj")]), ]

#Filter for significant genes
significant_genes <- filtered_data[filtered_data$padj < padj_threshold & abs(filtered_data$log2FoldChange) > 1, ]

#Separate upregulated and downregulated genes
upregulated_genes <- significant_genes[significant_genes$log2FoldChange > 0, ]
downregulated_genes <- significant_genes[significant_genes$log2FoldChange < 0, ]

#Create summary tables
upregulated_table <- upregulated_genes[, c("gene_id", "log2FoldChange", "pvalue", "padj")]
write.csv(upregulated_table, "upregulatedgenes_AvsD.csv", row.names = FALSE)
cat("Preview of Upregulated Genes:\n")
print(head(upregulated_table))

downregulated_table <- downregulated_genes[, c("gene_id", "log2FoldChange", "pvalue", "padj")]
write.csv(downregulated_table, "downregulatedgenes_AvsD.csv", row.names = FALSE)
cat("\nPreview of Downregulated Genes:\n")
print(head(downregulated_table))

