#code used to calculate pearson correlation between chromatin accessibility and gene expression

.libPaths("/users/tfurutan/R/4.3")
library(dplyr)
library(ggplot2)
pdf("/users/tfurutan/filter/plot/bulk_peak_hist.pdf")
data.rna <- readRDS("/data/training_data/RNA_bulk.rds")
data.atac <- readRDS("/users/tfurutan/filter/data/bulk_imputed_peak.rds")
genes.inspect <- rownames(data.rna)
print(length(genes.inspect))
data.rna <- as.matrix(data.rna)
data.atac <- as.matrix(data.atac)

# Check if the dimensions of both datasets match
if (!all(dim(data.rna) == dim(data.atac))) {
  stop("The dimensions of RNA and ATAC data do not match!")
}
print(dim(data.rna))
# Check for any NA values in the data
if (any(is.na(data.rna)) || any(is.na(data.atac))) {
  stop("NA values found in the datasets!")
}

corr.results <- data.frame(matrix(NA, nrow = length(genes.inspect), ncol = 1 ))
rownames(corr.results) <- genes.inspect
colnames(corr.results) <- c("Pearson Correlation")

for (gene in genes.inspect) {
  gene.rna <- data.rna[gene,]
  gene.atac <- data.atac[gene,]
  correlation = cor(gene.rna, gene.atac, method = "pearson")
  corr.results[gene, "Pearson Correlation"] = correlation
}
print(dim(corr.results))

corr.results$`Pearson Correlation` <- as.numeric(corr.results$`Pearson Correlation`)
na_count <- sum(is.na(corr.results$`Pearson Correlation`))
print(paste("Number of NA values in Pearson Correlation:", na_count))
ggplot(corr.results, aes(x = `Pearson Correlation`)) +
  geom_histogram(binwidth = 0.01, fill = "blue", color = "black") +
  labs(title = "Correlation Between Expression and Chromatin Accessibility",
       caption = paste0("genes=",length(genes.inspect), "experiments=", length(colnames(data.rna))),
       x = "Pearson Correlation",
       y = "Frequency") +
  theme_minimal()

write.csv(corr.results, "/users/tfurutan/filter/output/bulk_peak.csv")

bins <- seq(0, 1, by = 0.1)
corr.results$bins <- cut(corr.results$`Pearson Correlation`, breaks = bins, include.lowest = TRUE, right = FALSE)
count_table <- table(corr.results$bins)
write.csv(count_table, "/users/tfurutan/filter/output/bulk_peak_table.csv")