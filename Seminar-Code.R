# ---------------------------------------------------------
# 0. Load Packages
# ---------------------------------------------------------
suppressPackageStartupMessages({
  library(dplyr)
  library(qvalue)
  library(ggplot2)
})

set.seed(999)  # for reproducibility

# ---------------------------------------------------------
# 1. Global Simulation Parameters
# ---------------------------------------------------------
nGenes        <- 10000         # total number of features (genes)
propDE        <- 0.30          # proportion of truly DE genes (30% is quite high)
nSamplesGroup <- 10            # samples per group
muBase        <- 100           # baseline for gamma distribution
dispersion    <- 0.05          # NB dispersion parameter (lower means less overdispersion, easier DE detection)
fcLocation    <- log2(2.0)     # location of log2 fold-change for DE genes (center ~ 2-fold)
fcScale       <- 0.3           # spread in log2 fold-change

# We'll create group A and B library sizes (some variation around ~80k reads).
librarySizesA <- rpois(nSamplesGroup, 80000)
librarySizesB <- rpois(nSamplesGroup, 82000)

# ---------------------------------------------------------
# 2. Specify Which Genes Are Truly DE
# ---------------------------------------------------------
isDE <- rep(FALSE, nGenes)
deCount <- round(propDE * nGenes)
deIndices <- sample(seq_len(nGenes), size = deCount, replace = FALSE)
isDE[deIndices] <- TRUE

# ---------------------------------------------------------
# 3. Generate Baseline Means
#    - first from a gamma distribution for diversity
# ---------------------------------------------------------
# shape, rate chosen to yield a typical range around 'muBase'
baselineMeans <- rgamma(nGenes, shape = 2, rate = 2) * muBase

# ---------------------------------------------------------
# 4. Generate Fold Changes for DE Genes
#    - sample from a normal on log2 scale, then exponentiate
# ---------------------------------------------------------
log2fcVals <- rnorm(deCount, mean = fcLocation, sd = fcScale)
fcVals <- 2^log2fcVals  # exponentiate to get actual fold changes

# put them into a vector for all genes
foldChanges <- rep(1, nGenes)
foldChanges[isDE] <- fcVals

# means for group A vs group B
meansA <- baselineMeans
meansB <- baselineMeans * foldChanges

# ---------------------------------------------------------
# 5. Simulate NB Counts for Each Gene in Each Group
# ---------------------------------------------------------
simA <- matrix(0, nrow = nGenes, ncol = nSamplesGroup)
simB <- matrix(0, nrow = nGenes, ncol = nSamplesGroup)

# function to convert dispersion -> size param in NB
# size = 1/dispersion, mu=...
nbSize <- 1 / dispersion  

for (i in seq_len(nSamplesGroup)) {
  # scale means by library sizes
  muA_i <- meansA * (librarySizesA[i] / mean(librarySizesA))
  muB_i <- meansB * (librarySizesB[i] / mean(librarySizesB))
  
  simA[, i] <- rnbinom(nGenes, size = nbSize, mu = muA_i)
  simB[, i] <- rnbinom(nGenes, size = nbSize, mu = muB_i)
}

# combine into single matrix for convenience
countData <- cbind(simA, simB)
groupFactor <- factor(c(rep("A", nSamplesGroup), rep("B", nSamplesGroup)))

# ---------------------------------------------------------
# 6. Quick Differential Expression Test (Wilcoxon on log2-counts)
#    - Real usage might use edgeR or DESeq2
# ---------------------------------------------------------
pseudoCount <- 1
logCountsA  <- log2(simA + pseudoCount)
logCountsB  <- log2(simB + pseudoCount)

pVals <- sapply(seq_len(nGenes), function(g) {
  wilcox.test(logCountsA[g, ], logCountsB[g, ], exact = FALSE)$p.value
})

# approximate effect size (log2 FC)
foldChangeEst <- rowMeans(logCountsB) - rowMeans(logCountsA)

# ---------------------------------------------------------
# 7. Multiple Testing Corrections
# ---------------------------------------------------------
pAdjBH   <- p.adjust(pVals, method = "BH")   # Benjamini-Hochberg
pAdjBY   <- p.adjust(pVals, method = "BY")   # Benjamini-Yekutieli
qVals    <- qvalue(pVals)$qvalues            # Storey Q-values

alpha <- 0.05
sigBH     <- (pAdjBH <= alpha)
sigBY     <- (pAdjBY <= alpha)
sigStorey <- (qVals  <= alpha)

# ---------------------------------------------------------
# 8. Confusion Matrix & Performance Metrics
# ---------------------------------------------------------
# BH
TP_BH <- sum(sigBH &  isDE)
FP_BH <- sum(sigBH & !isDE)
TN_BH <- sum(!sigBH & !isDE)
FN_BH <- sum(!sigBH &  isDE)

# BY
TP_BY <- sum(sigBY &  isDE)
FP_BY <- sum(sigBY & !isDE)
TN_BY <- sum(!sigBY & !isDE)
FN_BY <- sum(!sigBY &  isDE)

# Storey
TP_Storey <- sum(sigStorey &  isDE)
FP_Storey <- sum(sigStorey & !isDE)
TN_Storey <- sum(!sigStorey & !isDE)
FN_Storey <- sum(!sigStorey &  isDE)

# Type I Error: FP / total true null
typeI_BH     <- FP_BH / sum(!isDE)
typeI_BY     <- FP_BY / sum(!isDE)
typeI_Storey <- FP_Storey / sum(!isDE)

# FDR: FP / (TP + FP)
fdr_BH       <- FP_BH / max(1, (TP_BH + FP_BH))
fdr_BY       <- FP_BY / max(1, (TP_BY + FP_BY))
fdr_Storey   <- FP_Storey / max(1, (TP_Storey + FP_Storey))

# Power: TP / total true DE
power_BH     <- TP_BH / sum(isDE)
power_BY     <- TP_BY / sum(isDE)
power_Storey <- TP_Storey / sum(isDE)

resultsDF <- data.frame(
  Method = c("BH","BY","StoreyQ"),
  TypeI  = c(typeI_BH, typeI_BY, typeI_Storey),
  FDR    = c(fdr_BH,   fdr_BY,   fdr_Storey),
  Power  = c(power_BH, power_BY, power_Storey),
  TP     = c(TP_BH,    TP_BY,    TP_Storey),
  FP     = c(FP_BH,    FP_BY,    FP_Storey),
  TN     = c(TN_BH,    TN_BY,    TN_Storey),
  FN     = c(FN_BH,    FN_BY,    FN_Storey)
)

cat("\n===== Performance Metrics =====\n")
print(resultsDF)

# ---------------------------------------------------------
# 9. Plots
# ---------------------------------------------------------
# (a) Volcano plot
volcanoDF <- data.frame(
  gene = 1:nGenes,
  log2FC_est = foldChangeEst,
  pvalue     = pVals,
  isDE       = isDE
)

ggplot(volcanoDF, aes(x = log2FC_est, y = -log10(pvalue), color = isDE)) +
  geom_point(alpha = 0.5) +
  theme_bw() +
  labs(title = "Volcano Plot (Approx)",
       x = "Estimated log2 Fold Change (B vs A)",
       y = "-log10(p-value)") +
  scale_color_manual(values = c("FALSE"="gray50","TRUE"="red"))

# (b) Confusion Matrix Heatmap for BH (example)
confMat_BH <- data.frame(
  category = factor(c("FN","TN","FP","TP"), levels=c("FN","TN","FP","TP")),
  value    = c(FN_BH, TN_BH, FP_BH, TP_BH)
)

ggplot(confMat_BH, aes(x="BH", y=category, fill=value)) +
  geom_tile() +
  geom_text(aes(label = value)) +
  scale_fill_gradient(low="white", high="blue") +
  theme_minimal() +
  labs(title="BH Method: Confusion Matrix Counts", x="Method", y="Category")

# (c) Boxplot of library sizes
libDF <- data.frame(
  libSize = c(librarySizesA, librarySizesB),
  group   = factor(c(rep("A", nSamplesGroup), rep("B", nSamplesGroup)))
)
ggplot(libDF, aes(x=group, y=libSize, fill=group)) +
  geom_boxplot() +
  theme_bw() +
  labs(title="Library Size Distribution")

# (d) Histogram of p-values
pDF <- data.frame(pvalue = pVals)
ggplot(pDF, aes(x=pvalue)) +
  geom_histogram(bins = 50, fill="lightblue", color="black") +
  theme_bw() +
  labs(title="Histogram of Raw P-values", x="p-value", y="Count")

## MA PLOT
meanLogCounts <- (rowMeans(logCountsA) + rowMeans(logCountsB)) / 2
log2FC_est    <- rowMeans(logCountsB) - rowMeans(logCountsA)

maDF <- data.frame(
  A    = meanLogCounts,
  M    = log2FC_est,
  isDE = isDE
)

ggplot(maDF, aes(x = A, y = M, color = isDE)) +
  geom_point(alpha=0.6) +
  labs(title = "MA Plot", x = "Average log2 Expression (A)", y = "log2 Fold Change (M)") +
  theme_bw()

##ROC
library(pROC)

# We have pVals from the test, isDE is the ground truth
# We'll do an example with BH corrected p-values, or raw p-values.
roc_obj <- roc(isDE, 1 - pVals)  # using "1 - pVals" so high = more likely DE
# Alternatively: roc_obj <- roc(isDE, -log10(pVals))

plot(roc_obj, col="blue", main="ROC Curve (Wilcoxon test)")
auc(roc_obj)

##PCA
allLogCounts <- cbind(logCountsA, logCountsB)  # genes x samples

# Perform PCA on transposed data (samples in rows)
pca <- prcomp(t(allLogCounts), scale. = TRUE)

# Make a df for plotting
pcaDF <- data.frame(
  PC1 = pca$x[, 1],
  PC2 = pca$x[, 2],
  group = rep(c("A","B"), each = nSamplesGroup)
)

ggplot(pcaDF, aes(x=PC1, y=PC2, color=group)) +
  geom_point(size=3) +
  theme_bw() +
  labs(title="PCA plot of Samples",
       x = "PC1", y = "PC2")


# If you only want top 100 DE genes
topGenes <- order(pVals, decreasing = FALSE)[1:100]  # sorted by p-value
heatmapData <- allLogCounts[topGenes, ]  # subset

# Option 1: base R
heatmap(as.matrix(heatmapData))

# Option 2: with ComplexHeatmap package
# install.packages("ComplexHeatmap")
# library(ComplexHeatmap)
# Heatmap(heatmapData, name = "log2 count",
#        row_names_gp = gpar(fontsize=4),
#        column_split = rep(c("A","B"), each=nSamplesGroup))


counts_long <- data.frame(
  value = c(as.vector(logCountsA), as.vector(logCountsB)),
  group = rep(c("A","B"), each = nGenes * nSamplesGroup)
)

ggplot(counts_long, aes(x=group, y=value, fill=group)) +
  geom_violin(trim=FALSE) +
  theme_bw() +
  labs(title="Violin Plot of log2(CPM/Counts)", y = "log2(count + 1)")


# 'foldChanges' is the true FC, 'foldChangeEst' is the approximate log2 FC from the data
compDF <- data.frame(
  trueFC_log2 = log2(foldChanges),
  estFC_log2  = foldChangeEst,
  isDE        = isDE
)

ggplot(compDF, aes(x=trueFC_log2, y=estFC_log2, color=isDE)) +
  geom_point(alpha=0.5) +
  geom_abline(slope=1, intercept=0, linetype="dashed", color="red") +
  theme_bw()+
  labs(title="Comparison of True vs Estimated log2 FC",
       x="True log2 FC", y="Estimated log2 FC")

#######
# END #
#######
