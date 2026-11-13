library(Seurat)
library(ggplot2)
library(patchwork)

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DESeq2")

library(DESeq2)

##############################################################################
# Source for the code: Alithea genomics: BRB-seq pipeline – downstream analysis best practices ("https://alitheagenomics.com/brb-seq-pipeline-downstream-analysis-best-practices/")

# Set the seed 
set.seed(42)

# Utilize data.table package to read in the brb-seq gene count matrix into R
library(data.table)
data.raw <- fread("C:/RU Medical Epigenomics/Omics_in_Practice/umi_counts.txt", sep="\t", header = T, data.table = F)
rownames(data.raw) <- data.raw[,1]
data.raw <- data.raw[,-1]

# Compute depth per sample
data.meta <- data.frame(samples =colnames(data.raw), depth = colSums(data.raw))

# Change column names to numerical values for clearer distinction in the PCA analysis
colnames(data.raw) <- c('D3_R3','D7_R1','D8_R2',
                        'D4_R3', 'D7_R2', '0h_R2',
                        'D7_R3', 'D8_R1', 'D1_R1',
                        'D3_R1', 'D3_R2', 'D2_R3', 
                        '0h_R3', '6h_R3', 'D1_R2',
                        'D2_R1', 'D4_R2', '6h_R1',
                        'D1_R3', 'D2_R2', 'D4_R1',
                        'D8_R3', '0h_R1', '6h_R2')

# Normalization using limma-voom
library(limma)
data.norm <- voom(counts = data.raw, normalize.method = "quantile")$E

# Perform the PCA
pca <- prcomp(t(data.norm), center = TRUE, scale. = FALSE)
data.pca <- data.frame(pca$x)

# Add sample labels
data.pca$Sample <- colnames(data.raw)

# Calculate % explained variance
var_explained <- 100 * pca$sdev^2 / sum(pca$sdev^2)

# Plot
library(ggplot2)
library(ggrepel)

 # The x-axis is PC1, the y axis is PC2 every point represents one sample, which are all labeled according to their time point and replicate
ggplot(data.pca, aes(x = PC1, y = PC2, color = sub("_.*", "", Sample))) +
  geom_point(shape = 1, size = 4, stroke = 1.2) +
  geom_text_repel(
    data = aggregate(cbind(PC1, PC2) ~ sub("_.*", "", Sample), data.pca, mean),
    aes(x = PC1, y = PC2, label = `sub("_.*", "", Sample)`), 
    segment.color = NA, 
    color = "black",
    nudge_x = 2.5,        
    nudge_y = 2.5,        
    box.padding = 1.3,
    point.padding = 0.5,
    max.overlaps = Inf
  ) +
  xlab(paste0("PC1 (", round(var_explained[1], 1), "%)")) +
  ylab(paste0("PC2 (", round(var_explained[2], 1), "%)")) +
  ggtitle("PCA of BRB-seq samples after voom normalization") +
  scale_color_manual(name = NULL, values = c("#FCC6BB", "#FAA18F", "#F87C63", "#F54927", "#F4320B", "#C82909", "#9C2007", "#701705")) +
  theme_bw()


##############################################################################
###### DESeq analysis ###### Source: https://olvtools.com/en/documents/deseq2
# Read in data again
data.DESeq <- fread("C:/RU Medical Epigenomics/Omics_in_Practice/umi_counts.txt", sep="\t", header = T, data.table = F)
rownames(data.DESeq) <- data.DESeq[,1]
data.DESeq <- data.DESeq[,-1]

# Compute depth per sample
DESeq.meta <- data.frame(samples =colnames(data.DESeq), depth = colSums(data.DESeq))

# Change column names to numerical values for clearer distinction in the PCA analysis
colnames(data.DESeq) <- c('D3_R3','D7_R1','D8_R2',
                        'D4_R3', 'D7_R2', '0h_R2',
                        'D7_R3', 'D8_R1', 'D1_R1',
                        'D3_R1', 'D3_R2', 'D2_R3', 
                        '0h_R3', '6h_R3', 'D1_R2',
                        'D2_R1', 'D4_R2', '6h_R1',
                        'D1_R3', 'D2_R2', 'D4_R1',
                        'D8_R3', '0h_R1', '6h_R2')

# Make coldata dataframe indicating the different conditions of the raw count table
coldata <- data.frame(condition = factor(c('D3', 'D7', 'D8',
                                           'D4', 'D7', '0h',
                                           'D7', 'D8', 'D1',
                                           'D3', 'D3', 'D2', 
                                           '0h', '6h', 'D1',
                                           'D2', 'D4', '6h',
                                           'D1', 'D2', 'D4',
                                           'D8', '0h', '6h')))

# Make a DESeq dataset from the raw count matrix
library(DESeq2)
dds <- DESeqDataSetFromMatrix(countData = data.DESeq, colData = coldata, design = ~ condition)

# Run the DESeq analysis
dds <- DESeq(dds)

### To visualize the volcano plot this source was used as a baseline: https://www.biostars.org/p/282295/ 
# Call ENSG library
library(org.Hs.eg.db)

# specify which timepoints to compare (e.g., D8 versus 0h)
res <- results(dds, contrast = c("condition", "0h", "D8"))
topT <- na.omit(as.data.frame(res)) # Remove genes with NA p-values
topT$gene <- rownames(topT) # Extract gene names from rownames

# Convert ALL ENSG IDs in topT to official gene symbols
clean_ensg <- sub("\\..*", "", topT$gene)
topT$symbol <- mapIds(org.Hs.eg.db, keys = clean_ensg, column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first")
topT$symbol <- ifelse(is.na(topT$symbol), topT$gene, topT$symbol) # Keep ENSG if no symbol exists

# Classify genes as significantly upregulated, downregulated, or not significant
topT$status <- "Not significant"
topT$status[topT$padj < 0.05 & topT$log2FoldChange > 1]  <- "Upregulated"
topT$status[topT$padj < 0.05 & topT$log2FoldChange < -1] <- "Downregulated"

# Select top 15 upregulated and top 15 downregulated genes (ranked by padj)
top_up <- head(subset(topT, status == "Upregulated")[order(subset(topT, status == "Upregulated")$padj), ], 15)
top_down <- head(subset(topT, status == "Downregulated")[order(subset(topT, status == "Downregulated")$padj), ], 15)
top_genes <- rbind(top_up, top_down)

# Isolate target gene of interest (PHB2) and prevent duplicate labeling if it's already in top_genes
tfoi <- subset(topT, symbol == "PHB2")
top_genes <- subset(top_genes, symbol != "PHB2")

# Create the volcano plot (matching the print-friendly style of the PCA)
ggplot(topT, aes(x = log2FoldChange, y = -log10(padj), color = status)) +
  geom_point(alpha = 0.6, size = 1.8) +
  geom_point(data = tfoi, color = "black", fill = "#FFD700", shape = 21, size = 2, stroke = 1.2) + # Highlight dot for PHB2
  geom_text_repel(
    data = top_genes,
    aes(label = symbol),
    color = "black",
    size = 3,
    box.padding = 0.5,
    point.padding = 0.3,
    min.segment.length = 0,
    segment.color = "grey40",
    segment.size = 0.3,
    max.overlaps = Inf
  ) +
  geom_label_repel(
    data = tfoi,
    aes(label = symbol),
    color = "black",
    fontface = "bold",
    size = 3,
    box.padding = 1,
    point.padding = 0.4,
    min.segment.length = 0,
    segment.color = "black",
    segment.size = 0.5
  ) +
  scale_color_manual(name = NULL, values = c("Not significant" = "grey70", "Upregulated" = "#D95F02", "Downregulated" = "#377EB8")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  xlab("Log2 (fold-change)") +
  ylab("-log10(adjusted p-value") +
  ggtitle("Volcano plot: 0h vs D8") +
  theme_bw()


##############################################################################
### to do this for all time points ###:
# Define the timepoints to compare against 0h
timepoints <- c("6h", "D1", "D2", "D3", "D4", "D7", "D8")

# Loop through all timepoints and combine results into one dataframe
all_res <- do.call(rbind, lapply(timepoints, function(tp) {
  res <- results(dds, contrast = c("condition", "0h", tp))
  df <- na.omit(as.data.frame(res))
  df$gene <- rownames(df)
  df$comparison <- paste0("0h vs ", tp)
  return(df)
}))

all_res$comparison <- factor(all_res$comparison, levels = paste0("0h vs ", timepoints))

# Convert ALL ENSG IDs to official gene symbols
clean_ensg <- sub("\\..*", "", all_res$gene)
symbol_map <- mapIds(org.Hs.eg.db, keys = unique(clean_ensg), column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first")
all_res$symbol <- symbol_map[clean_ensg]
all_res$symbol <- ifelse(is.na(all_res$symbol), all_res$gene, all_res$symbol)

# Classify significance per gene
all_res$status <- "Not significant"
all_res$status[all_res$padj < 0.05 & all_res$log2FoldChange > 1]  <- "Upregulated"
all_res$status[all_res$padj < 0.05 & all_res$log2FoldChange < -1] <- "Downregulated"

# Select top 15 upregulated and top 15 downregulated genes PER comparison (cleaner for multi-panel plots)
top_genes <- do.call(rbind, lapply(split(all_res, all_res$comparison), function(df) {
  up <- head(subset(df, status == "Upregulated")[order(subset(df, status == "Upregulated")$padj), ], 15)
  down <- head(subset(df, status == "Downregulated")[order(subset(df, status == "Downregulated")$padj), ], 15)
  return(rbind(up, down))
}))

# Combine top_genes and PHB2 into ONE dataframe so ggrepel prevents overlap between ALL labels
tfoi <- subset(all_res, symbol == "PHB2")
top_genes <- subset(top_genes, symbol != "PHB2")

top_genes$is_tfoi <- FALSE
tfoi$is_tfoi <- TRUE
all_labels <- rbind(top_genes, tfoi)


# Create the multi-panel volcano plot
ggplot(all_res, aes(x = log2FoldChange, y = -log10(padj), color = status)) +
  geom_point(alpha = 0.5, size = 1.0) +
  geom_point(data = tfoi, color = "black", fill = "#FFD700", shape = 21, size = 2.4, stroke = 0.8) +
  geom_text_repel(
    data = top_genes,
    aes(label = symbol),
    color = "black",
    size = 2.1,
    ylim = c(25, NA),      
    box.padding = 0.25,
    point.padding = 0.15,
    min.segment.length = 0,
    segment.color = "grey40",
    segment.size = 0.25,
    max.overlaps = Inf
  ) +
  geom_label_repel(
    data = tfoi,
    aes(label = symbol),
    color = "black",
    fill = "#FFF8DC",
    fontface = "bold",
    size = 2.1,
    label.padding = unit(0.12, "lines"),
    ylim = c(8, 24),       
    box.padding = 0.2,
    point.padding = 0.15,
    min.segment.length = 0,
    segment.color = "black",
    segment.size = 0.35
  ) +
  facet_wrap(~ comparison, nrow = 2) +
  scale_color_manual(name = NULL, values = c("Not significant" = "grey75", "Upregulated" = "#D95F02", "Downregulated" = "#377EB8")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", linewidth = 0.4) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", linewidth = 0.4) +
  xlab(expression(Log~fold~change)) +
  ylab(expression(-Log~padj)) +
  ggtitle("Volcano plots: 0h vs all timepoints") +
  theme_bw() +
  theme(legend.position = "bottom")



##############################################################################
##### Creation of heatmaps with pheatmap #####
# Extract variance-stabilized normalized counts from your DESeq2 object
vsd_mat <- assay(vst(dds, blind = TRUE))

# Reuse your existing symbol_map from the volcano plot step to name the rows
symbols <- symbol_map[sub("\\..*", "", rownames(vsd_mat))]
rownames(vsd_mat) <- ifelse(is.na(symbols), rownames(vsd_mat), symbols)

# Collect the top 30 DEGs (lowest padj) from EACH comparison against 0h
timepoints <- c("6h", "D1", "D2", "D3", "D4", "D7", "D8")

top_per_tp <- unlist(lapply(timepoints, function(tp) {
  res_tp <- na.omit(as.data.frame(results(dds, contrast = c("condition", "0h", tp))))
  res_tp$symbol <- symbol_map[sub("\\..*", "", rownames(res_tp))]
  res_tp$symbol <- ifelse(is.na(res_tp$symbol), rownames(res_tp), res_tp$symbol)
  head(res_tp[order(res_tp$padj), "symbol"], 30)
}))

# Combine all unique top DEGs across all timepoints and include PHB2
selected_genes <- unique(c(top_per_tp, "PHB2"))

# Match against your matrix
present <- intersect(selected_genes, rownames(vsd_mat))

# Compute row z-scores so each gene is shown relative to its own average
heat_mat <- vsd_mat[present, , drop = FALSE]
heat_z <- t(scale(t(heat_mat)))
heat_z <- heat_z[complete.cases(heat_z), , drop = FALSE]

# Sort columns chronologically from 0h to D8 (including R1, R2, R3 order)
time_order <- c("0h", "6h", "D1", "D2", "D3", "D4", "D7", "D8")
sample_time <- factor(sub("_.*", "", colnames(heat_z)), levels = time_order)
col_order <- order(sample_time, colnames(heat_z))

heat_z <- heat_z[, col_order]

# Create custom row labels where only PHB2 is bolded
custom_row_labels <- ifelse(
  rownames(heat_z) == "PHB2",
  expression(bold("PHB2")),
  rownames(heat_z)
)

# Prepare color-coded timepoint annotation matching your PCA palette
ann_col <- data.frame(
  Timepoint = sample_time[col_order],
  row.names = colnames(heat_z)
)

ann_colors <- list(
  Timepoint = c(
    "0h" = "#FCC6BB", "6h" = "#FAA18F", "D1" = "#F87C63", "D2" = "#F54927",
    "D3" = "#F4320B", "D4" = "#C82909", "D7" = "#9C2007", "D8" = "#701705"
  )
)

# Plot the clean, unified heatmap
pheatmap(
  heat_z,
  color = colorRampPalette(c("#67a9cf", "#f7f7f7", "#ef8a62"))(100), 
  breaks = seq(-2, 2, length.out = 101),                            
  border_color = NA, 
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  annotation_col = ann_col,
  annotation_colors = ann_colors,
  show_colnames = TRUE,
  labels_row = custom_row_labels,
  fontsize_row = 6.5, 
  legend_breaks = c(-2, -1, 0, 1, 2),
  legend_labels = c("-2 (Low)", "-1", "0 (Mean)", "1", "2 (High)"),
  main = paste0("Union of top 30 DEGs per timepoint vs 0h and PHB2 (n = ", nrow(heat_z), ")")
)



##############################################################################
### pheatmap for genes of interest ###
# Prepare normalized matrix (reusing symbol_map)
vsd_mat <- assay(vst(dds, blind = FALSE))
symbols <- symbol_map[sub("\\..*", "", rownames(vsd_mat))]
rownames(vsd_mat) <- ifelse(is.na(symbols), rownames(vsd_mat), symbols)

# Define the 3 gene sets (all including PHB2)
gene_sets <- list(
  "1_Nucleus_and_Chromatin_ChIPseq" = c(
    "PHB2", "PHB", "TFCP2", "MBD2", "RNF2", "MYC", "MYOD1", "MEF2A", "MEF2C", "MEF2D", 
    "HDAC1", "DNMT3A", "ESR1", "ESR2", "AR", "KPNA1", "KPNA2", "KPNA3", "KPNA4", 
    "ARFGEF3", "PGRMC1", "HES1", "TP53", "CDKN1A", "SMAD4", "STAT3", "TARDBP", "HNRNPA1"
  ),
  "2_Mitochondria_and_Mitophagy" = c(
    "PHB2", "PHB", "OPA1", "NDUFS1", "DNAJC19", "TAFAZZIN", "HAX1", "P4HB", "PINK1", "PRKN", 
    "OPTN", "CALCOCO2", "AURKA", "GOLPH3", "ULK1", "ZFYVE1", "WIPI1", "PARL", "PGAM5",
    "MAP1LC3A", "MAP1LC3B", "OMA1", "CYCS", "CASP9", "BAX", "BAK1", "STOML2", "LIX1", 
    "PARP1", "AIFM1", "PMAIP1", "BCL2L11", "BBC3", "BCL2", "RCHY1", "TYMS"
  ),
  "3_Signaling_and_Proteomics" = c(
    "PHB2", "AKT1", "HRAS", "KRAS", "NRAS", "RAF1", "BRAF", "MAP2K1", "MAP2K2", "MAPK1", "MAPK3",
    "EGFR", "IGF1R", "ERBB2", "MTOR", "RAD21", "PLK1", "PAK4", "YWHAH", "CDC42", "RSU1",
    "MITF", "RHOA", "DIAPH1", "TNFRSF1A", "CD14", "TLR4", "RACK1", "ITGB1", "PTK2", 
    "WEE1", "CCNB1", "CDC25C", "EIF4E", "ATP7A", "CCND1", "VEGFA", "TGFB1", "CCN2", 
    "IL11", "GALNT14", "IGFBP6", "SLAMF1", "SLAMF7", "ABL1", "PKM", "KEAP1", "NFE2L2", 
    "TNFAIP8L1", "TMBIM6"
  )
)

# Helper function to plot each thematic heatmap
plot_theme_heatmap <- function(genes, plot_title) {
  present <- intersect(unique(genes), rownames(vsd_mat))
  heat_z <- t(scale(t(vsd_mat[present, , drop = FALSE])))
  heat_z <- heat_z[complete.cases(heat_z), , drop = FALSE]
  
  # Chronological column order
  time_order <- c("0h", "6h", "D1", "D2", "D3", "D4", "D7", "D8")
  sample_time <- factor(sub("_.*", "", colnames(heat_z)), levels = time_order)
  col_order <- order(sample_time, colnames(heat_z))
  heat_z <- heat_z[, col_order]
  
  ann_col <- data.frame(Timepoint = sample_time[col_order], row.names = colnames(heat_z))
  ann_colors <- list(Timepoint = c(
    "0h" = "#FCC6BB", "6h" = "#FAA18F", "D1" = "#F87C63", "D2" = "#F54927",
    "D3" = "#F4320B", "D4" = "#C82909", "D7" = "#9C2007", "D8" = "#701705"
  ))
  
  custom_row_labels <- ifelse(rownames(heat_z) == "PHB2", expression(bold("PHB2")), rownames(heat_z))
  
  pheatmap(
    heat_z,
    color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
    breaks = seq(-2, 2, length.out = 101),
    border_color = NA,
    annotation_col = ann_col,
    annotation_colors = ann_colors,
    cluster_rows = TRUE,
    cluster_cols = FALSE,
    show_colnames = TRUE,
    labels_row = custom_row_labels,
    fontsize_row = 8.5,
    legend_breaks = c(-2, -1, 0, 1, 2),
    legend_labels = c("-2 (Low)", "-1", "0 (Mean)", "1", "2 (High)"),
    main = paste0(plot_title, " (n = ", nrow(heat_z), ")")
  )
}

# Run each plot one by one (click the blue arrows in RStudio Plots pane to switch between them)
plot_theme_heatmap(gene_sets[[1]], "Nucleus & Chromatin Regulators (ChIP-seq candidates)")
plot_theme_heatmap(gene_sets[[2]], "Mitochondria, Mitophagy & Apoptosis")
plot_theme_heatmap(gene_sets[[3]], "Signaling, Cell Cycle & Disease Regulators")


