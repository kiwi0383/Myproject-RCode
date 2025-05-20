#============================================================================
# Figure 2: Close Relationship Between LM.Sig and Immune Landscape in LUSC
#============================================================================

# Load necessary libraries
library(tibble)
library(pheatmap)
library(RColorBrewer)
library(utils)
library(openxlsx)
library(ggpubr)
library(ggExtra)
library(vegan)
library(dplyr)
library(purrr)
library(circlize)
library(corrplot)
library(Hmisc)
library(rstatix)

# Set working directory (adjust as needed)
setwd("E:/项目一")

#================================================
# Load Data
#================================================

# Load LUSC data
lusc <- read.csv('E:/项目一/data/fig2/lusc.csv', header = TRUE) %>%
  column_to_rownames(var = 'X')

# Load risk score data
riskscore <- openxlsx::read.xlsx('E:/项目一/data/fig1/df.xlsx') %>%
  select(2:14)

# Load ESTIMATE scores
estimate <- read.table("E:/TCGA/TCGA_Annovar/TCGA-LUSC/TCGA-LUSC_estimate_scores.txt", header = TRUE)
estimate$sample <- substring(estimate$ID, 1, 15)
write.csv(estimate, file = 'E:/项目一/data/fig2/estimate.csv')

#================================================
# 1. Relationship with ESTIMATE Scores
#================================================

# Merge risk score and ESTIMATE data
df <- merge(riskscore, estimate, by.y = 'sample')
df$Tumourpurity <- cos(0.6049872018 + 0.0001467884 * df$ESTIMATEScore)

# Perform Wilcoxon tests
wilcox.test(df$StromalScore ~ df$risk_level)
wilcox.test(df$ImmuneScore ~ df$risk_level)
wilcox.test(df$ESTIMATEScore ~ df$risk_level)
wilcox.test(df$Tumourpurity ~ df$risk_level)

# Scatter plot
p <- ggplot(df, aes(x = riskscore, y = Tumourpurity)) +
  geom_point(aes(color = risk_level, shape = risk_level), size = 1.8, alpha = 0.7) +
  scale_shape_manual(values = c(16, 18)) +
  geom_smooth(aes(color = risk_level), method = 'lm') +
  scale_color_manual(values = c("#C16E71", "#6e8fb2")) +
  scale_fill_manual(values = c("#C16E71", "#6e8fb2")) +
  ggpubr::stat_cor(aes(color = risk_level), label.x = 0.5, size = 5) +
  theme_bw(base_size = 15)

# Add marginal plots
ggMarginal(
  p,
  type = 'density',
  margins = 'both',
  size = 5,
  groupColour = TRUE,
  groupFill = TRUE
)

#================================================
# 2. Proportions of Signatures with Significant Differences
#================================================

# Load immune contexture data
immune.contexture <- read.csv('E:/项目一/data/fig2/immune contexture.csv')
immune.contexture.df <- fpkm[, colnames(fpkm) %in% immune.contexture$id]
immune.contexture.df$sample <- substring(fpkm$X, 1, 15)

# Merge with risk score data
df <- merge(immune.contexture.df, riskscore, by.y = 'sample')

# Perform Wilcoxon tests for each gene
a <- colnames(df)[2:99]
results_df <- data.frame(Gene = character(), P.Value = numeric(), stringsAsFactors = FALSE)
for (gene in a) {
  test_result <- wilcox.test(
    x = df[[gene]][df$risk_level == 'High'],
    y = df[[gene]][df$risk_level == 'Low'],
    alternative = "two.sided"
  )
  results_df <- rbind(results_df, data.frame(Gene = gene, P.Value = test_result$p.value))
}
results_df$adjusted_p_values <- p.adjust(results_df$P.Value, method = "fdr")

# Save significant genes
immune.contexture.sig.gene <- results_df[results_df$adjusted_p_values < 0.05, 1]
immune.contexture.sig.gene <- immune.contexture.sig.gene[c(-82)]
write.csv(results_df, file = 'E:/项目一/结果表格/immune.contexture.sig.gene.csv')


# Load immune gene list
immunegene_list1793 <- read.csv('E:/项目一/data/fig2/immunegene_list1793.csv')
immune.1793.df <- fpkm[, colnames(fpkm) %in% immunegene_list1793$id]
immune.1793.df$sample <- substring(fpkm$X, 1, 15)

# Merge with risk score data
df <- merge(immune.1793.df, riskscore, by.y = 'sample')

# Perform Wilcoxon tests for each gene
a <- colnames(df)[2:1315]
results_df <- data.frame(Gene = character(), P.Value = numeric(), stringsAsFactors = FALSE)
for (gene in a) {
  test_result <- wilcox.test(
    x = df[[gene]][df$risk_level == 'High'],
    y = df[[gene]][df$risk_level == 'Low'],
    alternative = "two.sided"
  )
  results_df <- rbind(results_df, data.frame(Gene = gene, P.Value = test_result$p.value))
}
results_df$adjusted_p_values <- p.adjust(results_df$P.Value, method = "fdr")

# Save significant genes
immune.1793.sig.gene <- results_df[results_df$adjusted_p_values < 0.05, 1]
immune.1793.sig.gene <- immune.1793.sig.gene[!is.na(immune.1793.sig.gene)]
write.csv(results_df, file = 'E:/项目一/结果表格/immune.1793.sig.gene.csv')


# Merge LUSC data with risk score data
df <- merge(lusc, riskscore, by = 'sample')

# Perform Wilcoxon tests for each gene
a <- colnames(df)[2:229]
results_df <- data.frame(Gene = character(), P.Value = numeric(), stringsAsFactors = FALSE)
for (gene in a) {
  test_result <- wilcox.test(
    x = df[[gene]][df$risk_level == 'High'],
    y = df[[gene]][df$risk_level == 'Low'],
    alternative = "two.sided"
  )
  results_df <- rbind(results_df, data.frame(Gene = gene, P.Value = test_result$p.value))
}
results_df$adjusted_p_values <- p.adjust(results_df$P.Value, method = "fdr")

# Save significant genes
lm.sig.gene <- results_df[results_df$adjusted_p_values < 0.05, 1]
write.csv(results_df, file = 'E:/项目一/结果表格/lm228.sig.gene.csv')

#================================================
# 3. Cell Heatmap
#================================================

# Load infiltration estimation data
infiltration_estimation <- openxlsx::read.xlsx('E:/项目一/data/fig2/infiltration_estimation_for_tcga(01样本).xlsx')

# Prepare data for heatmap
dfcell <- riskscore[, c('sample', 'riskscore', 'risk_level')]
dfcell$sample <- substring(dfcell$sample, 1, 12)
names(infiltration_estimation)[1] <- 'sample'
dfcell <- merge(dfcell, infiltration_estimation, by = 'sample')

# Perform Wilcoxon tests for each gene
df <- type.convert(dfcell)
df <- df[, c(1:9, 32:63, 65:121)]

a <- colnames(df)[4:98]
results_df <- data.frame(Gene = character(), P.Value = numeric(), stringsAsFactors = FALSE)
for (gene in a) {
  test_result <- wilcox.test(
    x = df[[gene]][df$risk_level == 'High'],
    y = df[[gene]][df$risk_level == 'Low'],
    alternative = "two.sided"
  )
  results_df <- rbind(results_df, data.frame(Gene = gene, P.Value = test_result$p.value))
}
results_df$adjusted_p_values <- p.adjust(results_df$P.Value, method = "fdr")

results_df$significance <- ifelse(results_df$adjusted_p_values < 0.001, "***",
                                  ifelse(results_df$adjusted_p_values < 0.01, "**",
                                         ifelse(results_df$adjusted_p_values < 0.05, "*", "NS")))

results_df=results_df[results_df$adjusted_p_values<0.05,]
rowanno=data.frame(sig=results_df$significance,method=results_df$method,row.names = results_df$Gene)


dfrt=arrange(df,riskscore)
dfrt=distinct(dfrt,sample,.keep_all = T)
rt=dfrt[,c(results_df$Gene)]
rownames(rt)=dfrt$sample
colnames(rt)=results_df$method1
group=data.frame(group=dfrt$risk_level,riskscore=dfrt$riskscore,row.names = dfrt$sample)
bk <- c(seq(-2,-0.1,by=0.01),seq(0,2,by=0.01))


ann_colors = list(
  group = c('High'="#C16E71",'Low'="#6e8fb2"),
  method=c("TIMER"='#CBE6E0',"CIBERSORT-ABS"='#FEF7CB',"QUANTISEQ"='#C1BED7',
           "MCPCOUNTER"='#8AB0D1' ,"XCELL"=  '#F2C896', "EPIC"='#F6E1D4' ),
  sig= c('*'= 'grey', '**'= "darkgrey","***"='black')
)

pheatmap(t(rt),cluster_rows =F,cluster_cols =F,show_colnames =F,show_rownames = T,
         fontfamily = "serif",
         color = c(colorRampPalette(colors = c("#248D89", "white"))(length(bk)/2),
                   colorRampPalette(colors = c("white", "#C26C35"))(length(bk)/2)),
         breaks=bk,
         annotation_colors = ann_colors,
         annotation_col =group,
         annotation_row =rowanno,
         scale = "row")


#============================================
# 4. Procrustes Analysis of 7-step Genes
#============================================

# Load 7-step gene data
folder_path <- "E:/STAD/乱七八糟/BLCA/7step"  # Replace with your folder path
csv_files <- list.files(path = folder_path, pattern = "\\.csv$", full.names = TRUE)
csv_data_list <- map(csv_files, read.csv, stringsAsFactors = FALSE)

# Extract all 7-step genes
all_7stepgenes <- c()
for (df in csv_data_list) {
  all_7stepgenes <- c(all_7stepgenes, df[["id"]])
}

# Load FPKM data and filter for 7-step genes
fpkm <- read.csv('E:/TCGA/TCGA_Annovar/TCGA-LUSC/TCGA-LUSC_transcriptome_profilingonts_mRNA_fpkm_annovar_format.csv')
df1 <- fpkm[, colnames(fpkm) %in% all_7stepgenes]
df1$sample <- substring(fpkm$X, 1, 15)

# Load LUSC data
lusc <- read.csv('E:/项目一/data/fig2/lusc.csv', header = TRUE) %>%
  column_to_rownames(var = 'X')

# Merge FPKM and LUSC data
df <- merge(df1, lusc, by = 'sample')
df <- distinct(df, sample, .keep_all = TRUE)

# Log-transform gene expression data
df[, 2:403] <- log(df[, 2:403] + 1)

# Extract 7-step genes and lactate metabolism genes
df7step <- df[, colnames(df) %in% all_7stepgenes]
dflm <- df[, colnames(df) %in% colnames(lusc)[-c(229)]]
rownames(df7step) <- df$sample
rownames(dflm) <- df$sample

# Calculate distance matrices
df7step.dist <- vegdist(df7step)
dflm.dist <- vegdist(dflm)

# Perform Mantel test
mantel_result <- mantel(df7step.dist, dflm.dist)
print(mantel_result)

# Perform Procrustes analysis
mds.s <- monoMDS(df7step.dist)
mds.e <- monoMDS(dflm.dist)
pro.s.e <- procrustes(mds.s, mds.e, symmetric = TRUE)
summary(pro.s.e)

# Perform Procrustes test with permutations
set.seed(1)
pro.s.e_t <- protest(mds.s, mds.e, permutations = 999)
print(pro.s.e_t)

# Plot Procrustes analysis results
Pro_Y <- cbind(data.frame(pro.s.e$Yrot), data.frame(pro.s.e$X))
Pro_X <- data.frame(pro.s.e$rotation)

ggplot(data = Pro_Y) +
  geom_segment(aes(x = X1, y = X2, xend = (X1 + MDS1) / 2, yend = (X2 + MDS2) / 2),
               arrow = arrow(length = unit(0, 'cm')), color = "#9BBB59", size = 0.5) +
  geom_segment(aes(x = (X1 + MDS1) / 2, y = (X2 + MDS2) / 2, xend = MDS1, yend = MDS2),
               color = "#957DB1", size = 0.5) +
  geom_point(aes(X1, X2), color = "#9BBB59", size = 2, shape = 16) +
  geom_point(aes(MDS1, MDS2), color = "#957DB1", size = 2, shape = 16) +
  theme(panel.grid = element_blank(),
        panel.background = element_rect(color = 'black', fill = 'transparent'),
        legend.key = element_rect(fill = 'transparent'),
        axis.ticks.length = unit(0.4, "lines"),
        axis.ticks = element_line(color = 'black'),
        axis.line = element_line(colour = "black"),
        axis.title.x = element_text(colour = 'black', size = 16),
        axis.title.y = element_text(colour = 'black', size = 16),
        axis.text = element_text(colour = 'black', size = 14)) +
  labs(x = 'Dimension 1', y = 'Dimension 2', color = '') +
  labs(title = "Correlation between 7-step genes and LM genes") +
  geom_vline(xintercept = 0, color = 'gray', linetype = 2, size = 0.3) +
  geom_hline(yintercept = 0, color = 'gray', linetype = 2, size = 0.3) +
  annotate('text', label = 'Mantel statistic:\nr = 0.4753, p-value = 0.001',
           x = -0.09, y = 0.05, size = 5, hjust = 0) +
  annotate('text', label = 'Procrustes analysis:\nM2 = 0.6297, p-value = 0.001',
           x = 0, y = -0.06, size = 5, hjust = 0) +
  theme(plot.title = element_text(size = 17, colour = "black", hjust = 0.8, face = "bold"))

# Perform Mantel test for risk score and 7-step genes
df <- merge(df, riskscore[, c('sample', 'riskscore')], by.y = 'sample')
rs <- data.frame(riskscore = df$riskscore, row.names = df$sample)
rs.dist <- vegdist(rs)
dfstep2 <- df[, colnames(df) %in% csv_data_list[[7]]$id]
rownames(dfstep2) <- df$sample
dfstep2.dist <- vegdist(dfstep2)  # Default method is Bray-Curtis
mantel_result <- mantel(rs.dist, dfstep2.dist)
print(mantel_result)

# Plot Mantel test results
coef <- data.frame(step = paste0("Step", 1:7), r = c(0.2683, 0.1749, 0.119, 0.1827, 0.465, 0.3746, 0.1106))
ggplot(coef, aes(x = step, y = r, group = NA)) +
  geom_line(size = 1, linetype = 'dashed') +
  geom_point(size = 6, shape = 21, colour = 'darkred', fill = 'pink') +
  labs(x = "Antitumor immune cycle", y = "Value") +
  theme_bw(base_size = 15) +
  annotate('text', label = '(Mantel test r; P value)',
           x = 1.5, y = 0.4, size = 5, hjust = 0) +
  theme(panel.grid.minor = element_blank())

#===============================================
# 5. Chord Diagram of Immune Scores and 7 Genes
#===============================================
library(circlize)
library(corrplot)
library(Hmisc)

# Merge risk score and ESTIMATE data
df <- merge(riskscore, estimate, by.y = 'sample')
df$Tumourpurity <- cos(0.6049872018 + 0.0001467884 * df$ESTIMATEScore)

# Perform correlation analysis
library(psych)  # For correlation analysis
library(circlize)  # For circular visualization

df <- type.convert(df)
gene <- df[, 4:10]
score <- df[, c('StromalScore', 'ImmuneScore', 'Tumourpurity')]

corr <- corr.test(gene, score, method = "spearman", adjust = "BH")
r <- corr$r  # Correlation matrix
r <- na.omit(r)
p <- corr$p
p <- na.omit(p)

gene_cor <- reshape2::melt(r)
gene_cor <- subset(gene_cor, value != 0)  

col = RColorBrewer::brewer.pal(n = 10,name = 'Set3')
names(col) = colnames(df)[c(4:10,15,16,18)]

chordDiagram(gene_cor,annotationTrack = c('grid', 'name', 'axis'),  
             grid.col = col,  
             col = colorRamp2(c(-0.5, 0, 0.5), c("#9CC5DA", 'white', "#D0908F"), transparency =0.5),
             annotationTrackHeight = c(0.05, 0.08), 
)


#========================================================================
# 6. Analysis of ICI Response Correlations with LM-related Genes and Risk Score
#========================================================================

#-----------------------------
# 0. Environment Setup
#-----------------------------
library(tibble)
library(dplyr)
library(data.table)
library(GEOquery)
library(pheatmap)
library(RColorBrewer)

#-----------------------------
# 6.1. Data Loading and Preprocessing
#-----------------------------
# Load count data
GSE126044_counts <- read.table('data/fig3/immunotherapy_cohort/GSE126044_counts.txt', row.names = NULL)

# Load clinical data
GSE126044_clinical <- read.xlsx('data/fig3/immunotherapy_cohort/GSE126044_clinical.xlsx')

# Load gene length data for TPM calculation
gene_length <- fread('data/fig3/hg19gene.txt', data.table = FALSE)
colnames(gene_length) <- c("gene_name", "Length")

# Filter genes of interest
GSE126044_counts <- GSE126044_counts %>%
  filter(row.names %in% gene) 
names(GSE126044_counts)[1] <- 'gene_name'

# Merge count data with gene length
count_length <- inner_join(gene_length, GSE126044_counts)

#-----------------------------
# 6.2. TPM Normalization
#-----------------------------
tpm_value <- count_length
for (i in 3:ncol(count_length)) {
  result <- round((count_length[, i] * 1000 * 1000000) / 
                    (count_length[, 2] * sum((count_length[, i] * 1000 / count_length[, 2]))), 3)  
  tpm_value[, i] <- result
}

# Prepare normalized expression matrix
df <- tpm_value[,3:18]
rownames(df) <- tpm_value$gene_name
df <- t(df) %>% as.data.frame()
df$id <- rownames(df)
df <- merge(df, GSE126044_clinical)
dftpm <- df

#-----------------------------
# 6.3. Risk Score Calculation
#-----------------------------
# Load coefficients from Cox model
dfcoef <- read.csv('results/fig1/multicox.csv')
coef <- dfcoef$estimate
names(coef) <- dfcoef$term

# Calculate risk score
dftpm <- dftpm %>%
  mutate(riskscore = CHEK2*coef[1] + LIPT1*coef[2] + TUFM*coef[3] + 
           NDUFA10*coef[4] + AGK*coef[5] + PNPLA2*coef[6] + GFM1*coef[7]) %>%
  mutate(risklevel = ifelse(riskscore > median(riskscore), 'high', 'low'))
names(dftpm)[228] <- 'response'
df <- dftpm

#-----------------------------
# 6.4 t-SNE Analysis
#-----------------------------
# Calculate distance matrix
dist <- vegdist(dftpm[,2:227], method = "bray", na.rm = TRUE)

# PERMANOVA test
Adonis <- adonis2(dist ~ pcoa_points$group,
                  distance = "bray",
                  permutations = 999)

# t-SNE dimensionality reduction
res <- Rtsne(dftpm[,2:227],
             dims = 2,
             max_iter = 1000,
             theta = 0.5,
             perplexity = 5,
             verbose = TRUE)

# Prepare plotting data
result <- as.data.frame(res$Y)
colnames(result) <- c("tSNE1", "tSNE2")
tsne_points <- data.frame(result, group = dftpm$response)

# Plot t-SNE results
color_palette <- c("#1597A5", "#6B5B95")
p1 <- ggplot(data = tsne_points, aes(x = tSNE1, y = tSNE2)) +
  theme_bw(base_size = 15) +
  geom_point(aes(color = group), shape = 19, size = 3, alpha = 0.8) +
  theme(panel.grid = element_blank()) +
  geom_vline(xintercept = 0, lty = "dashed", size = 1, color = 'grey50') +
  geom_hline(yintercept = 0, lty = "dashed", size = 1, color = 'grey50') +
  stat_ellipse(data = tsne_points,
               geom = "polygon", level = 0.95,
               linetype = 2, size = 0.5,
               aes(fill = group),
               alpha = 0.2) +
  scale_color_manual(values = color_palette) +
  scale_fill_manual(values = color_palette) +
  theme(panel.grid = element_blank())

# Add statistical annotation
p1 + annotate('text', 
              label = 'Permanova test:\np-value = 0.022',
              x = -150, y = 280, size = 5, hjust = 0)

#-----------------------------
# 6.5. Boxplot Visualization
#-----------------------------
library(rstatix)
col <- c("#1597A5", "#6B5B95")

# Create base boxplot
p <- ggplot(data = df, aes(x = response, y = riskscore, colour = response)) +
  geom_boxplot(aes(fill = response), alpha = 0.2, size = 0.6, width = 0.7) +
  geom_jitter(alpha = 0.3, size = 2, width = 0.1) +
  theme_bw(base_size = 15) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  scale_color_manual(values = col) +
  scale_fill_manual(values = col) +
  labs(x = "", y = 'Risk Score')

# Add statistical significance
stat.test <- df %>%
  wilcox_test(riskscore ~ response) %>%
  add_significance() %>%
  add_y_position()

# Final plot with p-value
p + stat_pvalue_manual(stat.test, label = "p = {p}", tip.length = 0.01) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)))

#-----------------------------
# 6.6. Gene Expression Heatmap
#-----------------------------
# Prepare data for heatmap
dfrt <- arrange(dftpm, response)
rt <- dfrt[, c(dfcoef$term)]
rownames(rt) <- dfrt$id

# Create risk score categories
breaks <- seq(-366, 80, by = 70)
riskscore_categories <- cut(dfrt$riskscore, breaks = breaks)
df$riskscore_cat <- as.factor(riskscore_categories)

# Prepare annotation data
group <- data.frame(group = dfrt$response,
                    riskscore = riskscore_categories,
                    row.names = dfrt$id)

# Color settings
riskscore_colors <- colorRampPalette(colors = c("#FFFFFF", "#000000"))(length(levels(df$riskscore_cat)))
bk <- c(seq(-2, -0.1, by = 0.01), seq(0, 2, by = 0.01))

ann_colors <- list(
  group = c('non-responder' = "#6FBFC8", 'responder' = "#A39ABE"),
  riskscore = setNames(riskscore_colors, levels(df$riskscore_cat))
)

# Generate heatmap
pheatmap(t(rt),
         cluster_rows = TRUE,
         cluster_cols = FALSE,
         show_colnames = TRUE,
         show_rownames = TRUE,
         fontfamily = "serif",
         color = c(colorRampPalette(colors = c("#248D89", "white"))(length(bk)/2),
                   colorRampPalette(colors = c("white", "#C26C35"))(length(bk)/2)),
         legend_breaks = seq(-2, 2, 2),
         breaks = bk,
         annotation_colors = ann_colors,
         annotation_col = group,
         scale = "row")
