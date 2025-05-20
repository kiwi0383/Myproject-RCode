#========================================================================
# Figure 4. Different driver genes and biological pathways in LM.Sig-based subgroups
#========================================================================

#-----------------------------
# 0. Environment Setup
#-----------------------------
library(maftools)
library(data.table)
library(openxlsx)
library(dplyr)
library(rstatix)
library(TCGAmutations)
library(tidyverse)
library(survival)
library(survminer)
library(clusterProfiler)
library(GSVA)
library(msigdbr)

#-----------------------------
# 1. Data Loading and Preparation
#-----------------------------
# Load TCGA LUSC MAF data
maf_data <- tcga_load(study = "LUSC")

# Prepare sample groups
riskscore <- read.xlsx("data/fig1/df.xlsx") %>%
  select(sample, risk_level, riskscore)

all_sample_name <- riskscore$sample
high_risk_samples <- riskscore[riskscore$risk_level == 'High', "sample"]
low_risk_samples <- riskscore[riskscore$risk_level == 'Low', "sample"]

# Match TCGA barcodes with our samples
maf_sample <- data.frame(
  fullid = maf_data@clinical.data[["Tumor_Sample_Barcode"]],
  idtcga = substring(maf_data@clinical.data[["Tumor_Sample_Barcode"]], 1, 15)
)

# Get overlapping samples
all_fullid <- maf_sample[maf_sample$idtcga %in% all_sample_name, 1]
high_fullid <- maf_sample[maf_sample$idtcga %in% high_risk_samples, 1]
low_fullid <- maf_sample[maf_sample$idtcga %in% low_risk_samples, 1]

# Create MAF subsets
maf_high <- subsetMaf(maf = maf_data, tsb = high_fullid)
maf_low <- subsetMaf(maf = maf_data, tsb = low_fullid)

#-----------------------------
# 2. Mutation Landscape Overview
#-----------------------------
# Summary plot
plotmafSummary(maf = maf_data, 
               rmOutlier = FALSE,
               fs = 1.1,
               titleSize = c(1.4, 1.2),
               addStat = 'median', 
               dashboard = TRUE, 
               titvRaw = FALSE)

# Define variant classification colors
vc_cols <- RColorBrewer::brewer.pal(n = 8, name = 'Paired')
names(vc_cols) <- c(
  'Frame_Shift_Del',
  'Missense_Mutation',
  'Nonsense_Mutation',
  'Multi_Hit',
  'Frame_Shift_Ins',
  'In_Frame_Ins',
  'Splice_Site',
  'In_Frame_Del'
)

# Oncoplots for high and low risk groups
oncoplot(maf = maf_high, 
         colors = vc_cols, 
         top = 20,
         draw_titv = TRUE,
         fontSize = 0.9,
         legendFontSize = 1.5)

oncoplot(maf = maf_low, 
         colors = vc_cols, 
         top = 20,
         fontSize = 0.9,
         legendFontSize = 1.5)

#-----------------------------
# 3. Mutation Interactions
#-----------------------------
# Co-occurrence and mutual exclusivity analysis
somaticInteractions(maf_high, 
                    fontSize = 0.9, 
                    showCounts = FALSE,
                    top = 15,
                    showSum = FALSE,
                    pvalue = c(0.05, 0.1))

somaticInteractions(maf_low, 
                    fontSize = 0.9,
                    top = 15, 
                    showSum = FALSE,
                    pvalue = c(0.05, 0.1))

#-----------------------------
# 4. Tumor Mutational Burden (TMB) Analysis
#-----------------------------
# Calculate TMB
tmb_high <- tmb(maf_high, captureSize = 38, logScale = TRUE)
tmb_low <- tmb(maf_low, captureSize = 38, logScale = TRUE)

tmb_high$group <- "High"
tmb_low$group <- "Low"

tmb_df <- rbind(tmb_high, tmb_low)

# Statistical test
wilcox.test(total_perMB_log ~ group, data = tmb_df)

# Visualization
tmb_plot <- ggplot(data = tmb_df, 
                   aes(x = group, y = total_perMB_log, colour = group)) +
  geom_boxplot(aes(fill = group), alpha = 0.3, size = 0.7, width = 0.5) + 
  geom_jitter(alpha = 0.3, size = 2, width = 0.1) +
  theme_bw(base_size = 15) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  scale_color_manual(values = c("#7B87B9", "#CC9982")) +
  scale_fill_manual(values = c("#7B87B9", "#CC9982")) +
  labs(x = "", y = "TMB (log10)")

# Add p-value
stat_test <- tmb_df %>% 
  wilcox_test(total_perMB_log ~ group) %>% 
  add_significance() %>% 
  add_y_position()

tmb_plot + 
  stat_pvalue_manual(stat_test, label = "p = {p}", tip.length = 0.01) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)))

#-----------------------------
# 5. Survival Analysis
#-----------------------------
# Gene-level survival analysis
mafSurvival(maf = maf_high,
            genes = 'MUC16',
            time = 'CDR_OS.time',
            Status = 'CDR_OS')

# Batch survival analysis for top mutated genes
survGroup(maf = maf_low,
          top = 50, 
          geneSetSize = 1, 
          time = 'CDR_OS.time',
          Status = 'CDR_OS',
          verbose = FALSE)

df1=subsetMaf(maf_high,genes = c("ZNF804B"))
mutant=df1@clinical.data$Tumor_Sample_Barcode_min
km1=data.frame(id=maf_high@clinical.data$Tumor_Sample_Barcode_min,
               OS=maf_high@clinical.data$CDR_OS,
               OS.time=maf_high@clinical.data$CDR_OS.time)
km1$group=ifelse(km1$id %in% mutant,"high_mutant","high_wt")

df2=subsetMaf(maf_low,genes = c("ZNF804B"))
mutant=df2@clinical.data$Tumor_Sample_Barcode_min
km2=data.frame(id=maf_low@clinical.data$Tumor_Sample_Barcode_min,
               OS=maf_low@clinical.data$CDR_OS,
               OS.time=maf_low@clinical.data$CDR_OS.time)
km2$group=ifelse(km2$id %in% mutant,"low_mutant","low_wt")

km_PCLO=rbind(km1,km2)

fit <- survfit(Surv(OS.time, OS) ~group, data = km_PCLO)
ggsurvplot(fit, data = km_PCLO, pval = T)

ggsurvplot(
  fit,
  data = km_PCLO,
  censor.shape="|", censor.size = 4,
  conf.int = F,
  pval = TRUE,
  palette = "npg",
  surv.median.line = "hv",
  ggtheme = theme_survminer(base_size = 20,base_line_size = 0.4),
  legend = "right",
  xlab = "OS_time(days)",
  ylab = "Survival probablity",
  title = "ZNF804B",
  break.x.by = 1000,
  break.y.by = 0.2)

#-----------------------------
# 6. Gene Set Enrichment Analysis (GSEA)
#-----------------------------
# Load expression data
fpkm <- read.csv('data/TCGA-LUSC_fpkm.csv')
fpkm$sample <- substring(fpkm$X, 1, 15)
expr_data <- merge(fpkm, riskscore[, c('sample', 'risk_level', 'riskscore')], 
                   by = 'sample')

# Differential expression analysis
diff_res <- diff_analysis(exprset = t(expr_data[, 3:19533]), 
                          group = expr_data$risk_level, 
                          is_count = FALSE)

# Prepare gene list for GSEA
gene_df <- data.frame(
  logFC = diff_res$deg_limma$logFC,
  SYMBOL = rownames(diff_res$deg_limma)
)

gene_df <- merge(gene_df, 
                 bitr(gene_df$SYMBOL, 
                      fromType = "SYMBOL", 
                      toType = "ENTREZID", 
                      OrgDb = "org.Hs.eg.db"),
                 by = "SYMBOL")

gene_list <- sort(gene_df$logFC, decreasing = TRUE)
names(gene_list) <- gene_df$SYMBOL

# Load Hallmark gene sets
hallmarks <- read.gmt("data/h.all.v2024.1.Hs.symbols.gmt")

# Run GSEA
gsea_res <- GSEA(geneList = gene_list, 
                 TERM2GENE = hallmarks)

# Visualization
dotplot(gsea_res, 
        showCategory = 20,
        split = ".sign",
        font.size = 10,
        label_format = 20) +
  facet_grid(~.sign) +
  scale_y_discrete(labels = function(x) substr(x, 10, nchar(x)))


geneSetID = c('HALLMARK_G2M_CHECKPOINT',
              'HALLMARK_E2F_TARGETS', 'HALLMARK_MYC_TARGETS_V1',"HALLMARK_MYC_TARGETS_V2",
              'HALLMARK_DNA_REPAIR','HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION',
              'HALLMARK_TNFA_SIGNALING_VIA_NFKB',
              'HALLMARK_IL6_JAK_STAT3_SIGNALING')

gseaNb(object = y, geneSetID = geneSetID, subPlot = 2,addPval = T,
       pvalX = 0.02, pvalY = -0.05,rmHt = T,termWidth = 20, base_size = 12,
       legend.position = c(0.8,0.8),curveCol = c("#009D73","#5BB3E4","#E59F24","#000000",
                                                 "#B19CD9","#E15759","#00796B","darkgrey") )

labels=rt[order(rt$FDR,decreasing =T),"Term"]
rt$Term = factor(rt$Term,levels=labels)
yd$RichFactor=as.numeric(substring(yd$leading_edge,6,7))/100
yd<-arrange(yd,RichFactor)
yd$adp<-log10(yd$p.adjust)

yd$ID=gsub("^HALLMARK_", "", rownames(yd))
ggplot(yd,aes(NES,reorder(ID,NES), fill=adp))+
  geom_bar(stat='identity')+
  #scale_fill_viridis(begin = 0, end = 1, option = "D")+
  scale_fill_distiller(palette = "YlOrRd",direction = 1)+
  xlab("NES") + ylab("Pathways") + 
  theme(axis.text.x=element_text(color="black", size=10),
        axis.text.y=element_text(color="black", size=10)) + 
  theme_classic()+scale_y_discrete(labels=function(x) stri_sub(x,10))+
  theme(legend.position = "right") +  
  theme(axis.text = element_text(size = 10)) +
  theme(axis.title.x = element_text(vjust = 0, size =10, face = "plain")) + 
  theme(axis.title.y = element_text(vjust = 0, size = 10, face = "plain")) +   
  theme(axis.line = element_line(color = "black",linewidth = 0.8)) 

ggplot(yd,aes(x=reorder(ID,NES),y=NES,fill=-log10(pvalue)))+
  geom_bar(stat='identity',alpha=0.6)+
  coord_flip()+
  scale_fill_distiller(palette = "GnBu",direction =1)+
  xlab("Pathways") + ylab("NES") + 
  geom_text(aes(y=0,label=ID,hjust=as.numeric(NES>0)),size=4)+
  theme_bw(base_size = 15)+
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        panel.grid.major.y=element_blank())

#-----------------------------
# 7. Single Sample GSEA (ssGSEA)
#-----------------------------
# Prepare expression matrix
expr_matrix <- t(expr_data[, 3:19532])
rownames(expr_matrix) <- expr_data$sample

# Get Hallmark gene sets
h_sets <- msigdbr(species = "Homo sapiens", category = "H") %>%
  split(x = .$gene_symbol, f = .$gs_name)

# Run ssGSEA
gsva_res <- gsva(expr = as.matrix(expr_matrix),
                 gset.idx.list = h_sets,
                 method = "ssgsea")

# Compare pathway activity between groups
pathway_results <- data.frame(
  pathway = character(),
  P.Value = numeric(),
  stringsAsFactors = FALSE
)

for (path in colnames(gsva_res)) {
  test_res <- wilcox.test(
    x = gsva_res[path, expr_data$risk_level == 'High'],
    y = gsva_res[path, expr_data$risk_level == 'Low'],
    alternative = "two.sided"
  )
  pathway_results <- rbind(pathway_results, 
                           data.frame(pathway = path, 
                                      P.Value = test_res$p.value))
}

pathway_results$adj.P.Val <- p.adjust(pathway_results$P.Value, method = "fdr")
sig_pathways <- pathway_results[pathway_results$adj.P.Val < 0.05, ]

#visualization
results_df$adjusted_p_values <- p.adjust(results_df$P.Value, method = "fdr")
result<-results_df[results_df$adjusted_p_values < 0.05, ]

hallmarks = c('G2M_CHECKPOINT',
              'E2F_TARGETS', 'MYC_TARGETS_V1',"MYC_TARGETS_V2",
              'DNA_REPAIR','EPITHELIAL_MESENCHYMAL_TRANSITION',
              'TNFA_SIGNALING_VIA_NFKB',
              'IL6_JAK_STAT3_SIGNALING') 


plots <- list()

for (hallmark in hallmarks) {
  p <- ggplot(data = gsva_data, aes_string(x = "risk_level", y = hallmark, colour = "risk_level")) +
    geom_boxplot(aes_string(fill = "risk_level"), alpha = 0.1, size = 0.6, width = 0.7) +
    geom_jitter(alpha = 0.3, size = 2, width = 0.1) +
    theme_bw(base_size = 15) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          legend.position ='none'
          #       plot.title = element_text(size = 8)
    ) +
    scale_color_manual(values = col) +
    scale_fill_manual(values = col) +
    labs(x = "", y =hallmark)
  #+ggtitle(hallmark)
  
  stat.test <- gsva_data %>%
    wilcox_test(as.formula(paste(hallmark, "~ risk_level"))) %>%
    add_significance() %>%
    add_y_position()
  
  p <- p + stat_pvalue_manual(stat.test, label = "p = {p}", tip.length = 0.01) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)))
  
  plots[[hallmark]] <- p
}

library(ggpubr)
ggarrange(plotlist = plots[c(hallmarks)], ncol = 4,nrow=2)