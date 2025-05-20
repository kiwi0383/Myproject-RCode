#================================================
# Figure 1. Development of prognostic model based on the LM-related genes.
#================================================

# Load necessary libraries
library(openxlsx)
library(easyTCGA)
library(dplyr)
library(data.table)
library(R.utils)
library(survival)
library(survminer)
library(openxlsx)
library(My.stepwise)
library(Hmisc)
library(glmnet)
library(randomForestSRC)
library(GOSemSim)
library(ggplot2)
library(pheatmap)
library(AnnoProbe)
library(GEOquery)
library(ggsignif)
library(corrplot)
library(pheatmap)


# Set working directory (adjust as needed)
setwd("E:/项目一")

#================================================
# 1.Data Preparation
#================================================

# Load data
fpkm <- read.csv('E:/TCGA/TCGA_Annovar/TCGA-LUSC/TCGA-LUSC_transcriptome_profilingonts_mRNA_fpkm_annovar_format.csv')
lmgene <- read.xlsx('E:/项目一/LM.gene.xlsx')
gene <- lmgene$lmgene

# Filter FPKM data for genes of interest
lusc <- fpkm[, colnames(fpkm) %in% gene]
rownames(lusc) <- fpkm$X

# Save filtered data
write.csv(lusc, file = 'E:/项目一/data/fig2/lusc.csv')

#================================================
# 2.Survival Analysis and Cox Regression
#================================================

# Load survival data
a <- tcga_cdr_surv %>%
  filter(`cancer type abbreviation` == 'LUSC')

# Prepare data for survival analysis
lusc$sample <- substring(rownames(lusc), 1, 15)
lusc <- distinct(lusc, sample, .keep_all = TRUE)
osdata <- merge(a[, c(1, 26, 27)], lusc, by.y = 'sample')

# Single-variable Cox regression
covariates <- colnames(osdata)[-c(1, 2, 3)]
results <- data.frame(covariate = character(0),
                      beta = numeric(0),
                      HR.confint.lower = numeric(0),
                      HR.confint.upper = numeric(0),
                      wald.test = numeric(0),
                      p.value = numeric(0))

for (covariate in covariates) {
  formula <- as.formula(paste("Surv(OS.time, OS) ~", covariate))
  model <- coxph(formula, data = osdata)
  summary <- summary(model)
  p.value <- signif(summary$wald["pvalue"], digits = 2)
  wald.test <- signif(summary$wald["test"], digits = 2)
  beta <- signif(summary$coef[1], digits = 2)
  HR.confint.lower <- signif(summary$conf.int["lower .95"], 2)
  HR.confint.upper <- signif(summary$conf.int["upper .95"], 2)
  results <- rbind(results, data.frame(covariate = covariate,
                                       beta = beta,
                                       HR.confint.lower = HR.confint.lower,
                                       HR.confint.upper = HR.confint.upper,
                                       wald.test = wald.test,
                                       p.value = p.value))
}

# Save single-variable Cox regression results
write.csv(results, file = 'E:/项目一/结果表格/单变量cox.csv')
cox.sig.gene <- results[results$p.value < 0.05, 1]

# Lasso regression
osdata <- osdata %>% filter(OS.time > 0)
x <- as.matrix(osdata[, cox.sig.gene])
y <- Surv(osdata$OS.time, osdata$OS)

set.seed(1234)
alpha1_fit <- glmnet(x, y, alpha = 1, family = "cox", nlambda = 100)
plot(alpha1_fit, xvar = "lambda", label = TRUE)

# Cross-validation for Lasso
cvfit <- cv.glmnet(x, y, family = "cox", alpha = 1)
plot(cvfit, xvar = "lambda", label = TRUE)
cvfit$lambda.min

# Extract Lasso results
coef.min <- coef(cvfit, s = "lambda.min")
active.min <- which(coef.min != 0)
geneids <- rownames(coef.min)[active.min]

lasso.result <- data.frame(coef = coef.min@x)
lasso.result$gene <- rownames(lasso.result)
write.csv(lasso.result, file = 'E:/项目一/结果表格/lasso_coef15gene.csv')

#================================================
# 3.Random Survival Forest (RSF) and Gene Friends
#================================================

# RSF analysis
df <- osdata[, c('OS', 'OS.time', geneids)]
fit.RSF <- rfsrc(Surv(OS.time, OS) ~ ., data = df, ntree = 1000, nodesize = 10,
                 splitrule = 'logrank', importance = TRUE, proximity = TRUE, forest = TRUE, seed = 1234)

# Plot RSF results
plot(fit.RSF)

# Extract variable importance
importance_gene <- data.frame(fit.RSF$importance) %>%
  rownames_to_column("gene") %>%
  arrange(-fit.RSF.importance) %>%
  head(10)

# Save RSF results
write.csv(importance_gene, file = 'E:/项目一/结果表格/rsf_importance.csv')

# Gene Friends analysis
rt <- bitr(gene, fromType = 'SYMBOL', toType = 'ENTREZID', OrgDb = 'org.Hs.eg.db')
rt$ENTREZID <- as.character(rt$ENTREZID)

bp <- godata('org.Hs.eg.db', ont = "BP", computeIC = FALSE)
cc <- godata('org.Hs.eg.db', ont = "CC", computeIC = FALSE)
mf <- godata('org.Hs.eg.db', ont = "MF", computeIC = FALSE)

simbp <- mgeneSim(rt$ENTREZID, semData = bp, measure = "Wang", drop = NULL, combine = "BMA")
simcc <- mgeneSim(rt$ENTREZID, semData = cc, measure = "Wang", drop = NULL, combine = "BMA")
simmf <- mgeneSim(rt$ENTREZID, semData = mf, measure = "Wang", drop = NULL, combine = "BMA")
fsim <- (simmf * simcc * simbp)^(1/3)
colnames(fsim) <- rt$SYMBOL
rownames(fsim) <- rt$SYMBOL

# Extract top gene friends
gene.friends <- names(sort(rowMeans(fsim), decreasing = TRUE)[1:10])
intersect.gene.new <- intersect(importance_gene$gene, gene.friends)

#================================================
# 4.Multivariate Cox Regression
#================================================

# Prepare data for multivariate Cox regression
df <- osdata[, c('OS', 'OS.time', intersect.gene.new)]
model <- coxph(Surv(OS.time, OS) ~ ., data = df)
summary(model)

# Save multivariate Cox regression results
multicox.result <- broom::tidy(model)
write.csv(multicox.result, file = 'E:/项目一/结果表格/multicox.csv')

# Calculate risk score
df$riskscore <- predict(model, newdata = df, type = 'risk')
df$risk_level <- as.vector(ifelse(df$riskscore > median(df$riskscore), "High", "Low"))
write.csv(df, file = 'E:/项目一/结果表格/tcga-riskscore.csv')

# Survival analysis based on risk score
fit <- survfit(Surv(OS.time, OS) ~ risk_level, data = df)
km1 <- ggsurvplot(fit, data = df, censor.shape = "|", censor.size = 4, conf.int = TRUE,
                  conf.int.style = "ribbon", conf.int.alpha = 0.2, pval = TRUE,
                  ggtheme = theme_survminer(), palette = c("#b76d7b", "#87B5C1"),
                  surv.median.line = "hv", legend = "top", xlab = "OS_time(days)",
                  ylab = "Survival probability", title = "TCGA", break.x.by = 1000, break.y.by = 0.2)



#================================================
# 5.Merge Clinical Data with Risk Score Data
#================================================

# Load clinical data
clinical <- fread(file = "E:/项目一/TCGA-LUSC.clinical.tsv.gz", header = TRUE)

# Process clinical data
clinical <- clinical %>%
  mutate(
    stage = case_when(
      ajcc_pathologic_stage.diagnoses %in% c('Stage I', 'Stage IA', 'Stage IB') ~ 'Stage I',
      ajcc_pathologic_stage.diagnoses %in% c('Stage II', 'Stage IIA', 'Stage IIB') ~ 'Stage II',
      ajcc_pathologic_stage.diagnoses %in% c('Stage III', 'Stage IIIA', 'Stage IIIB', 'Stage IIIC') ~ 'Stage III',
      ajcc_pathologic_stage.diagnoses == 'Stage IV' ~ 'Stage IV'
    ),
    tstage = case_when(
      ajcc_pathologic_t.diagnoses %in% c('T1', 'T1a', 'T1b') ~ 'T1',
      ajcc_pathologic_t.diagnoses %in% c('T2', 'T2a', 'T2b') ~ 'T2',
      TRUE ~ ajcc_pathologic_t.diagnoses
    ),
    nstage = ajcc_pathologic_n.diagnoses,
    mstage = case_when(
      ajcc_pathologic_m.diagnoses %in% c('M1', 'M1a', 'M1b') ~ 'M1',
      TRUE ~ ajcc_pathologic_m.diagnoses
    ),
    age = ifelse(age_at_index.demographic >= 65, '>=65', '<65'),
    sample = substring(sample, 1, 15)
  )

# Load risk score data (assuming df is already defined and contains risk score)
df_clinical <- df
df_clinical <- merge(df_clinical, clinical[, c(1, 91:95)], by = 'sample')
df_clinical$status <- ifelse(df_clinical$OS == 1, 'Dead', 'Alive')

# Define colors
colors <- c("#68d0ce", "#4955d0", "#f4a4c9", "#e09a2a", "#de6a73")
          
# Create violin plot
p <- ggplot(na.omit(df_clinical), aes(x = stage, y = riskscore)) +  
  geom_violin(aes(fill = stage), color = 'white', alpha = 0.5) +
  geom_boxplot(width = 0.1, fill = "white", alpha = 0.6) +
  scale_fill_manual(values = colors) +
  labs(x = "", y = "Risk Score") +
  theme_bw(base_size = 17)

# Add significance annotations
p <- p + geom_signif(
  comparisons = list(
    c('Stage I', 'Stage IV'),
    c('Stage II', 'Stage IV'),
    c('Stage III', 'Stage IV'),
    c('Stage I', 'Stage III')
  ),
  test = "wilcox.test",
  y_position = c(2.2, 2.3, 2.4, 2.5),
  textsize = 4
) + theme(legend.position = "")


#================================================
# 6.Expression of 7 Genes in Tumor and Normal Tissues
#================================================

# Add group information to df
df$group <- as.factor(ifelse(substring(df$sample, 14, 15) == '01', 'Tumor', 'Normal'))

# Perform Wilcoxon test for each gene
results_df <- data.frame(Gene = character(), P.Value = numeric(), stringsAsFactors = FALSE)

for (gene in intersect.gene.new) {
  test_result <- wilcox.test(
    x = df[[gene]][df$group == 'Tumor'],
    y = df[[gene]][df$group == 'Normal'],
    alternative = "two.sided"
  )
  results_df <- rbind(results_df, data.frame(Gene = gene, P.Value = test_result$p.value))
}

# Adjust p-values for multiple testing
results_df$adjusted_p_values <- p.adjust(results_df$P.Value, method = "fdr")

# Save results
write.csv(results_df, file = 'E:/项目一/结果表格/7gene_肿瘤癌旁检验结果.csv')

# Create boxplot for gene expression
df1 <- df[, c('group', intersect.gene.new)]
df11 <- melt(df1, id.vars = c("group"))
colnames(df11) <- c("group", "Gene", "Expression")
df11$group <- factor(df11$group)
df11$Expression <- log(df11$Expression)

ch1 <- ggplot(df11[df11$Expression > -1 & df11$Expression < 5, ], aes(x = Gene, y = Expression, fill = group, color = group)) +
  geom_boxplot(outlier.size = 0.4, position = position_dodge(1), size = 0.3, alpha = 0.3) +
  theme_classic(base_size = 15) +
  theme(panel.border = element_rect(color = 'black', fill = NA, size = 0.5)) +
  theme(axis.text.x = element_text(angle = 90))

p1 <- ch1 + scale_color_manual(values = c("#1597A5", "#6B5B95")) +
  scale_fill_manual(values = c("#1597A5", "#6B5B95"))

#================================================
# 7.Correlation Heatmap of 7 Lactate Metabolism Genes
#================================================

# Load necessary libraries
library(linkET)
library(cols4all)

# Calculate correlation matrix
dfrt <- df[, c(intersect.gene.new)]
M <- cor(dfrt)
p <- cor.mtest(dfrt)$p

# Create correlation heatmap
p1 <- corrplot(
  M,
  method = 'square',
  order = 'alphabet',
  type = 'lower',
  tl.col = 'black',
  tl.pos = 'ld',
  tl.cex = 0.8,
  tl.srt = 45,
  is.corr = FALSE,
  col = rev(COL2('PiYG', 10)),
  col.lim = c(-1, 1),
  cl.pos = 'b',
  p.mat = p,
  insig = "label_sig",
  sig.level = c(0.001, 0.01, 0.05),
  pch.cex = 1.5,
  addCoef.col = 'black'
)


#================================================
# 8.Heatmap of Risk Score and Clinical Information
#================================================

# Arrange data by risk score
dfrt <- arrange(df_clinical, riskscore)
dfrt <- distinct(dfrt, sample, .keep_all = TRUE)
dfrt <- dfrt[dfrt$mstage != '', ]
rt <- dfrt[, 4:10]
rownames(rt) <- dfrt$sample

# Prepare annotation data
group <- dfrt[, c(3, 11, 13:18)]
rownames(group) <- dfrt$sample

# Define annotation colors
ann_colors = list(
  status = c('Dead'="#F8D09D",'Alive'='#FBE6C7'),
  stage=c('Stage I'='#C4D8EB','Stage II'='#9CC5DA','Stage III'='#679EC2','Stage IV'='#4384B1'),
  mstage=c('M0'='#F9DEEB','M1'='#E9B3D1','MX'='#D272A5'),
  nstage=c('N0'='#D9D9EA','N1'='#BABAD8','N2'='#9B97C1','N3'='#7B79B1','NX'='#644E9A'),
  tstage=c('T1'='#E3EFCE','T2'='#B2D685','T3'='#7AB444','T4'='#4E8D37'),
  age=c('<65'="#FC8D59",'>=65'="#FDBB84")
)

pheatmap(t(rt),cluster_rows =F,cluster_cols =F,show_colnames =F,show_rownames = T,
         # fontsize_row =9,
         #fontsize=10,
         #fontface="italic",fontfamily= "Arial",
         #font = c(8, "italic"),
         fontfamily = "serif",
         color = c(colorRampPalette(colors = c("#248D89", "white"))(length(bk)/2),
                   colorRampPalette(colors = c("white", "#C26C35"))(length(bk)/2)),
         legend_breaks=seq(-2,2,2),
         breaks=bk,
         annotation_colors = ann_colors,
         annotation_col =group,
         scale = "row")
                       
#============================================================================
# 9.Validation of Risk Score Stratification in External Cohorts
#============================================================================
options(timeout=100000)
getOption('timeout')

#GSE73403
gset=AnnoProbe::geoChina('GSE73403')
pdata <- pData(gset[[1]])
GSE73403_osdata<-pdata[,c(2,39,43)]
colnames(GSE73403_osdata)=c('id','OS','OS.time')
GSE73403_osdata=type.convert(GSE73403_osdata)
GSE73403_osdata$OS.time=GSE73403_osdata$OS.time*365
write.csv(GSE73403_osdata,"E:/项目一/验证集数据/GSE73403/GSE73403_osdata.csv")


exp <- exprs(gset[[1]])
exp<-as.data.frame(exp)

GSE73403=openxlsx::read.xlsx("E:/项目一/验证集数据/GSE73403/GSE73403_exp.xlsx",sheet=2,rowNames=F)
GSE73403=GSE73403%>%
  filter(gene %in% intersect.gene.new) %>%
  t()
colnames(GSE73403)=GSE73403[1,]
GSE73403=GSE73403[2:70,]
GSE73403=as.data.frame(GSE73403)
GSE73403$id=rownames(GSE73403)

GSE73403=merge(GSE73403,GSE73403_osdata,by.y='id')
GSE73403=type.convert(GSE73403)

GSE73403$riskscore=GSE73403$CHEK2*coef[1]+
  GSE73403$LIPT1*coef[2]+
  GSE73403$NDUFA10*coef[4]+
  GSE73403$AGK*coef[5]+
  GSE73403$PNPLA2*coef[6]+GSE73403$GFM1*coef[7]

GSE73403$risk_level <- as.vector(ifelse(GSE73403$riskscore > median(GSE73403$riskscore), "High", "Low"))  

fit <- survfit(Surv(OS.time, OS) ~ risk_level, data =GSE73403)
res.cut <- surv_cutpoint(GSE73403, time = "OS.time", event = "OS",
                         variables = c("riskscore") 
)
summary(res.cut)
res.cat <- surv_categorize(res.cut)
fit <- survfit(Surv(OS.time,OS) ~riskscore, data = res.cat)
km2=ggsurvplot(
  fit,
  data = res.cat,
  censor.shape="|", censor.size = 4,
  conf.int = TRUE,
  conf.int.style = "ribbon",
  conf.int.alpha = 0.2,
  pval = TRUE,
  ggtheme = theme_survminer(), 
  palette = c("#b76d7b", "#87B5C1") ,
  surv.median.line = "hv",
  legend = "top",
  xlab = "OS_time(days)",
  ylab = "Survival probablity",
  title = "GSE73403",
  break.x.by = 1000,
  break.y.by = 0.2
)
write.csv(GSE73403,file='E:/项目一/验证集数据/GSE73403/GSE73403_riskscore.csv')

#GSE37745
gset=AnnoProbe::geoChina('GSE37745')
pdata <- pData(gset[[1]])
GSE37745_osdata<-pdata[,c(2,45,43)]
colnames(GSE37745_osdata)=c('id','OS','OS.time')
GSE37745_osdata=GSE37745_osdata %>%
  mutate(OS=case_when(OS=='yes'~1,OS=='no'~0))%>%
  type.convert()

write.csv(GSE37745_osdata,"E:/项目一/验证集数据/GSE37745/GSE37745_osdata.csv")

exp <- exprs(gset[[1]])
exp<-as.data.frame(exp)
write.csv(exp, file = "E:/项目一/验证集数据/GSE37745/GSE37745_exp.csv")

GSE37745=openxlsx::read.xlsx("E:/项目一/验证集数据/GSE37745/GSE37745_exp.xlsx",sheet=3,rowNames=F)
GSE37745=GSE37745%>%
  filter(gene %in% intersect.gene.new) %>%
  t()
colnames(GSE37745)=GSE37745[1,]
GSE37745=GSE37745[2:197,]
GSE37745=as.data.frame(GSE37745)
GSE37745$id=rownames(GSE37745)

GSE37745=merge(GSE37745,GSE37745_osdata,by.y='id')
GSE37745=type.convert(GSE37745)

GSE37745$riskscore=
  GSE37745$CHEK*coef[1]+
  GSE37745$LIPT1*coef[2]+
  GSE37745$TUFM*coef[3]+
  GSE37745$NDUFA10*coef[4]+
  GSE37745$AGK*coef[5]+
  GSE37745$PNPLA2*coef[6]+
  GSE37745$GFM1*coef[7]

GSE37745$risk_level <- as.vector(ifelse(GSE37745$riskscore > median(GSE37745$riskscore), "High", "Low"))  
fit <- survfit(Surv(OS.time, OS) ~ risk_level, data =GSE37745)
ggsurvplot(fit,data = GSE37745, pval = TRUE)

res.cut <- surv_cutpoint(GSE37745, time = "OS.time", event = "OS",
                         variables = c("riskscore") )
summary(res.cut)
res.cat <- surv_categorize(res.cut)
fit <- survfit(Surv(OS.time,OS) ~riskscore, data = res.cat)
km3=ggsurvplot(
  fit,
  data = res.cat,
  censor.shape="|", censor.size = 4,
  conf.int = TRUE,
  conf.int.style = "ribbon",
  conf.int.alpha = 0.2,
  pval = TRUE,
  ggtheme = theme_survminer(), 
  palette = c("#b76d7b", "#87B5C1") ,
  surv.median.line = "hv",
  legend = "top",
  xlab = "OS_time(days)",
  ylab = "Survival probablity",
  title = "GSE37745",
  break.x.by = 1000,
  break.y.by = 0.2
)
write.csv(GSE37745,file='E:/项目一/验证集数据/GSE37745/GSE37745_riskscore.csv')
