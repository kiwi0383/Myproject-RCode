#-------------------------------------------------------------------------------
#Linkages between the intratumor microbiota, LM.Sig and immunity
#-------------------------------------------------------------------------------

#===================================
#1.Data Preparation
#===================================
library(dplyr)
library(readr)
library(openxlsx)
library(vegan)
library(ggplot2)
library(rstatix)
library(ggpubr)

lusc_group<-openxlsx::read.xlsx('E:/项目一/data/fig1/df.xlsx') %>%
  select(2:14) %>%
  mutate(sample=paste0(sample,'A'))

micro<-read_tsv("E:/STAD/数据/mycobiome-master/Final_files/count_data_genus_raw_WIS_overlapping_fungi_bacteria_14494samples.tsv")
meta<-read_tsv("E:/STAD/数据/mycobiome-master/Final_files/metadata_genus_WIS_overlapping_fungi_bacteria_14494samples.tsv")

lusc_sampleid=meta %>%
  filter(tcga_sample_id %in% lusc_group$sample)%>%
  filter(experimental_strategy=='RNA-Seq') %>%
  select('tcga_sample_id','sampleid')

lusc_sampleid<-lusc_sampleid %>%
  filter(substring(tcga_sample_id,14,15)=='01')

lusc_micro=merge(lusc_sampleid,micro,by.y='sampleid')
names(lusc_group)[1]='tcga_sample_id'
lusc_micro=merge(lusc_micro,lusc_group,by.y='tcga_sample_id')

lusc_micro <- distinct(lusc_micro,sampleid, .keep_all = TRUE)
rownames(lusc_micro)=lusc_micro$sampleid

#=========================================
#2、Differential Microbiome Analysis
#=========================================

#-------------------------------------------------------------------------------
#2.1 Wilcoxon Rank Sum Test
#-------------------------------------------------------------------------------
lusc_micro=read.xlsx('E:/项目一/data/fig3/lusc_micro.xlsx')

microbiota<-colnames(lusc_micro[,2:1781])
results_df <- data.frame(microbiota = character(),  
                         P.Value = numeric(),  
                         stringsAsFactors = FALSE)
for (micro in microbiota) {  
  test_result <- wilcox.test(x = lusc_micro[[micro]][lusc_micro$risk_level == 'High'],  
                             y = lusc_micro[[micro]][lusc_micro$risk_level == 'Low'],  
                             alternative = "two.sided")  
  results_df <- rbind(results_df, data.frame(microbiota = micro, P.Value = test_result$p.value))  
}  

results_df$adjusted_p_values <- p.adjust(results_df$P.Value, method = "fdr")

bacteria0.05=results_df[results_df$adjusted_p_values<0.05,1]
bacteria0.01=results_df[results_df$P.Value<0.01,1]
bacteria0.05 <- Filter(function(x) !is.na(x), bacteria0.05)

wilcox.micro=results_df[results_df$adjusted_p_values<0.05,]
write.xlsx(wilcox.micro,file='E:/项目一/结果表格/wilcox.micro.xlsx')

#-------------------------------------------------------------------------------
#2.2 PCoA Analysis
#-------------------------------------------------------------------------------
dist = vegdist(lusc_micro[,c(bacteria0.05)], method="bray",na.rm = T)

dist <- as.matrix(dist)
pcoa = cmdscale(dist, eig=T)
eig = summary(eigenvals(pcoa))
axis = paste0("PCoA", 1:ncol(eig))
eig = data.frame(Axis = axis, t(eig)[, -3])

pco1 = round(eig[1, 3] * 100, 2)
pco2 = round(eig[2, 3] * 100, 2)

xlab = paste0("PCoA1 (",pco1,"%)")
ylab = paste0("PCoA2 (",pco2,"%)")

pcoa_points = as.data.frame(pcoa$points)
pcoa_points = data.frame(pcoa_points, group = lusc_micro$risk_level)

color=c("#129392","#8d5a85") 
p1<-ggplot(data=pcoa_points,aes(x=V1,y=V2))+
  theme_bw()+
  geom_point(aes(color = group), shape = 19, size=2)+
  theme(panel.grid = element_blank())+
  geom_vline(xintercept = 0,lty="dashed", size = 1, color = 'grey50')+
  geom_hline(yintercept = 0,lty="dashed", size = 1, color = 'grey50')+
  stat_ellipse(aes(color = group),level=0.95)+
  scale_color_manual(values = color) +
  scale_fill_manual(values = color)+
  theme(axis.title.x=element_text(size=15),
        axis.title.y=element_text(size=15,angle=90),
        
        panel.grid=element_blank())
p1

Adonis <- adonis2(dist~ pcoa_points$group,
                  distance = "bray",
                  permutations = 999)
Adonis
adonis <- paste0("adonis R2: ",round(Adonis$R2,2), 
                 "; P-value: ", Adonis$`Pr(>F)`)
p2 <- p1 + labs(subtitle = adonis)
p2+theme_bw(base_size = 15)+theme(panel.grid=element_blank())

#-------------------------------------------------------------------------------
#2.3lefse
#-------------------------------------------------------------------------------
otu<-lusc_micro[,2:1781]
otu<-t(otu) %>%
  as.data.frame()
group=lusc_micro[,c(1,1782)]
names(group)=c('sample','group')
otu=rbind(t(group),otu)
write.table(otu,file='E:/项目一/data/fig3/lefse.txt')
#画图
ggplot(results_df, aes(lda, reorder(micro, lda),fill = group)) +
  geom_bar(stat = "identity", color = "black", linewidth = 0.1,alpha=1) +
  #geom_text(aes(x = ifelse(lda > 0, -1, 1),label = micro),color = "black",size=4) +  
  #scale_fill_viridis(begin = 0, end = 0.6, option = "D",alpha=0.6) +  
  scale_fill_manual(values = c("#aabcdb","#e2d1b5")) +
  labs(x = "LDA SCORE",y="micro") +  theme_classic() +ggtitle("lefse")+ 
  theme(legend.position = "right") +  
  theme(axis.text = element_text(size = 12)) +
  theme(axis.title.x = element_text(vjust = 0, size = 15, face = "plain")) + # 调整xlab字体大小
  theme(axis.title.y = element_text(vjust = 0, size = 15, face = "plain")) + # 调整ylab字体大小  
  theme(plot.title = element_text(size = 15, face = "bold")) +
  theme(axis.line = element_line(color = "black",linewidth = 0.8)) # 调整坐标轴样式

#-------------------------------------------------------------------------------
#2.4 DESeq2 Analysis
#-------------------------------------------------------------------------------
library(DESeq2)
otu<-lusc_micro[,2:1781]
otu<-t(otu) %>%
  as.data.frame()
colData=data.frame(group=as.factor(lusc_micro$risk_level),row.names = lusc_micro$id1)

dds <- DESeqDataSetFromMatrix(countData = otu,
                              colData = colData,
                              design = ~group)
dds<-DESeq(dds)
resultsNames(dds)
res<-results(dds,contrast=c("group",'Low','High'))

resOrdered<-res[order(res$padj),]

#火山图
DEG_deseq2<-na.omit(DEG)
write.xlsx(DEG_deseq2,file='E:/项目一/结果表格/deseq2.micro.xlsx')

logFC=0.585
P.Value=0.05
k1<-(DEG_deseq2$pvalue<P.Value)&(DEG_deseq2$log2FoldChange<-logFC)
k2<-(DEG_deseq2$pvalue<P.Value)&(DEG_deseq2$log2FoldChange>logFC)
DEG_deseq2<-mutate(DEG_deseq2,change=ifelse(k1,"down",ifelse(k2,"up","stable")))
table(DEG_deseq2$change)

DEG_deseq2 <- DEG_deseq2 |>  
  mutate(change = case_when(log2FoldChange > 0.585 & padj < 0.05 ~ "Up",
                            abs(log2FoldChange) <  0.585 | padj > 0.05 ~ "None",
                            log2FoldChange < - 0.585 & padj < 0.05 ~ "Down"))
DEG_deseq2$change |> as.factor()
table(DEG_deseq2$change)

p<-ggplot(data=DEG_deseq2,
          aes(x=log2FoldChange,
              y=-log10(pvalue)))+
  geom_point(alpha=0.4,size=3.5,
             aes(color=change))+
  ylab("-log10(Pvalue)")+
  scale_color_manual(values=c("blue4","grey","red3"))+
  geom_vline(xintercept=c(-logFC,logFC),lty=4,col="black",lwd=0.8)+
  geom_hline(yintercept=-log10(P.Value),lty=4,col="black",lwd=0.8)+
  theme_bw(base_size = 15)+ ylim(c(-1, 17)) 
p

DEG<-as.data.frame(resOrdered)             

DEG$micro=rownames(DEG)
deseq.micro=DEG[DEG$padj<0.05,7]
deseq.micro<- Filter(function(x) !is.na(x), deseq.micro)
wilc.micro=bacteria0.05


#-------------------------------------------------------------------------------
#2.5 Venn plot of the collection of differential microorganisms
#-------------------------------------------------------------------------------

if(!requireNamespace("ggvenn")){remotes::install_github("yanlinlin82/ggvenn")}
library(ggvenn)
mylist=list(deseq2.micro=deseq2.micro,wilcox.micro=wilcox.micro,lefse.micro=lefse.micro)
p2 <- ggvenn(data = mylist[1:3], 
             columns = NULL, 
             show_elements = FALSE, 
             show_percentage = FALSE, 
             fill_color = c("#b6b6d4","#b0cde6","#b9e1f5") , 
             fill_alpha = 0.3,
             stroke_color ="black", 
             stroke_alpha = 0.5, 
             stroke_size = 0.5,
             set_name_color ="black", 
             set_name_size = 6,
             text_color ="black", 
             text_size = 5 )
p2

#===============================================================================
#3、Microbial abundance heatmap
#===============================================================================
library(pheatmap)
library(RColorBrewer)

rownames(lusc_micro)=lusc_micro$id1
otu<-arrange(lusc_micro,riskscore)
rt=otu[,c(bacteria0.05)]
grouprt=data.frame(group=otu$risk_level,row.names = rownames(otu))

bk <- c(seq(-2,-0.1,by=0.01),seq(0,2,by=0.01))
ann_colors = list(
  #cluster = c('cluster 1'="#4E79A7",'cluster 2'='#F28E2B','normal'='#b76d7b'),
  group=c("High"="#129392",'Low'="#8d5a85")
)

col_1=viridis::mako(13,direction=-1)[1:10]
col_2=scico::scico(10,direction=-1,palette='acton')
mycolor=c(rev(col_1[1:8]),'white',col_2[1:8])
pheatmap(t(rt),cluster_rows =F,cluster_cols =F,show_colnames =F,show_rownames = F,
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
         annotation_col =grouprt,
         scale = "row")


#===============================================================================
#4、Microbiome-Immune Mediation Analysis
#===============================================================================

library(mediation)

data=read.xlsx('E:/项目一/结果表格/wilcox.micro.xlsx')
wilcox.micro=data[data$adjusted_p_values<0.05,1]

data=read.xlsx('E:/项目一/结果表格/deseq2.micro.xlsx')
deseq2.micro=data[data$padj<0.05,7]

data=read.xlsx('E:/项目一/data/fig3/lefse.micro.xlsx')
lefse.micro=data$micro

a=intersect(deseq2.micro,wilcox.micro)
intersect.micro=intersect(a,lefse.micro)

data=read.csv('E:/项目一/结果表格/fig1/multicox.csv')
lmgene=data$term


lusc=read.csv('E:/项目一/data/fig2/lusc.csv',header = T) %>%
  column_to_rownames(.,var='X')

lusc$sample=substring(lusc$sample,1,12)
lusc=distinct(lusc,sample,.keep_all = T)

#读取反卷积免疫细胞丰度文件
infiltration_estimation=openxlsx::read.xlsx('E:/项目一/data/fig2/infiltration_estimation_for_tcga(01样本).xlsx')
riskscore<-openxlsx::read.xlsx('E:/项目一/data/fig1/df.xlsx') %>%
  select(2:14)
dfcell=riskscore[,c('sample','riskscore','risk_level')]
dfcell$sample=substring(dfcell$sample,1,12)
names(infiltration_estimation)[1]='sample'
dfcell=merge(dfcell,infiltration_estimation,by='sample')

df=type.convert(dfcell)
df=df[,c(1:9,32:63,65:121)]

data=merge(lusc[,c(lmgene,'sample')],df)

lusc_micro=read.xlsx('E:/项目一/data/fig3/lusc_micro.xlsx')
lusc_micro=merge(lusc_micro,micro_clin[,1:2])
lusc_micro$sample=substring(lusc_micro$id,1,12)
micro=lusc_micro[,c(intersect.micro,'sample')]

dfmediation=merge(micro,data)
dfmediation=distinct(dfmediation,sample,.keep_all = T)
dfmediation=type.convert(dfmediation)

colnames(dfmediation)
cells=colnames(dfmediation)[c(20:114)]
names(dfmediation)[c(2:10)]=c('Terrabacter','Flammeovirga','Cyanothece','Acidibacillus','Lachnoclostridium',
                              'Gallibacterium','Paraburkholderia','Streptococcus','Gemmata')
write.xlsx(dfmediation,"E:/项目一/data/fig3/dfmediation.xlsx")

#微生物-lmgene-cell
dfmediation=read.xlsx("E:/项目一/data/fig3/dfmediation.xlsx")

microname=colnames(dfmediation)[c(2:10)]
cells=colnames(dfmediation)[c(20:114)]
cell1=cells[c(2,3,42,46)]

#先找出和微生物显著相关的lmgene
results_df <- data.frame(Gene = character(),  
                         cor = numeric(), p=numeric(), 
                         stringsAsFactors = FALSE)
for (gene in lmgene) {  
  test_result <- cor.test(x = dfmediation[[gene]],  
                          y = dfmediation$Lachnoclostridium)  
  results_df <- rbind(results_df, data.frame(Gene = gene, cor = test_result$estimate,p=test_result$p.value))  
} 

lmgene.sig=results_df[results_df$p<0.05,1]
microname='Lachnoclostridium'
cell1


results_list <- data.frame(micro=character(),lmgene = character(),  
                           cell1 = character(),  
                           ACME = numeric(), ADM=numeric(), ACMEp = numeric(), ADMp=numeric(), Prop.mediated.p=numeric(),
                           stringsAsFactors = T)

for(microbiome in microname) {
  for (cell in cell1) {
    # 内层循环处理每个lmgene
    for (gene in lmgene.sig) {
      # 构建模型M
      model.M <- lm(as.formula(paste(gene,'~',microbiome)), data = dfmediation)
      
      # 构建模型Y
      model.Y <- lm(as.formula(paste(cell,'~',microbiome,'+',gene)), data = dfmediation)
      
      # 执行中介分析
      results <- mediate(model.M, model.Y, treat=microbiome, mediator=gene, boot=TRUE, sims=500)
      
      # 将结果添加到results_list数据框中
      results_list <- rbind(results_list, data.frame(micro=microbiome,
                                                     lmgene= gene,
                                                     cell1 = cell, 
                                                     ACME = results$d.avg,  # （Average Causal Mediation Effect）
                                                     ADE = results$z.avg,  # 直接效应
                                                     ACMEp = results$d.avg.p,
                                                     ADEp = results$z.avg.p,
                                                     Prop.mediated.p=results$n.avg.p))
    }
  }
}
# 查看结果

allresult=results_list
allresult=rbind(allresult,results_list)
write.xlsx(allresult,"E:/项目一/data/fig3/微生物-lmgene-免疫.xlsx")

write.xlsx(results_list,"E:/项目一/data/fig3/Lachnoclostridium-lmgene-免疫.xlsx")

#sankey plot
data=read.xlsx("E:/项目一/data/fig3/Lachnoclostridium-lmgene-免疫.xlsx",sheet=2)
sankey_data <- data %>%
  make_long(micro, lmgene, cell)
sankey_data$node <- factor(sankey_data$node)
sankey_data$next_node <- factor(sankey_data$next_node)
sankey_data$x <- factor(sankey_data$x)

others <-  c("Lachnoclostridium"="#e2d1b5")
color_palette <- colorRampPalette(c("#944D55","#B5555A","#D05C5F","#E46665","#F57B70","#FC9780",
                                    "#FCAD98","#FDC5B2","#FDD9CA","#FFECE5","grey90"))

color_palette <- colorRampPalette(c("#b79bac","#e2d1b8","#9193b1"))


colors <- c(others, setNames(color_palette(6), unique(data$lmgene)), 
            setNames(color_palette(4), cell1))
ggplot(sankey_data, 
       aes(x = x, 
           next_x = next_x, 
           node = node, 
           next_node = next_node, 
           fill = factor(node))) +
  geom_sankey(width = 0.3, 
              smooth = 5,
              space = 15,
              na.rm = TRUE,
              position = "identity",
              flow.alpha = 0.4, 
              node.color = "transparent") +
  geom_sankey_text(aes(label = node),
                   width = 0.3, 
                   space = 15,
                   position = "identity", 
                   size = 4, 
                   color = "black", 
                   hjust =0) +
  theme_sankey(base_size = 10) +
  scale_fill_manual(values = colors) +
  theme(legend.position = "none")

#spearman test
attach(dfmediation)
cor.test(GFM1,Lachnoclostridium,method="spearman",exact = F)
cor.test(GFM1,`T.cell.CD8_TIMER`,method="spearman",exact = F)
cor.test(`T.cell.CD8_TIMER`,Lachnoclostridium,method="spearman",exact = F)

cor.test(PNPLA2,Lachnoclostridium,method="spearman",exact = F)
cor.test(PNPLA2,`T.cell.CD4_TIMER`,method="spearman",exact = F)
cor.test(`T.cell.CD4_TIMER`,Lachnoclostridium,method="spearman",exact = F)

cor.test(GFM1,Lachnoclostridium,method="spearman",exact = F)
cor.test(GFM1,`NK.cell_MCPCOUNTER`,method="spearman",exact = F)
cor.test(`NK.cell_MCPCOUNTER`,Lachnoclostridium,method="spearman",exact = F)

cor.test(CHEK2,Lachnoclostridium,method="spearman",exact = F)
cor.test(CHEK2,`Myeloid.dendritic.cell_MCPCOUNTER`,method="spearman",exact = F)
cor.test(`Myeloid.dendritic.cell_MCPCOUNTER`,Lachnoclostridium,method="spearman",exact = F)

detach(dfmediation)


#===============================================================================
#5、Microbiome-Immune Correlations
#===============================================================================

#5.1 ESTIMATE \ TIDE & lachnoclostridium
riskscore<-openxlsx::read.xlsx('E:/项目一/data/fig1/df.xlsx') %>%
  dplyr::select(2:14)
estimate<-read.table("E:/TCGA/TCGA_Annovar/TCGA-LUSC/TCGA-LUSC_estimate_scores.txt",header = T)
estimate$sample=substring(estimate$ID,1,15)
df1=merge(riskscore,estimate,by.y='sample')
df1$Tumourpurity = cos (0.6049872018+0.0001467884 *df1$ESTIMATEScore)

tide=read.csv('E:/项目一/data/fig2/LUSC-TIDE.csv')
tide$sample=substring(tide$Patient,1,15)
df2=merge(riskscore,tide,by.y='sample')

df=merge(df1,df2)
df$sample=substring(df$sample,1,12)
df=merge(df,dfmediation[,c('sample','Lachnoclostridium')])

df1=df[,c(4:11,15:18,22,23,25:27,29,34)]
colnames(df1)

results_df <- data.frame(Gene = character(),  
                         cor = numeric(), p=numeric(), 
                         stringsAsFactors = FALSE)
for (gene in colnames(df1)[c(1:18)]) {  
  test_result <- cor.test(x = df1[[gene]],  
                          y = df1$Lachnoclostridium,method = 'spearman',exact = F)  
  # 将检验结果存储在列表中，键名为基因名  
  results_df <- rbind(results_df, data.frame(Gene = gene, cor = test_result$estimate,p=test_result$p.value))  
}  

results_df$adjusted_p_values <- p.adjust(results_df$p, method = "BH")

results_df$significance <- ifelse(results_df$adjusted_p_values < 0.001, "***",
                                  ifelse(results_df$adjusted_p_values < 0.01, "** ",
                                         ifelse(results_df$adjusted_p_values < 0.05, "* ",
                                                "NS")))

results_df$cor1<-cut(abs(results_df$cor),#绝对值
                     breaks=c(0,0.3,0.5,0.7,0.9,1),
                     labels=c("<0.3","0.3-0.5","0.5-0.7","0.7-0.9",">0.9"),
                     right=FALSE)#right=FALSE表示表示区间为左闭右开
results_df$pvalue1<-cut(results_df$pvalue,
                        breaks=c(0,0.001,0.01,0.05,1),
                        labels=c("<0.001","<0.01","<0.05",">0.05"),
                        right=FALSE)

dat=results_df[order(results_df$cor),]
dat$Gene=factor(dat$Gene,levels=dat$Gene)
ggplot(dat,aes(x=cor,y=Gene,color=significance))+
  scale_color_manual(name="significance",
                     values=c("#D4AF37","#6C5B7B","darkgrey"))+
  geom_segment(aes(x=0,y=Gene,xend=cor,yend=Gene),size=1)+
  geom_point(aes(size=cor1))+
  geom_point(aes(size=cor1)) +
  #theme_bw(base_size = 15)+
  theme_classic(base_size = 15)+
  labs(size="Cor")


#===============================================================================
#6、lachnoclostridium & GSVA score of LM-genes
#===============================================================================
# Run GSVA with LM gene signature
library(dplyr) 
library(data.table) 
library(GSVA)
library(GSEABase) 

lusc=read.csv('E:/项目一/data/fig2/lusc.csv',header = T) %>%
  column_to_rownames(.,var='X')

lmgenelist=list(gene=colnames(lusc)[1:228])

data <- log(fpkm[,2:19532]+0.01)
fpkm=column_to_rownames(fpkm,'X')

gsvaP <- gsvaParam(exprData = as.matrix(t(fpkm[,2:19531])),  
                   geneSets = lmgenelist,
                   minSize = 1, maxSize = Inf, kcdf = "Gaussian", tau = 1, maxDiff = TRUE, absRanking = FALSE)

gsva_data <- gsva(gsvaP)
gsva_data <- t(gsva_data)

gsva_data=gsva_data %>%
  as.data.frame() %>%
  mutate(sample=substring(rownames(.),1,12))

df=merge(gsva_data,dfmediation)
names(df)[2]='GSVA'

# Correlate GSVA scores with microbiome
library(ggpubr)
library(ggExtra)
p=ggscatter(df[df$Lachnoclostridium<15000,],x='Lachnoclostridium',y='GSVA',
            size=1.3,add='reg.line',
            add.params=list(color="#734C81",fill="#D7C3DF",size=1),conf.int = T)+
  stat_cor(method='spearman',label.x = 1,label.y=0.5,size=5)+theme_bw(base_size = 15)

ggMarginal(
  p,
  type =c('density'),margins = 'both',size = 5,colour ="#734C81",fill="#D7C3DF",
  groupColour = F,groupFill = F) 