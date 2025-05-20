#========================================================================
# Figure 3: Evaluation of LM.Sig at Single-Cell Resolution
#========================================================================

#-----------------------------
# 0. Environment Setup
#-----------------------------
library(Seurat)
library(dplyr)
library(patchwork)
library(ggplot2)
library(data.table)
library(irGSEA)
library(harmony)
library(COSG)
library(ClusterGVis)
library(org.Hs.eg.db)
library(GSVA)
library(GSEABase)
library(msigdbr)

# Install required packages if missing
if (!require("promises")) install.packages("promises", type = "source")
if (!require("qs")) BiocManager::install("qs", force = TRUE)
if (!require("jjAnno")) devtools::install_github("junjunlab/jjAnno")

#-----------------------------
# 1. Data Loading and Preprocessing
#-----------------------------
# Set data directory (use relative path for portability)
data_dir <- "data/GSE148071_RAW/"
samples <- list.files(data_dir)

# Read and merge multiple samples
sceList <- lapply(samples, function(sample) {
  ct <- fread(file.path(data_dir, sample), data.table = FALSE)
  rownames(ct) <- ct[,1]
  colnames(ct) <- paste(gsub(".*_([Pp][0-9]+)_.*", "\\1", sample),
                        colnames(ct), sep = '_')
  ct <- ct[,-1]
  CreateSeuratObject(counts = ct,
                     project = sample,
                     min.cells = 5,
                     min.features = 500)
})

# Merge all samples
sce.all <- merge(x = sceList[[1]], 
                 y = sceList[-1],
                 add.cell.ids = samples)

# Check merged object
print(dim(sce.all[["RNA"]]$counts))
table(sce.all$orig.ident) # Cell counts per sample

# Save merged data
save(sce.all, file = "data/sce.all.RData")

#-----------------------------
# 2. Quality Control
#-----------------------------
seurat.data <- sce.all

# Calculate mitochondrial percentage
seurat.data[["percent.mt"]] <- PercentageFeatureSet(seurat.data, pattern = "^MT-")

# Visualize QC metrics
VlnPlot(seurat.data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"))
plot1 <- FeatureScatter(seurat.data, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(seurat.data, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# Filter cells
seurat.data <- subset(seurat.data,
                      subset = nFeature_RNA > 200 & nFeature_RNA < 5000 &
                        nCount_RNA > 500 & percent.mt < 20)

# Normalization and variable features
seurat.data <- NormalizeData(seurat.data)
seurat.data <- FindVariableFeatures(seurat.data,
                                    selection.method = "vst",
                                    nfeatures = 1000)
seurat.data <- ScaleData(seurat.data, features = VariableFeatures(seurat.data))

#-----------------------------
# 3. Dimensionality Reduction and Clustering
#-----------------------------
# Initial PCA/UMAP
seurat.data <- seurat.data %>% 
  RunPCA(npcs = 30, verbose = FALSE) %>% 
  RunUMAP(reduction = "pca", dims = 1:30, verbose = FALSE)

# Batch effect visualization
p1.compare <- wrap_plots(
  DimPlot(seurat.data, reduction = "pca", group.by = "orig.ident") + 
    NoAxes() + ggtitle("Before_PCA"),
  DimPlot(seurat.data, reduction = "umap", group.by = "orig.ident") + 
    NoAxes() + ggtitle("Before_UMAP"),
  guides = "collect"
)

# Harmony integration
seurat.data <- seurat.data %>% 
  RunHarmony("orig.ident", plot_convergence = TRUE) %>% 
  RunUMAP(reduction = "harmony", dims = 1:20, verbose = FALSE) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:20)

# Post-harmony visualization
p2.compare <- wrap_plots(
  DimPlot(seurat.data, reduction = "harmony", group.by = "orig.ident") + 
    NoAxes() + ggtitle("After_PCA (harmony)"),
  DimPlot(seurat.data, reduction = "umap", group.by = "orig.ident") + 
    NoAxes() + ggtitle("After_UMAP"),
  guides = "collect"
)

wrap_plots(p1.compare, p2.compare, ncol = 1)

# Cluster at multiple resolutions
for (res in c(0.05, 0.1, 0.3, 0.5, 0.8, 1, 1.2, 1.4, 1.5, 2)) {
  seurat.data <- FindClusters(seurat.data, resolution = res, algorithm = 1)
}

# Visualize clustering results
cluster_umap <- wrap_plots(ncol = 5,lapply(c(0.05, 0.1, 0.3, 0.5, 0.8), function(res) {
  DimPlot(seurat.data, reduction = "umap", 
          group.by = paste0("RNA_snn_res.", res), label = TRUE) + NoAxes()
})
cluster_umap

# Set final resolution (0.5)
Idents(seurat.data) <- "RNA_snn_res.0.5"

#-----------------------------
# 4. Cell Type Annotation
#-----------------------------
# Marker gene identification using COSG
marker_cosg <- COSG::cosg(
  seurat.data,
  groups = 'all',
  assay = 'RNA',
  slot = 'data',
  mu = 1,
  expressed_pct = 0.1,
  remove_lowly_expressed = TRUE,
  n_genes_user = 200)

# Find all markers
pbmc.markers <- FindAllMarkers(seurat.data,
                               only.pos = TRUE,
                               min.pct = 0.25,
                               logfc.threshold = 0.25)

#T
DotPlot(seurat.data,features=c('CD3D','CD3E', 'CD8A','PTPRC','CD4','CD2'))+RotatedAxis()
#B
DotPlot(seurat.data,features=c('CD19', 'CD79A', 'MS4A1', 'CD20'))+RotatedAxis()
#Plasma cells 
DotPlot(seurat.data,features=c('IGHG1','MZB1','SDC1','CD79A'))+RotatedAxis()

#nk
DotPlot(seurat.data,features=c('FGFBP2', 'FCG3RA', 'CX3CR1'))+RotatedAxis()
#Fibroblasts
DotPlot(seurat.data,features=c('FGF7','MME','GSN','LUM','DCN'))+RotatedAxis()

#mast
DotPlot(seurat.data,features=c('CPA3', 'CST3', 'KIT', 'TPSAB1', 'TPSB2', 'MS4A2'))+RotatedAxis()
#Endothelial
DotPlot(seurat.data,features=c('PECAM1', 'VWF','PLVAP'))+RotatedAxis()
#Epithelial
DotPlot(seurat.data,features=c('EPCAM','CDH1','TJP1','TPH1','KRT8'))+RotatedAxis()
#Monocyte:S100A8, S100A9
DotPlot(seurat.data,features=c('S100A8', 'S100A9','CD68','LYZ','CD14','FCGR3A'))+RotatedAxis()
#Macrophage:APOE, C1QA, C1QB
DotPlot(seurat.data,features=c('APOE', 'C1QA', 'C1QB','CSF1R','CD68','CD163','CD14'))+RotatedAxis()
#alveolar
DotPlot(seurat.data,features=c('CLDN18'))+RotatedAxis()
#squamous cell
DotPlot(seurat.data,features=c('CLCA2','LGALS7','SERPINB5','SOX2','SPRR3'))+RotatedAxis()
#Proliferative cell
DotPlot(seurat.data,features=c('CCND1','Cdc20','Mki67','TOP2A'))+RotatedAxis()

# Cell type annotation based on marker genes
seurat.data <- RenameIdents(seurat.data,
                            "0" = "Epithelial",
                            "1" = "Epithelial",
                            "2" = "Macrophage",
                            "3" = "Proliferating cell",
                            "4" = "Epithelial",
                            "5" = "Monocyte",
                            "6" = "Basal cell",
                            "7" = "T cell",
                            "8" = "Fibroblast",
                            "9" = "Alveolar cell",
                            "10" = "B cell",
                            "11" = "Basal cell",
                            "12" = "Fibroblast",
                            "13" = "Endothelial",
                            "14" = "NK Cells",
                            "15" = "Mast cell",
                            "16" = "Ciliated cell",
                            "17" = "Fibroblast",
                            "18" = "Proliferating cell",
                            "19" = "NK Cells")

seurat.data$celltype <- Idents(seurat.data)

# Visualization
my_colors <- c("#CAB2D6", "#00A087FF", "#4682B4FF", "#F39B7FFF", "#8491B4FF",
               "#E64B35FF", "#91D1C2FF", "#87CEEBFF", "#FFDAB9FF", "#B09C85FF",
               "#FDBF6F", "#FF7F00", "#FFDAB9FF", "#6A3D9A", "#FFD700FF",
               "#90EE90FF", "#FFA07AFF", "#87CEEBFF", "#FA8072FF", "#4682B4FF")

DimPlot(seurat.data, reduction = "umap", group.by = "celltype", label = TRUE) +
  scale_color_manual(values = my_colors) +
  theme(panel.border = element_rect(fill = NA, color = "black", size = 1),
        legend.position = "right") +
  ggtitle("Cell Type UMAP")

#-----------------------------
# 5. Risk Score Analysis
#-----------------------------
# Load risk score coefficients
dfcoef <- read.csv("results/fig1/multicox.csv")
gene7 <- dfcoef$term
coef <- dfcoef$estimate

# Visualize gene expression
FeaturePlot(seurat.data, features = gene7)
DotPlot(seurat.data, features = gene7) + RotatedAxis()

# Calculate risk score
lactate_expression <- FetchData(seurat.data, vars = gene7)
lactate_expression$riskscore <- as.matrix(lactate_expression) %*% coef

# Add to metadata
seurat.data$risk_score <- lactate_expression$riskscore
seurat.data$risk_level <- ifelse(seurat.data$risk_score > median(seurat.data$risk_score), 
                                 "High", "Low")

# Visualization
FeaturePlot(seurat.data, features = "risk_score") +
  scale_color_gradientn(colors = c("blue", "white"))

VlnPlot(seurat.data, features = gene7, group.by = "risk_level")

# Save processed data
save(seurat.data, file = "data/seurat.data.RData")


#-----------------------------
# 6. Heatmap showing the expression of five marker genes in each cell type
#-----------------------------
library(ClusterGVis)
library(org.Hs.eg.db)
library(ggplot2)
BiocManager::install("jjAnno",force = TRUE)
devtools::install_github("junjunlab/jjAnno")
library(jjAnno)

pbmc.markers.all <- Seurat::FindAllMarkers(seurat.data,
                                           only.pos = TRUE,
                                           min.pct = 0.25,
                                           logfc.threshold = 0.25)
save(pbmc.markers.all,file='E:/项目一/结果表格/pbmc.markers.all.RData')

pbmc.markers <- pbmc.markers.all %>%  
  dplyr::group_by(cluster) %>%  dplyr::top_n(n = 20, wt = avg_log2FC)


st.data1 <- prepareDataFromscRNA(object = seurat.data,
                                 diffData = pbmc.markers,
                                 showAverage=TRUE)

enrich <- enrichCluster(object = st.data1,
                        OrgDb = org.Hs.eg.db,
                        type = "BP",
                        organism = "hsa",
                        pvalueCutoff = 0.5,
                        topn = 3,
                        seed = 5201314)
head(enrich, 3)


pbmc.markers1 <- pbmc.markers.all %>%  
  dplyr::group_by(cluster) %>%  
  dplyr::top_n(n = 5, wt = avg_log2FC)

table(Idents(seurat.data))
sce_down<-subset(seurat.data,downsample=1000)
table(Idents(sce_down))

sce_down<-subset(seurat.data,subset=nFeature_RNA>200 & nFeature_RNA<5000 &
                   nCount_RNA> 500 & percent.mt<10)
table(Idents(sce_down))

st.data <- prepareDataFromscRNA(object = sce_down,
                                diffData = pbmc.markers1,
                                showAverage = FALSE)

visCluster(object = st.data,
           plot.type = "line")

pbmc.markers2 <- pbmc.markers.all %>%  
  dplyr::group_by(cluster) %>%  
  dplyr::top_n(n = 3, wt = avg_log2FC)

visCluster(object = st.data,
           plot.type = "heatmap",
           markGenes = unique(pbmc.markers2$gene),
           column_title_rot = 45,
           cluster.order = 1:13)

pdf('tmp4.pdf',height = 8,width =14,onefile = F)
visCluster(object = st.data,
           plot.type = "both",
           column_title_rot = 45,
           show_row_dend = T,
           markGenes = unique(pbmc.markers2$gene),
           markGenes.side = "left",
           annoTerm.data = enrich,add.line = F,
           line.side = "left",
           cluster.order = c(1:13),
           sample.col=jjAnno::useMyCol("paired",n=13),
           go.col = rep(jjAnno::useMyCol("paired",n = 13),each = 3),
           add.bar = T)
dev.off()

pdf('tmp5.pdf',height=10,width=8,onefile=F)
visCluster(object=st.data,
           plot.type="heatmap",
           markGenes=unique(pbmc.markers1$gene),
           column_title_rot=45,
           cluster.order=1:9,
           sample.col=jjAnno::useMyCol("paired",n=13))
dev.off()


#-----------------------------
# 7. Pathway Analysis (GSVA)
#-----------------------------
# Prepare gene sets
h.human <- msigdbr(species = "Homo sapiens", category = "H")
h.sets <- split(h.human$gene_symbol, h.human$gs_name)

# Aggregate expression by cell type
av <- AggregateExpression(seurat.data, group.by = "celltype", assays = "RNA")[[1]]

# Run GSVA
gsva_res <- gsva(as.matrix(av), h.sets, 
                 min.sz = 1, max.sz = Inf, 
                 kcdf = "Poisson", parallel.sz = 4)

# Visualization
rownames(gsva_res) <- gsub("^HALLMARK_", "", rownames(gsva_res))
pheatmap(gsva_res,
         cluster_rows = FALSE,
         cluster_cols = TRUE,
         show_colnames = TRUE,
         fontfamily = "serif",
         color = colorRampPalette(c("#248D89", "white", "#C26C35"))(50),
         annotation_col = data.frame(row.names = colnames(gsva_res),
                                     CellType = colnames(gsva_res)))

#-----------------------------
# 7. Single-Cell Metabolism Analysis
#-----------------------------
# Subset cells for computational efficiency
set.seed(123)
sample_CB <- seurat.data@meta.data %>% 
  group_by(celltype) %>% 
  sample_frac(0.3)
sce_down <- subset(seurat.data, cells = sample_CB$CB)


sc.metabolism.SeuratV5<-function(obj,method="ssGSEA",imputation=F,ncores=2,
                                 metabolism.type="KEGG")
{
  countexp <- obj@assays$RNA$counts
  countexp <- data.frame(as.matrix(countexp))
  signatures_KEGG_metab<-system.file("data","KEGG_metabolism_nc.gmt",
                                     package="scMetabolism")
  signatures_REACTOME_metab<-system.file("data","REACTOME_metabolism.gmt",
                                         package="scMetabolism")
  if(metabolism.type=="KEGG"){
    gmtFile<-signatures_KEGG_metab
    cat("Yourchoiceis:KEGG\n")
  }
  if(metabolism.type=="REACTOME"){
    gmtFile<-signatures_REACTOME_metab
    cat("Yourchoiceis:REACTOME\n")
  }
  if(imputation==F){
    countexp2<-countexp
  }
  if(imputation==T){
    cat("Startimputation...\n")
    cat("Citation:GeorgeC.Linderman,JunZhao,YuvalKluger.Zero-preservingimputationofscRNA-seqdatausinglow-rankapproximation.bioRxiv.doi:https://doi.org/10.1101/397588\n")
    result.completed<-alra(as.matrix(countexp))
    countexp2<-result.completed[[3]]
    row.names(countexp2)<-row.names(countexp)
  }
  cat("Startquantifythemetabolismactivity...\n")
  if(method=="VISION"){
    library(VISION)
    n.umi<-colSums(countexp2)
    scaled_counts<-t(t(countexp2)/n.umi)*median(n.umi)
    vis<-Vision(scaled_counts,signatures=gmtFile)
    options(mc.cores=ncores)
    vis<-analyze(vis)
    signature_exp<-data.frame(t(vis@SigScores))
  }
  if(method=="AUCell"){
    library(AUCell)
    library(GSEABase)
    cells_rankings<-AUCell_buildRankings(as.matrix(countexp2),
                                         nCores=ncores,plotStats=F)
    geneSets<-getGmt(gmtFile)
    cells_AUC<-AUCell_calcAUC(geneSets,cells_rankings)
    signature_exp<-data.frame(getAUC(cells_AUC))
  }
  if(method=="ssGSEA"){
    library(GSVA)
    library(GSEABase)
    geneSets<-getGmt(gmtFile)
    gsvap<-ssgseaParam(exprData = as.matrix(countexp2),geneSets,
                       minSize = 1, maxSize = Inf)
    gsva_es<-gsva(gsvap)
    signature_exp<-data.frame(gsva_es)
  }
  
  if(method=="GSVA"){
    library(GSVA)
    library(GSEABase)
    geneSets<-getGmt(gmtFile)
    gsvaP<- gsvaParam(
      exprData = as.matrix(countexp2),
      geneSets,
      minSize = 1, maxSize = Inf, kcdf = "Poisson",parallel.sz=ncores)
    gsva_es <- gsva(gsvaP)
    signature_exp<-data.frame(gsva_es)
  }
  cat("\nPleaseCite:\nYingchengWu,QiangGao,etal.CancerDiscovery.2021.\nhttps://pubmed.ncbi.nlm.nih.gov/34417225/\n\n")
  obj@assays$METABOLISM$score<-signature_exp
  obj
}



DimPlot.metabolismV5<-function(obj,pathway,dimention.reduction.type="umap",dimention.reduction.run=T,
                               size=1)
{
  cat("\nPleaseCite:\nYingchengWu,QiangGao,etal.CancerDiscovery.2021.\nhttps://pubmed.ncbi.nlm.nih.gov/34417225/\n\n")
  if(dimention.reduction.type=="umap"){
    if(dimention.reduction.run==T)
      obj<-Seurat::RunUMAP(obj,reduction="pca",dims=1:40)
    umap.loc<-obj@reductions$umap@cell.embeddings
    row.names(umap.loc)<-colnames(obj)
    signature_exp<-obj@assays$METABOLISM$score
    input.pathway<-pathway
    signature_ggplot<-data.frame(umap.loc,t(signature_exp[input.pathway,
    ]))
    library(wesanderson)
    pal<-wes_palette("Zissou1",100,type="continuous")
    library(ggplot2)
    plot<-ggplot(data=signature_ggplot,aes(x=umap_1,
                                           y=umap_2,color=signature_ggplot[,3]))+geom_point(size=size)+
      scale_fill_gradientn(colours=pal)+scale_color_gradientn(colours=pal)+
      #labs(color=input.pathway)+xlab("UMAP1")+ylab("UMAP2")+theme_bw(base_size = 15)+
      theme(panel.border = element_rect(fill = NA, color = "black", size = 1, linetype = "solid"),
            panel.grid.major=element_blank(),
            panel.grid.minor=element_blank(),panel.background=element_blank(),
            legend.position = "right") +
      labs(title = input.pathway)
  }
  if(dimention.reduction.type=="tsne"){
    if(dimention.reduction.run==T)
      obj<-Seurat::RunTSNE(obj,reduction="pca",dims=1:40)
    tsne.loc<-obj@reductions$tsne@cell.embeddings
    row.names(tsne.loc)<-colnames(obj)
    signature_exp<-obj@assays$METABOLISM$score
    input.pathway<-pathway
    signature_ggplot<-data.frame(tsne.loc,t(signature_exp[input.pathway,
    ]))
    pal<-wes_palette("Zissou1",100,type="continuous")
    library(ggplot2)
    plot<-ggplot(data=signature_ggplot,aes(x=tSNE_1,
                                           y=tSNE_2,color=signature_ggplot[,3]))+geom_point(size=size)+
      scale_fill_gradientn(colours=pal)+scale_color_gradientn(colours=pal)+
      labs(color=input.pathway)+xlab("tSNE1")+ylab("tSNE2")+
      theme(aspect.ratio=1)+theme(panel.grid.major=element_blank(),
                                  panel.grid.minor=element_blank(),panel.background=element_blank(),
                                  axis.line=element_line(colour="black"))
  }
  plot
}

# Run metabolism analysis
metabolism_res <- sc.metabolism.SeuratV5(
  obj = sce_down,
  method = "ssGSEA",
  imputation = FALSE,
  ncores = 4,
  metabolism.type = "KEGG")

# Visualize metabolic pathways
metabolic_pathways <- c("Glycolysis / Gluconeogenesis",
                        "Citrate cycle (TCA cycle)",
                        "Oxidative phosphorylation",
                        "Glutathione metabolism")

lapply(metabolic_pathways, function(pathway) {
  DimPlot.metabolismV5(obj = metabolism_res,
                       pathway = pathway,
                       dimention.reduction.type = "umap",
                       dimention.reduction.run = FALSE,
                       size = 1)
})