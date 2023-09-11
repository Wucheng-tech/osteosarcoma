#####GSE162454 OS Untreated
library(Seurat)
setwd("./GSE162454/GSE162454_RAW/")
p0<- CreateSeuratObject(Read10X(data.dir = "./OS1"), project = "OS1_Tumor",min.features =300, min.cells=10)
p1 <- CreateSeuratObject(Read10X(data.dir = "./OS2"), project = "OS2_Tumor",min.features =300, min.cells=10)
p2 <- CreateSeuratObject(Read10X(data.dir = "./OS3"), project = "OS3_Tumor",min.features =300, min.cells=10)
p3 <- CreateSeuratObject(Read10X(data.dir = "./OS4"), project = "OS4_Tumor",min.features =300, min.cells=10)
p4 <- CreateSeuratObject(Read10X(data.dir = "./OS5"), project = "OS5_Tumor",min.features =300, min.cells=10)
p5 <- CreateSeuratObject(Read10X(data.dir = "./OS6"), project = "OS6_Tumor",min.features =300, min.cells=10)
ifnb.list<- list(p0,p1,p2,p3,p4,p5)
ifnb.list <- lapply(X = ifnb.list, FUN = function(x) {
    x <- NormalizeData(x, normalization.method = "LogNormalize", scale.factor = 10000)
})
OS <-merge(p0,y=c(p1,p2,p3,p4,p5))
OS_GSE162454<-OS
OS_GSE162454$donor_ID = substr(OS_GSE162454@meta.data[,1],1,4)
OS_GSE162454$donor_status = "Untreated"
OS_GSE162454$GSE_ID = "GSE162454"
saveRDS(OS_GSE162454, file="data_GSE162454.rds")
#
######GSE152048 OS Chemo
library(Seurat)
setwd("./GSE152048/primary")
p0<- CreateSeuratObject(Read10X(data.dir = "./BC2"), project = "BC2_Chemo",min.features =300, min.cells=10)
p1 <- CreateSeuratObject(Read10X(data.dir = "./BC3"), project = "BC3_Chemo",min.features =300, min.cells=10)
p2 <- CreateSeuratObject(Read10X(data.dir = "./BC5"), project = "BC5_Chemo",min.features =300, min.cells=10)
p3 <- CreateSeuratObject(Read10X(data.dir = "./BC6"), project = "BC6_Chemo",min.features =300, min.cells=10)
p4 <- CreateSeuratObject(Read10X(data.dir = "./BC16"), project = "BC16_Chemo",min.features =300, min.cells=10)
p5 <- CreateSeuratObject(Read10X(data.dir = "./BC21"), project = "BC21_Chemo",min.features =300, min.cells=10)
p6 <- CreateSeuratObject(Read10X(data.dir = "./BC22"), project = "BC22_Chemo",min.features =300, min.cells=10)
ifnb.list<- list(p0,p1,p2,p3,p4,p5,p6)
ifnb.list <- lapply(X = ifnb.list, FUN = function(x) {
    x <- NormalizeData(x, normalization.method = "LogNormalize", scale.factor = 10000)
})
OS <-merge(p0,y=c(p1,p2,p3,p4,p5,p6))
OS_GSE162454<-OS
OS_GSE162454$donor_ID = substr(OS_GSE162454@meta.data[,1],1,4)
OS_GSE162454$donor_status = "Chemo"
OS_GSE162454$GSE_ID = "GSE152048"
saveRDS(OS_GSE162454, file="data_GSE162454_pri.rds")
##Recurrent
library(Seurat)
setwd("./GSE152048/recurret")
p0<- CreateSeuratObject(Read10X(data.dir = "./BC11"), project = "BC11_recurrent",min.features =300, min.cells=10)
p1 <- CreateSeuratObject(Read10X(data.dir = "./BC20"), project = "BC20_recurrent",min.features =300, min.cells=10)
ifnb.list<- list(p0,p1)
ifnb.list <- lapply(X = ifnb.list, FUN = function(x) {
    x <- NormalizeData(x, normalization.method = "LogNormalize", scale.factor = 10000)
})
OS <-merge(p0,y=c(p1))
OS_GSE162454<-OS
OS_GSE162454$donor_ID = substr(OS_GSE162454@meta.data[,1],1,4)
OS_GSE162454$donor_status = "rechem"
OS_GSE162454$GSE_ID = "GSE152048"
saveRDS(OS_GSE162454, file="data_GSE162454_re.rds")

##Lung Metastasis
library(Seurat)
setwd("./GSE152048/Metastasis")
p0<- CreateSeuratObject(Read10X(data.dir = "./BC10"), project = "BC10_metastasis",min.features =300, min.cells=10)
p1 <- CreateSeuratObject(Read10X(data.dir = "./BC17"), project = "BC17_metastasis",min.features =300, min.cells=10)
ifnb.list<- list(p0,p1)
ifnb.list <- lapply(X = ifnb.list, FUN = function(x) {
    x <- NormalizeData(x, normalization.method = "LogNormalize", scale.factor = 10000)
})
OS <-merge(p0,y=c(p1))
OS_GSE162454<-OS
OS_GSE162454$donor_ID = substr(OS_GSE162454@meta.data[,1],1,4)
OS_GSE162454$donor_status = "Metastasis"
OS_GSE162454$GSE_ID = "GSE152048"
saveRDS(OS_GSE162454, file="data_GSE162454_Me.rds")

##########
######PRCA merge
####
library(Seurat)
data_GSE162454<-readRDS("./GSE162454/GSE162454_RAW/data_GSE162454.rds")
data_GSE162454_pri<-readRDS("./GSE152048/primary/data_GSE162454_pri.rds")
data_GSE162454_re<-readRDS("./GSE152048/recurret/data_GSE162454_re.rds")
data_GSE162454_Me<-readRDS("./GSE152048/Metastasis/data_GSE162454_Me.rds")
ifnb.list<- list(data_GSE162454,data_GSE162454_pri,data_GSE162454_re,data_GSE162454_Me)
ifnb.list <- lapply(X = ifnb.list, FUN = function(x) {
    x[["percent.mt"]]<- PercentageFeatureSet(x, pattern = "^(M|m)(T|t)-")
    x[["percent.RP"]]<- PercentageFeatureSet(x, pattern = "^RPL|^RPS")
    x <- subset(x, subset = nFeature_RNA > 500  & nFeature_RNA < 4000 & percent.mt < 50 & percent.RP < 50 )
})
ifnb.list <- lapply(X = ifnb.list, FUN = function(x) {
    x <- NormalizeData(x, normalization.method = "LogNormalize", scale.factor = 10000)
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})
features <- SelectIntegrationFeatures(object.list = ifnb.list)
ifnb.list <- lapply(X = ifnb.list, FUN = function(x) {
    x <- ScaleData(x, features = features, verbose = FALSE)
    x <- RunPCA(x, features = features, verbose = FALSE)
})##
anchors <- FindIntegrationAnchors(object.list = ifnb.list, anchor.features = features, reduction = "rpca")
combined <- IntegrateData(anchorset = anchors)
DefaultAssay(combined) <- "integrated"
# Run the standard workflow for visualization and clustering
combined <- ScaleData(combined, verbose = FALSE)
combined <- RunPCA(combined, npcs = 30, verbose = FALSE)
combined <- RunUMAP(combined, reduction = "pca", dims = 1:30)
combined <- FindNeighbors(combined, reduction = "pca", dims = 1:30)
combined <- FindClusters(combined, resolution = 0.3)
setwd("./Analysis")
p0 <- DimPlot(object = combined, reduction = "umap",label=TRUE,raster=FALSE)
p1 <- DimPlot(object = combined, reduction = "umap",group.by ="donor_ID",label=TRUE,raster=FALSE)
p2 <- DimPlot(object = combined, reduction = "umap",group.by ="donor_status",label=TRUE,raster=FALSE)
pdf("Umap_Sample_combined_Nor1.pdf",width=22,height=20)
CombinePlots(plots  = list(p0, p1, p2))
dev.off()
saveRDS(combined, file="data_rPCA_combined.rds")
############
#combined <- readRDS("data_rPCA_combined.rds")
OS <- combined
DefaultAssay(OS) <- "RNA"
OS.markers <- FindAllMarkers(object = OS, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.table(OS.markers,"markers_r0.3.txt",row.names=TRUE,col.names=TRUE)
write.table(OS@meta.data,"metadata_r0.3.txt",row.names=TRUE,col.names=TRUE,sep="\t")

#OS <- NormalizeData(object = OS, normalization.method = "LogNormalize", scale.factor = 10000)
#saveRDS(OS, file="OS1.rds")
###25029 features across 119283 samples within 2 assays 
######
#########
########识别细胞类型
# COL2A1,SOX9 ACAN       Chondroblastic
# CTSK,CKB,CSTB          Osteoclast
# CD3D NKG7 CD3E         immune  TIL 
# CD74  LYZ              Myeloid
# C1QA C1QB APOC1 FCGR3A Myeloid
# SFRP2  CXCL12               MSC
# VWF PECAM1            Endothelial
# CD79A CD79B           B
# ACTA2                  CAF
#######
library(Seurat)
setwd("./Analysis/Gene_exp")
pdf("umap_T.pdf",width=20,height=10)
print(FeaturePlot(OS, features =c("CD3D","NKG7"), slot = "data",cols = c("grey", "red"),blend.threshold = 1,reduction = "umap",raster=FALSE))
dev.off()
pdf("umap_Mye.pdf",width=20,height=10)
print(FeaturePlot(OS, features =c("CD74","FCGR3A"), slot = "data",cols = c("grey", "red"),blend.threshold = 1,reduction = "umap",raster=FALSE))
dev.off()
print(FeaturePlot(OS, features =c("COL2A1","SOX9"), slot = "data",cols = c("grey", "red"),blend.threshold = 1,reduction = "umap",raster=FALSE))
dev.off()
pdf("umap_Osteoclast.pdf",width=20,height=10)
print(FeaturePlot(OS, features =c("CTSK","CKB"), slot = "data",cols = c("grey", "red"),blend.threshold = 1,reduction = "umap",raster=FALSE))
dev.off()
pdf("umap_Endothelial.pdf",width=20,height=10)
print(FeaturePlot(OS, features =c("VWF","PECAM1"), slot = "data",cols = c("grey", "red"),blend.threshold = 1,reduction = "umap",raster=FALSE))
dev.off()
pdf("umap_B.pdf",width=20,height=10)
print(FeaturePlot(OS, features =c("CD79A","CD79B"), slot = "data",cols = c("grey", "red"),blend.threshold = 1,reduction = "umap",raster=FALSE))
dev.off()
pdf("umap_MSC.pdf",width=20,height=10)
print(FeaturePlot(OS, features =c("SFRP2","CXCL12"), slot = "data",cols = c("grey", "red"),blend.threshold = 1,reduction = "umap",raster=FALSE))
dev.off()
pdf("umap_CAF.pdf",width=20,height=10)
print(FeaturePlot(OS, features =c("ACTA2","CXCL12"), slot = "data",cols = c("grey", "red"),blend.threshold = 1,reduction = "umap",raster=FALSE))
dev.off()
df("Umap_Chondroblastic.pdf",width=20,height=10)
print(FeaturePlot(OS, features =c("COL2A1","SOX9"), slot = "data",cols = c("grey", "red"),blend.threshold = 1,reduction = "umap",raster=FALSE))
dev.off()
pdf("Umap_Osteoblasts.pdf",width=20,height=10)
print(FeaturePlot(OS, features =c("ALPL","RUNX2"), slot = "data",cols = c("grey", "red"),blend.threshold = 1,reduction = "umap",raster=FALSE))
dev.off()
null device 
#
#####
#########
#####
gene <-c("COL1A1","CDH11","RUNX2","ALPL","LUM","COL2A1","SOX9","ACAN","ISG15","MMP13","ACTA2","RGS5","SFRP2","CXCL12","MME","CTSK","CKB","CSTB","CD3D","NKG7","CD3E","CD74","LYZ","C1QA","C1QB","CD79A","CD79B","VWF","PECAM1")
library(Seurat)
library(patchwork)
DefaultAssay(OS) <- "RNA"
#OS$seurat_clusters <- factor(OS$seurat_clusters, levels =c("4","11","17","13","6","10","19","14","16","21","25","2","9","20","22","23","7","5","15","3","26","27","0","1","18","8","24","12"))
OS$seurat_clusters <- factor(OS$seurat_clusters, levels =c("1","2","6","7","8","13","15","4","9","16","17","20","14","5","3","19","0","11","12","22","18","21","10"))
plots <- VlnPlot(OS, features =gene,  group.by = "seurat_clusters", pt.size = 0, combine = FALSE)
pdf("Gene1.pdf",width=20,height=50)
wrap_plots(plots = plots, ncol = 1)
dev.off() 
###
######
###
gene <-c("COL1A1","CDH11","RUNX2","ALPL","LUM","COL2A1","SOX9","ACAN","ISG15","MMP13","ACTA2","RGS5","SFRP2","CXCL12","MME","CTSK","CKB","CSTB","CD3D","NKG7","CD3E","CD74","LYZ","C1QA","C1QB","CD79A","CD79B","VWF","PECAM1")
OS_data1 <-OS[["RNA"]]@data[gene,]
OS_data <-t(scale(t(as.matrix(OS_data1))))
type <-unique(OS@meta.data[,"seurat_clusters"])
bb <-NULL
for(i in 1: length(type)){
ind <-which(OS@meta.data[,"seurat_clusters"]==type[i])
aa <-rowMeans(as.matrix(OS_data[,ind]))
bb <-cbind(bb,aa)
}
colnames(bb) <-type
library(reshape2)
library(ggplot2)
bb[bb>4] <-3.8
data_melt<-melt(bb)
data_melt$Var2 <-factor(data_melt$Var2,levels=c("1","2","6","7","13","15","8","4","16","9","14","17","20","5","3","19","0","12","22","11","18","21","10"))
p <- ggplot(data = data_melt) + geom_tile(aes(x=Var2,y=Var1, fill = value)) +theme_classic() +  theme(axis.ticks = element_blank(), axis.line = element_blank()) + xlab('row name') + ylab('column name')
p1 <- p+scale_fill_gradient2('z-score',low = 'blue', high = 'red', mid = 'white') 
pdf("heatmap_gene_celltype.pdf",width=10,height=15)			
p1
dev.off()
###附图fig s1 B
marker <-read.table("./Analysis/Gene_exp/markers_r0.3.txt")
gene1 <-NULL
for(i in c("1","2","6","7","13","15","8","4","16","9","14","17","20","5","3","19","0","12","22","11","18","21","10")){
index <- which(marker$cluster==i)
gene2 <-marker[index,"gene"][1:10]
gene1 <-c(gene1,gene2)
}
gene <-unique(gene1)
gene <- gene[!is.na(gene)] 
gene <-gene[-grep("^RPS|RPL|MT",gene)]
OS_data1 <-OS[["RNA"]]@data[gene,]
OS_data <-t(scale(t(as.matrix(OS_data1))))
type <-unique(OS@meta.data[,"seurat_clusters"])
bb <-NULL
for(i in 1: length(type)){
ind <-which(OS@meta.data[,"seurat_clusters"]==type[i])
aa <-rowMeans(as.matrix(OS_data[,ind]))
bb <-cbind(bb,aa)
}
colnames(bb) <-type
library(reshape2)
library(ggplot2)
bb[bb>4] <-3.8
data_melt<-melt(bb)
data_melt$Var2 <-factor(data_melt$Var2,levels=c("1","2","6","7","13","15","8","4","16","9","14","17","20","5","3","19","0","12","22","11","18","21","10"))
p <- ggplot(data = data_melt) + geom_tile(aes(x=Var2,y=Var1, fill = value)) +theme_classic() +  theme(axis.ticks = element_blank(), axis.line = element_blank()) + xlab('row name') + ylab('column name')
p1 <- p+scale_fill_gradient2('z-score',low = 'blue', high = 'red', mid = 'white') 
pdf("heatmap_gene_celltype1.pdf",width=10,height=25)			
p1
dev.off()
####fig 1 B
current.cluster.ids <- c("1","2","6","7","13","15","8","4","16","9","14","17","20","5","3","19","0","12","22","11","18","21","10")
new.cluster.ids <- c("Osteoblasts","Osteoblasts","Osteoblasts","Osteoblasts","Osteoblasts","Osteoblasts","Chondroblastic",
"Fibroblast","Fibroblast","Pericyte","MSC","MSC","MSC","Osteoclast","TIL","TIL","Myeloid","Myeloid","Myeloid","Myeloid","B cell","B cell","Endothelial cell")
OS@meta.data$celltype <- plyr::mapvalues(x = OS@meta.data[,"seurat_clusters"], from = current.cluster.ids, to = new.cluster.ids)

setwd("./Analysis/")
UmapPlot<- function(obj) {
  library(ggplot2)
  library(dplyr)
  UMAP = obj@reductions$umap@cell.embeddings %>% as.data.frame() %>% cbind(cell_type = obj@meta.data$seurat_clusters) 
  allcolour=c("#DC143C","#0000FF","#20B2AA","#FFA500","#9370DB","#98FB98","#F08080","#1E90FF","#7CFC00","#FFFF00","#808000","#FF00FF","#FA8072","#7B68EE","#9400D3","#800080","#A0522D","#D2B48C","#D2691E","#87CEEB","#40E0D0","#5F9EA0","#FF1493","#0000CD","#008B8B","#FFE4B5","#8A2BE2","#228B22","#E9967A","#4682B4","#32CD32","#F0E68C","#FFFFE0","#EE82EE","#FF6347","#6A5ACD","#9932CC","#8B008B","#8B4513","#DEB887")
  my36colors <-c('#E5D2DD', '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87','#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#8C549C', '#585658', '#9FA3A8', '#E0D4CA', '#5F3D69', '#C5DEBA', '#58A4C3', '#E4C755', '#F7F398','#AA9A59', '#E63863', '#E39A35', '#C1E6F3', '#6778AE', '#91D0BE', '#B53E2B','#712820', '#DCC1DD', '#CCE0F5',  '#CCC9E6', '#625D9E', '#68A180', '#3A6963','#968175')
  p <- ggplot(UMAP,aes(x= UMAP_1, y = UMAP_2 ,color = cell_type)) +  geom_point(size=0.2, alpha =1 )  +  scale_color_manual(values = my36colors)
  p2 <- p  +theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),  panel.border = element_blank(),  axis.title = element_blank(),  axis.text = element_blank(), axis.ticks = element_blank(),panel.background = element_rect(fill = 'white'), plot.background=element_rect(fill="white"))
  p3 <- p2 + theme(legend.title = element_blank(),legend.key=element_rect(fill='white'), legend.text = element_text(size=20),legend.key.size=unit(1,'cm') ) + guides(color = guide_legend(override.aes = list(size=5))) 
#  p4 <- p3 + geom_segment(aes(x = min(TSNE$tSNE_1) , y = min(TSNE$tSNE_2) ,xend = min(TSNE$tSNE_1) +3, yend = min(TSNE$tSNE_2) ),colour = "black", size=1,arrow = arrow(length = unit(0.3,"cm")))+ geom_segment(aes(x = min(TSNE$tSNE_1), y = min(TSNE$tSNE_2),xend = min(TSNE$tSNE_1) , yend = min(TSNE$tSNE_2) + 3),colour = "black", size=1,arrow = arrow(length = unit(0.3,"cm"))) + annotate("text", x = min(TSNE$tSNE_1) +1.5, y = min(TSNE$tSNE_2) -1, label = "tSNE_1",color="black",size = 3, fontface="bold" ) + annotate("text", x = min(TSNE$tSNE_1) -1, y = min(TSNE$tSNE_2) + 1.5, label = "tSNE_2",color="black",size = 3, fontface="bold" ,angle=90) 
  cell_type_med <- UMAP %>%group_by(cell_type) %>%summarise(UMAP_1 = median(UMAP_1),UMAP_2 = median(UMAP_2))
  library(ggrepel)
  p5 <- p3 + geom_label_repel(aes(label=cell_type), fontface="bold",data = cell_type_med, point.padding=unit(1.5, "lines")) + theme(legend.position = "none")                  
  return(p5)
}
pdf("UMAP_Seurat.pdf",width=5,height=5)
UmapPlot(OS)
dev.off()

##
OS@meta.data$celltype <-factor(OS@meta.data$celltype,levels=c("Osteoblasts","Chondroblastic","Fibroblast","Pericyte","MSC","Osteoclast","TIL","Myeloid","B cell","Endothelial cell"))
UmapPlot<- function(obj) {
  library(ggplot2)
  library(dplyr)
  UMAP = obj@reductions$umap@cell.embeddings %>% as.data.frame() %>% cbind(cell_type = obj@meta.data$celltype) 
  allcolour=c("#DC143C","#0000FF","#20B2AA","#FFA500","#9370DB","#98FB98","#F08080","#1E90FF","#7CFC00","#FFFF00","#808000","#FF00FF","#FA8072","#7B68EE","#9400D3","#800080","#A0522D","#D2B48C","#D2691E","#87CEEB","#40E0D0","#5F9EA0","#FF1493","#0000CD","#008B8B","#FFE4B5","#8A2BE2","#228B22","#E9967A","#4682B4","#32CD32","#F0E68C","#FFFFE0","#EE82EE","#FF6347","#6A5ACD","#9932CC","#8B008B","#8B4513","#DEB887")
  my36colors <-c('#E5D2DD', '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87','#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#8C549C', '#585658', '#9FA3A8', '#E0D4CA', '#5F3D69', '#C5DEBA', '#58A4C3', '#E4C755', '#F7F398','#AA9A59', '#E63863', '#E39A35', '#C1E6F3', '#6778AE', '#91D0BE', '#B53E2B','#712820', '#DCC1DD', '#CCE0F5',  '#CCC9E6', '#625D9E', '#68A180', '#3A6963','#968175')
  UMAP$cell_type <-factor(UMAP$cell_type,levels=c("Osteoblasts","Chondroblastic","Fibroblast","Pericyte","MSC","Osteoclast","TIL","Myeloid","B cell","Endothelial cell"))
  p <- ggplot(UMAP,aes(x= UMAP_1, y = UMAP_2 ,color = cell_type)) +  geom_point(size=0.2, alpha =1 )  +  scale_color_manual(values = my36colors[1:10])
  p2 <- p  +theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),  panel.border = element_blank(),  axis.title = element_blank(),  axis.text = element_blank(), axis.ticks = element_blank(),panel.background = element_rect(fill = 'white'), plot.background=element_rect(fill="white"))
  p3 <- p2 + theme(legend.title = element_blank(),legend.key=element_rect(fill='white'), legend.text = element_text(size=20),legend.key.size=unit(1,'cm') ) + guides(color = guide_legend(override.aes = list(size=5))) 
#  p4 <- p3 + geom_segment(aes(x = min(TSNE$tSNE_1) , y = min(TSNE$tSNE_2) ,xend = min(TSNE$tSNE_1) +3, yend = min(TSNE$tSNE_2) ),colour = "black", size=1,arrow = arrow(length = unit(0.3,"cm")))+ geom_segment(aes(x = min(TSNE$tSNE_1), y = min(TSNE$tSNE_2),xend = min(TSNE$tSNE_1) , yend = min(TSNE$tSNE_2) + 3),colour = "black", size=1,arrow = arrow(length = unit(0.3,"cm"))) + annotate("text", x = min(TSNE$tSNE_1) +1.5, y = min(TSNE$tSNE_2) -1, label = "tSNE_1",color="black",size = 3, fontface="bold" ) + annotate("text", x = min(TSNE$tSNE_1) -1, y = min(TSNE$tSNE_2) + 1.5, label = "tSNE_2",color="black",size = 3, fontface="bold" ,angle=90) 
  cell_type_med <- UMAP %>%group_by(cell_type) %>%summarise(UMAP_1 = median(UMAP_1),UMAP_2 = median(UMAP_2))
  library(ggrepel)
  p5 <- p3 + geom_label_repel(aes(label=cell_type), fontface="bold",data = cell_type_med, point.padding=unit(1.5, "lines")) + theme(legend.position = "none")                  
  return(p5)
}
pdf("UMAP_celltype.pdf",width=5,height=5)
UmapPlot(OS)
dev.off()
###fig 1s A
current.cluster.ids <- c("Untreated","Chemo","Rechem","Metastasis")
new.cluster.ids <- c("Tumor","Chemo","Recurrent","Metastasis")
OS@meta.data$donor_status1 <- plyr::mapvalues(x = OS@meta.data[,"donor_status"], from = current.cluster.ids, to = new.cluster.ids)
OS@meta.data$donor_status1 <-factor(OS@meta.data$donor_status1,levels=c("Tumor","Chemo","Recurrent","Metastasis"))
UmapPlot<- function(obj) {
  library(ggplot2)
  library(dplyr)
  UMAP = obj@reductions$umap@cell.embeddings %>% as.data.frame() %>% cbind(cell_type = obj@meta.data$donor_status1) 
  allcolour=c("#DC143C","#0000FF","#20B2AA","#FFA500","#9370DB","#98FB98","#F08080","#1E90FF","#7CFC00","#FFFF00","#808000","#FF00FF","#FA8072","#7B68EE","#9400D3","#800080","#A0522D","#D2B48C","#D2691E","#87CEEB","#40E0D0","#5F9EA0","#FF1493","#0000CD","#008B8B","#FFE4B5","#8A2BE2","#228B22","#E9967A","#4682B4","#32CD32","#F0E68C","#FFFFE0","#EE82EE","#FF6347","#6A5ACD","#9932CC","#8B008B","#8B4513","#DEB887")
  my36colors <-c('#E5D2DD', '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87','#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#8C549C', '#585658', '#9FA3A8', '#E0D4CA', '#5F3D69', '#C5DEBA', '#58A4C3', '#E4C755', '#F7F398','#AA9A59', '#E63863', '#E39A35', '#C1E6F3', '#6778AE', '#91D0BE', '#B53E2B','#712820', '#DCC1DD', '#CCE0F5',  '#CCC9E6', '#625D9E', '#68A180', '#3A6963','#968175')
  p <- ggplot(UMAP,aes(x= UMAP_1, y = UMAP_2 ,color = cell_type)) +  geom_point(size=0.2, alpha =1 )  +  scale_color_manual(values = my36colors[1:10])
  p2 <- p  +theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),  panel.border = element_blank(),  axis.title = element_blank(),  axis.text = element_blank(), axis.ticks = element_blank(),panel.background = element_rect(fill = 'white'), plot.background=element_rect(fill="white"))
  p3 <- p2 + theme(legend.title = element_blank(),legend.key=element_rect(fill='white'), legend.text = element_text(size=20),legend.key.size=unit(1,'cm') ) + guides(color = guide_legend(override.aes = list(size=5))) 
#  p4 <- p3 + geom_segment(aes(x = min(TSNE$tSNE_1) , y = min(TSNE$tSNE_2) ,xend = min(TSNE$tSNE_1) +3, yend = min(TSNE$tSNE_2) ),colour = "black", size=1,arrow = arrow(length = unit(0.3,"cm")))+ geom_segment(aes(x = min(TSNE$tSNE_1), y = min(TSNE$tSNE_2),xend = min(TSNE$tSNE_1) , yend = min(TSNE$tSNE_2) + 3),colour = "black", size=1,arrow = arrow(length = unit(0.3,"cm"))) + annotate("text", x = min(TSNE$tSNE_1) +1.5, y = min(TSNE$tSNE_2) -1, label = "tSNE_1",color="black",size = 3, fontface="bold" ) + annotate("text", x = min(TSNE$tSNE_1) -1, y = min(TSNE$tSNE_2) + 1.5, label = "tSNE_2",color="black",size = 3, fontface="bold" ,angle=90) 
  cell_type_med <- UMAP %>%group_by(cell_type) %>%summarise(UMAP_1 = median(UMAP_1),UMAP_2 = median(UMAP_2))
  library(ggrepel)
  p5 <- p3 + geom_label_repel(aes(label=cell_type), fontface="bold",data = cell_type_med, point.padding=unit(1.5, "lines")) + theme(legend.position = "none")                  
  return(p5)
}
pdf("UMAP_condition.pdf",width=5,height=5)
UmapPlot(OS)
dev.off()

###比例 fig 1D
setwd("./Analysis")
bb <-NULL
sample <- unique(OS@meta.data[,"donor_status"])
for(i in 1:length(sample)){
aa <-table(OS@meta.data[which(OS@meta.data[,"donor_status"]==sample[i]),"celltype"])
bb <-cbind(bb,aa)
}
colnames(bb) <-sample
cc <-cbind(bb[,1]/sum(bb[,1]),bb[,2]/sum(bb[,2]),bb[,3]/sum(bb[,3]),bb[,4]/sum(bb[,4]))
colnames(cc) <-colnames(bb)
library(reshape2)
library(ggplot2)
library(RColorBrewer)
colorr <- brewer.pal(4,"Set2")
df <-melt(cc)
df$Var2 <-factor(df$Var2,levels=c("Metastasis","rechem","Chemo","Tumor"))
df$Var1 <-factor(df$Var1,levels=c("Osteoblasts","Chondroblastic","Fibroblast","Pericyte","MSC","Osteoclast","TIL","Myeloid","B cell","Endothelial cell"))
pdf(file="box_precent_new.pdf",width=5,height=10)
my36colors <-c('#E5D2DD', '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87','#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#8C549C', '#585658', '#9FA3A8', '#E0D4CA', '#5F3D69', '#C5DEBA', '#58A4C3', '#E4C755', '#F7F398','#AA9A59', '#E63863', '#E39A35', '#C1E6F3', '#6778AE', '#91D0BE', '#B53E2B','#712820', '#DCC1DD', '#CCE0F5',  '#CCC9E6', '#625D9E', '#68A180', '#3A6963','#968175')
ggplot(df,aes(x=Var2,y=value,fill=Var1))+scale_fill_manual(values=my36colors[1:10])+geom_bar(stat="identity")+theme(legend.title=element_blank()) 
dev.off()
######## fig s1 C
###比例
bb <-NULL
sample <- unique(OS@meta.data[,"donor_ID"])
for(i in 1:length(sample)){
aa <-table(OS@meta.data[which(OS@meta.data[,"donor_ID"]==sample[i]),"celltype"])
bb <-cbind(bb,aa)
}
colnames(bb) <-sample
cc <-NULL
for(i in 1:length(sample)){
aa <-bb[,i]/sum(bb[,i])
cc <-cbind(cc,aa)
}
colnames(cc) <-colnames(bb)
library(reshape2)
library(ggplot2)
library(RColorBrewer)
colorr <- brewer.pal(4,"Set2")
df <-melt(cc)
df$Var2 <-factor(df$Var2,levels=rev(sample))
df$Var1 <-factor(df$Var1,levels=c("Osteoblasts","Chondroblastic","Fibroblast","Pericyte","MSC","Osteoclast","TIL","Myeloid","B cell","Endothelial cell"))
pdf(file="box_precent_new1.pdf",width=20,height=10)
ggplot(df,aes(x=Var2,y=value,fill=Var1))+scale_fill_manual(values=my36colors[1:10])+geom_bar(stat="identity")+theme(legend.title=element_blank()) 
dev.off()
########饼图
pdf(file="Pie_precent.pdf",width=20,height=10)
pie(table(OS@meta.data$celltype)[c("Osteoblasts","Chondroblastic","Fibroblast","Pericyte","MSC","Osteoclast","TIL","Myeloid","B cell","Endothelial cell")])
dev.off()
########box图
current.cluster.ids <- c("Untreated","Chemo","rechem","Metastasis")
new.cluster.ids <- c("Tumor","Chemo","Re_Me","Re_Me")
OS@meta.data$donor_status2 <- plyr::mapvalues(x = OS@meta.data[,"donor_status"], from = current.cluster.ids, to = new.cluster.ids)
library(reshape2)
library(ggplot2)
library(RColorBrewer)
dat <- OS@meta.data[,c("donor_status2","celltype","donor_ID")]
dat1 <- dat[which(dat$donor_status2 =="Tumor"),]
bb <-NULL
sample <- unique(dat1[,"donor_ID"])
for(i in 1:length(sample)){
aa <-table(dat1[which(dat1[,"donor_ID"]==sample[i]),"celltype"])
bb <-cbind(bb,aa)}
colnames(bb) <-sample
cc <-t(t(bb)/as.vector(table(OS$donor_ID)[sample]))
dd1 <-melt(cc)
dd1$status <-"Tumor"
dat1 <- dat[which(dat$donor_status2 =="Chemo"),]
bb <-NULL
sample <- unique(dat1[,"donor_ID"])
for(i in 1:length(sample)){
aa <-table(dat1[which(dat1[,"donor_ID"]==sample[i]),"celltype"])
bb <-cbind(bb,aa)}
colnames(bb) <-sample
cc <-t(t(bb)/as.vector(table(OS$donor_ID)[sample]))
dd2 <-melt(cc)
dd2$status <-"Chemo"
#
dat1 <- dat[which(dat$donor_status2 =="Re_Me"),]
bb <-NULL
sample <- unique(dat1[,"donor_ID"])
for(i in 1:length(sample)){
aa <-table(dat1[which(dat1[,"donor_ID"]==sample[i]),"celltype"])
bb <-cbind(bb,aa)}
colnames(bb) <-sample
cc <-t(t(bb)/as.vector(table(OS$donor_ID)[sample]))
dd3 <-melt(cc)
dd3$status <-"Re_Me"
df <-rbind(dd1,dd2,dd3)
levels <-c("Osteoblasts","Chondroblastic","Fibroblast","Pericyte","MSC","Osteoclast","TIL","Myeloid","B cell","Endothelial cell")
df$status = factor(df$status, levels=c("Tumor","Chemo","Re_Me"))
dodge <- position_dodge(width = 0.8)
pdf(file="box_precent_stage.pdf",width=10,height=4)
ggplot(df, aes(x=factor(Var1,levels=levels), y=value,fill = status))  + stat_boxplot(geom="errorbar",width=0.3,position = dodge,size=2)+geom_boxplot(width=0.8,outlier.colour = NA)+scale_fill_manual(values=my36colors[1:15])+
ylim(0,0.8)+theme(legend.title=element_blank())+theme(axis.text.x = element_text(size = 12,angle = 45, hjust = 1, vjust = 1, face = "bold"))+labs(x="Cluster", y = "overlap")
dev.off()
##
current.cluster.ids <- c("Untreated","Chemo","rechem","Metastasis")
new.cluster.ids <- c("Tumor","Chemo","Recurrent","Metastasis")
OS@meta.data$donor_status3 <- plyr::mapvalues(x = OS@meta.data[,"donor_status"], from = current.cluster.ids, to = new.cluster.ids)
library(reshape2)
library(ggplot2)
library(RColorBrewer)
dat <- OS@meta.data[,c("donor_status3","celltype","donor_ID")]
dat1 <- dat[which(dat$donor_status3 =="Tumor"),]
bb <-NULL
sample <- unique(dat1[,"donor_ID"])
for(i in 1:length(sample)){
aa <-table(dat1[which(dat1[,"donor_ID"]==sample[i]),"celltype"])
bb <-cbind(bb,aa)}
colnames(bb) <-sample
cc <-t(t(bb)/as.vector(table(OS$donor_ID)[sample]))
dd1 <-melt(cc)
dd1$status <-"Tumor"
dat1 <- dat[which(dat$donor_status3 =="Chemo"),]
bb <-NULL
sample <- unique(dat1[,"donor_ID"])
for(i in 1:length(sample)){
aa <-table(dat1[which(dat1[,"donor_ID"]==sample[i]),"celltype"])
bb <-cbind(bb,aa)}
colnames(bb) <-sample
cc <-t(t(bb)/as.vector(table(OS$donor_ID)[sample]))
dd2 <-melt(cc)
dd2$status <-"Chemo"
#
dat1 <- dat[which(dat$donor_status3 =="Recurrent"),]
bb <-NULL
sample <- unique(dat1[,"donor_ID"])
for(i in 1:length(sample)){
aa <-table(dat1[which(dat1[,"donor_ID"]==sample[i]),"celltype"])
bb <-cbind(bb,aa)}
colnames(bb) <-sample
cc <-t(t(bb)/as.vector(table(OS$donor_ID)[sample]))
dd3 <-melt(cc)
dd3$status <-"Recurrent"

dat1 <- dat[which(dat$donor_status3 =="Metastasis"),]
bb <-NULL
sample <- unique(dat1[,"donor_ID"])
for(i in 1:length(sample)){
aa <-table(dat1[which(dat1[,"donor_ID"]==sample[i]),"celltype"])
bb <-cbind(bb,aa)}
colnames(bb) <-sample
cc <-t(t(bb)/as.vector(table(OS$donor_ID)[sample]))
dd4 <-melt(cc)
dd4$status <-"Metastasis"

df <-rbind(dd1,dd2,dd3,dd4)
levels <-c("Osteoblasts","Chondroblastic","Fibroblast","Pericyte","MSC","Osteoclast","TIL","Myeloid","B cell","Endothelial cell")
df$status = factor(df$status, levels=c("Tumor","Chemo","Recurrent","Metastasis"))
dodge <- position_dodge(width = 0.8)
pdf(file="box_precent_stage1.pdf",width=10,height=4)
ggplot(df, aes(x=factor(Var1,levels=levels), y=value,fill = status))  + stat_boxplot(geom="errorbar",width=0.3,position = dodge,size=2)+geom_boxplot(width=0.8,outlier.colour = NA)+scale_fill_manual(values=my36colors[1:15])+
ylim(0,0.8)+theme(legend.title=element_blank())+theme(axis.text.x = element_text(size = 12,angle = 45, hjust = 1, vjust = 1, face = "bold"))+labs(x="Cluster", y = "overlap")
dev.off()

#############t.test/wilcox
bb <-NULL
for(i in c("Osteoblasts","Chondroblastic","Fibroblast","Pericyte","MSC","Osteoclast","TIL","Myeloid","B cell","Endothelial cell")){
aa <- c(wilcox.test(df[which(df[,1]==i),][which(df[which(df[,1]==i),][,4]=="Tumor"),3],df[which(df[,1]==i),][which(df[which(df[,1]==i),][,4]=="Chemo"),3])$p.value,
wilcox.test(df[which(df[,1]==i),][which(df[which(df[,1]==i),][,4]=="Tumor"),3],df[which(df[,1]==i),][which(df[which(df[,1]==i),][,4]=="Re_Me"),3])$p.value,
wilcox.test(df[which(df[,1]==i),][which(df[which(df[,1]==i),][,4]=="Chemo"),3],df[which(df[,1]==i),][which(df[which(df[,1]==i),][,4]=="Re_Me"),3])$p.value)
bb <-rbind(bb,aa)
}
bb
         [,1]       [,2]        [,3]
aa 0.53379953 0.91428571 0.230303030
aa 0.28131745 0.60952381 0.107395172
aa 0.53379953 0.76190476 0.648484848
aa 0.73076923 0.91428571 0.787878788
aa 1.00000000 0.60952381 0.315151515
aa 0.29487179 0.23952554 0.104984533
aa 0.13752914 0.06666667 0.648484848
aa 0.44522145 0.06666667 0.006060606
aa 0.36596737 0.11428571 0.527272727
aa 0.01398601 0.25714286 0.012121212






#########
#############
##################
########### part 2
#####  Osteoblasts
library(Seurat)
setwd("./Analysis/Osteoblasts")
sub_OS <-subset(OS,idents= c("1","2","6","7","13","15"))
DefaultAssay(sub_OS) <- "integrated"
# Run the standard workflow for visualization and clustering
#sub_OS <- ScaleData(sub_OS, verbose = FALSE)
sub_OS <- RunPCA(sub_OS, npcs = 30, verbose = FALSE)
sub_OS <- RunUMAP(sub_OS, reduction = "pca", dims = 1:30)
sub_OS <- FindNeighbors(sub_OS, reduction = "pca", dims = 1:30)
sub_OS <- FindClusters(sub_OS, resolution = 0.1)
p0 <- DimPlot(object = sub_OS, reduction = "umap",label=TRUE)
p1<- DimPlot(object = sub_OS, reduction = "umap",group.by ="donor_ID",label=TRUE)
p2 <- DimPlot(object = sub_OS, reduction = "umap",group.by ="donor_status",label=TRUE)
pdf("Umap_sub_OS.pdf",width=22,height=20)
CombinePlots(plots  = list(p0, p1, p2))
dev.off()
saveRDS(sub_OS, file="data_rPCA_sub_OS.rds")
#sub_OS <-readRDS("data_rPCA_sub_OS.rds")
DefaultAssay(sub_OS) <- "RNA"
OS.markers <- FindAllMarkers(object = sub_OS, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.table(OS.markers,"markers_r0.1.txt",row.names=TRUE,col.names=TRUE)
#
UmapPlot<- function(obj) {
  library(ggplot2)
  library(dplyr)
  UMAP = obj@reductions$umap@cell.embeddings %>% as.data.frame() %>% cbind(cell_type = obj@meta.data$seurat_clusters) 
  allcolour=c("#DC143C","#0000FF","#20B2AA","#FFA500","#9370DB","#98FB98","#F08080","#1E90FF","#7CFC00","#FFFF00","#808000","#FF00FF","#FA8072","#7B68EE","#9400D3","#800080","#A0522D","#D2B48C","#D2691E","#87CEEB","#40E0D0","#5F9EA0","#FF1493","#0000CD","#008B8B","#FFE4B5","#8A2BE2","#228B22","#E9967A","#4682B4","#32CD32","#F0E68C","#FFFFE0","#EE82EE","#FF6347","#6A5ACD","#9932CC","#8B008B","#8B4513","#DEB887")
  my36colors <-c('#E5D2DD', '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87','#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#8C549C', '#585658', '#9FA3A8', '#E0D4CA', '#5F3D69', '#C5DEBA', '#58A4C3', '#E4C755', '#F7F398','#AA9A59', '#E63863', '#E39A35', '#C1E6F3', '#6778AE', '#91D0BE', '#B53E2B','#712820', '#DCC1DD', '#CCE0F5',  '#CCC9E6', '#625D9E', '#68A180', '#3A6963','#968175')
  p <- ggplot(UMAP,aes(x= UMAP_1, y = UMAP_2 ,color = cell_type)) +  geom_point(size=0.2, alpha =1 )  +  scale_color_manual(values = my36colors)
  p2 <- p  +theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),  panel.border = element_blank(),  axis.title = element_blank(),  axis.text = element_blank(), axis.ticks = element_blank(),panel.background = element_rect(fill = 'white'), plot.background=element_rect(fill="white"))
  p3 <- p2 + theme(legend.title = element_blank(),legend.key=element_rect(fill='white'), legend.text = element_text(size=20),legend.key.size=unit(1,'cm') ) + guides(color = guide_legend(override.aes = list(size=5))) 
#  p4 <- p3 + geom_segment(aes(x = min(TSNE$tSNE_1) , y = min(TSNE$tSNE_2) ,xend = min(TSNE$tSNE_1) +3, yend = min(TSNE$tSNE_2) ),colour = "black", size=1,arrow = arrow(length = unit(0.3,"cm")))+ geom_segment(aes(x = min(TSNE$tSNE_1), y = min(TSNE$tSNE_2),xend = min(TSNE$tSNE_1) , yend = min(TSNE$tSNE_2) + 3),colour = "black", size=1,arrow = arrow(length = unit(0.3,"cm"))) + annotate("text", x = min(TSNE$tSNE_1) +1.5, y = min(TSNE$tSNE_2) -1, label = "tSNE_1",color="black",size = 3, fontface="bold" ) + annotate("text", x = min(TSNE$tSNE_1) -1, y = min(TSNE$tSNE_2) + 1.5, label = "tSNE_2",color="black",size = 3, fontface="bold" ,angle=90) 
  cell_type_med <- UMAP %>%group_by(cell_type) %>%summarise(UMAP_1 = median(UMAP_1),UMAP_2 = median(UMAP_2))
  library(ggrepel)
  p5 <- p3 + geom_label_repel(aes(label=cell_type), fontface="bold",data = cell_type_med, point.padding=unit(1.5, "lines")) + theme(legend.position = "none")                  
  return(p5)
}
pdf("UMAP_Seurat.pdf",width=4,height=4)
UmapPlot(sub_OS)
dev.off()
#
UmapPlot<- function(obj) {
  library(ggplot2)
  library(dplyr)
  UMAP = obj@reductions$umap@cell.embeddings %>% as.data.frame() %>% cbind(cell_type = obj@meta.data$donor_status1) 
  allcolour=c("#DC143C","#0000FF","#20B2AA","#FFA500","#9370DB","#98FB98","#F08080","#1E90FF","#7CFC00","#FFFF00","#808000","#FF00FF","#FA8072","#7B68EE","#9400D3","#800080","#A0522D","#D2B48C","#D2691E","#87CEEB","#40E0D0","#5F9EA0","#FF1493","#0000CD","#008B8B","#FFE4B5","#8A2BE2","#228B22","#E9967A","#4682B4","#32CD32","#F0E68C","#FFFFE0","#EE82EE","#FF6347","#6A5ACD","#9932CC","#8B008B","#8B4513","#DEB887")
  my36colors <-c('#E5D2DD', '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87','#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#8C549C', '#585658', '#9FA3A8', '#E0D4CA', '#5F3D69', '#C5DEBA', '#58A4C3', '#E4C755', '#F7F398','#AA9A59', '#E63863', '#E39A35', '#C1E6F3', '#6778AE', '#91D0BE', '#B53E2B','#712820', '#DCC1DD', '#CCE0F5',  '#CCC9E6', '#625D9E', '#68A180', '#3A6963','#968175')
  p <- ggplot(UMAP,aes(x= UMAP_1, y = UMAP_2 ,color = cell_type)) +  geom_point(size=0.2, alpha =1 )  +  scale_color_manual(values = my36colors)
  p2 <- p  +theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),  panel.border = element_blank(),  axis.title = element_blank(),  axis.text = element_blank(), axis.ticks = element_blank(),panel.background = element_rect(fill = 'white'), plot.background=element_rect(fill="white"))
  p3 <- p2 + theme(legend.title = element_blank(),legend.key=element_rect(fill='white'), legend.text = element_text(size=20),legend.key.size=unit(1,'cm') ) + guides(color = guide_legend(override.aes = list(size=5))) 
#  p4 <- p3 + geom_segment(aes(x = min(TSNE$tSNE_1) , y = min(TSNE$tSNE_2) ,xend = min(TSNE$tSNE_1) +3, yend = min(TSNE$tSNE_2) ),colour = "black", size=1,arrow = arrow(length = unit(0.3,"cm")))+ geom_segment(aes(x = min(TSNE$tSNE_1), y = min(TSNE$tSNE_2),xend = min(TSNE$tSNE_1) , yend = min(TSNE$tSNE_2) + 3),colour = "black", size=1,arrow = arrow(length = unit(0.3,"cm"))) + annotate("text", x = min(TSNE$tSNE_1) +1.5, y = min(TSNE$tSNE_2) -1, label = "tSNE_1",color="black",size = 3, fontface="bold" ) + annotate("text", x = min(TSNE$tSNE_1) -1, y = min(TSNE$tSNE_2) + 1.5, label = "tSNE_2",color="black",size = 3, fontface="bold" ,angle=90) 
  cell_type_med <- UMAP %>%group_by(cell_type) %>%summarise(UMAP_1 = median(UMAP_1),UMAP_2 = median(UMAP_2))
  library(ggrepel)
  p5 <- p3 + geom_label_repel(aes(label=cell_type), fontface="bold",data = cell_type_med, point.padding=unit(1.5, "lines")) + theme(legend.position = "none")                  
  return(p5)
}
pdf("UMAP_condition.pdf",width=4,height=4)
UmapPlot(sub_OS)
dev.off()

######
###
##
bb <-NULL
sample <- unique(sub_OS@meta.data[,"donor_ID"])
for(i in 1:length(sample)){
aa <-table(sub_OS@meta.data[which(sub_OS@meta.data[,"donor_ID"]==sample[i]),"seurat_clusters"])
bb <-cbind(bb,aa)
}
colnames(bb) <-sample
cc <-NULL
for(i in 1:7){
aa <-bb[i,]/sum(bb[i,])
cc <-rbind(cc,aa)
}
rownames(cc) <-rownames(bb)
library(reshape2)
library(ggplot2)
library(RColorBrewer)
colorr <- brewer.pal(4,"Set2")
df <-melt(cc)
df$Var2 <-factor(df$Var2,levels=colnames(cc))
df$Var1 <-factor(df$Var1,levels=c("0","1","2","3","4","5","6"))
pdf(file="box_precent_sample.pdf",width=10,height=5)
ggplot(df,aes(x=Var1,y=value,fill=Var2))+scale_fill_manual(values=my36colors[1:17])+geom_bar(stat="identity")+theme(legend.title=element_blank()) 
dev.off()
##
bb <-NULL
sample <- unique(sub_OS@meta.data[,"donor_status2"])
for(i in 1:length(sample)){
aa <-table(sub_OS@meta.data[which(sub_OS@meta.data[,"donor_status2"]==sample[i]),"seurat_clusters"])
bb <-cbind(bb,aa)
}
colnames(bb) <-sample
cc <-NULL
for(i in 1:7){
aa <-bb[i,]/sum(bb[i,])
cc <-rbind(cc,aa)
}
rownames(cc) <-rownames(bb)
library(reshape2)
library(ggplot2)
library(RColorBrewer)
colorr <- brewer.pal(4,"Set2")
df <-melt(cc)
df$Var2 <-factor(df$Var2,levels=colnames(cc))
df$Var1 <-factor(df$Var1,levels=c("0","1","2","3","4","5","6"))
pdf(file="box_precent_condition.pdf",width=10,height=5)
ggplot(df,aes(x=Var1,y=value,fill=Var2))+scale_fill_manual(values=my36colors[1:3])+geom_bar(stat="identity")+theme(legend.title=element_blank()) 
dev.off()
###
library(reshape2)
library(ggplot2)
library(RColorBrewer)
dat <- sub_OS@meta.data[,c("donor_status2","seurat_clusters","donor_ID")]
dat1 <- dat[which(dat$donor_status2 =="Tumor"),]
bb <-NULL
sample <- unique(dat1[,"donor_ID"])
for(i in 1:length(sample)){
aa <-table(dat1[which(dat1[,"donor_ID"]==sample[i]),"seurat_clusters"])
bb <-cbind(bb,aa)}
colnames(bb) <-sample
cc <-t(t(bb)/as.vector(table(OS$donor_ID)[sample]))
dd1 <-melt(cc)
dd1$status <-"Tumor"
dat1 <- dat[which(dat$donor_status2 =="Chemo"),]
bb <-NULL
sample <- unique(dat1[,"donor_ID"])
for(i in 1:length(sample)){
aa <-table(dat1[which(dat1[,"donor_ID"]==sample[i]),"seurat_clusters"])
bb <-cbind(bb,aa)}
colnames(bb) <-sample
cc <-t(t(bb)/as.vector(table(OS$donor_ID)[sample]))
dd2 <-melt(cc)
dd2$status <-"Chemo"
#
dat1 <- dat[which(dat$donor_status2 =="Re_Me"),]
bb <-NULL
sample <- unique(dat1[,"donor_ID"])
for(i in 1:length(sample)){
aa <-table(dat1[which(dat1[,"donor_ID"]==sample[i]),"seurat_clusters"])
bb <-cbind(bb,aa)}
colnames(bb) <-sample
cc <-t(t(bb)/as.vector(table(OS$donor_ID)[sample]))
dd3 <-melt(cc)
dd3$status <-"Re_Me"
df <-rbind(dd1,dd2,dd3)
levels <-c("0","1","2","3","4","5","6")
df$status = factor(df$status, levels=c("Tumor","Chemo","Re_Me"))
dodge <- position_dodge(width = 0.8)
pdf(file="box_precent_stage.pdf",width=10,height=4)
ggplot(df, aes(x=factor(Var1,levels=levels), y=value,fill = status))  + stat_boxplot(geom="errorbar",width=0.3,position = dodge,size=2)+geom_boxplot(width=0.8,outlier.colour = NA)+scale_fill_manual(values=my36colors[1:15])+
ylim(0,0.3)+theme(legend.title=element_blank())+theme(axis.text.x = element_text(size = 12,angle = 45, hjust = 1, vjust = 1, face = "bold"))+labs(x="Cluster", y = "overlap")
dev.off()
#########
bb <-NULL
for(i in c(0:6)){
aa <- c(wilcox.test(df[which(df[,1]==i),][which(df[which(df[,1]==i),][,4]=="Tumor"),3],df[which(df[,1]==i),][which(df[which(df[,1]==i),][,4]=="Chemo"),3])$p.value,
wilcox.test(df[which(df[,1]==i),][which(df[which(df[,1]==i),][,4]=="Tumor"),3],df[which(df[,1]==i),][which(df[which(df[,1]==i),][,4]=="Re_Me"),3])$p.value,
wilcox.test(df[which(df[,1]==i),][which(df[which(df[,1]==i),][,4]=="Chemo"),3],df[which(df[,1]==i),][which(df[which(df[,1]==i),][,4]=="Re_Me"),3])$p.value)
bb <-rbind(bb,aa)
}
bb
aa 0.03496503 0.6095238 0.23030303
aa 0.01398601 0.7619048 0.07272727
aa 0.23426573 0.7619048 0.31515152
aa 0.03310613 0.7619048 0.21626897
aa 0.33722298 0.5093103 0.15481324
aa 0.03725991 0.2571429 0.43909761
aa 0.44040070 0.3074342 0.88853912
#### Go and KEGG function enrichment
marker <-read.table("markers_r0.1.txt")
library(ggplot2)
library(clusterProfiler)
library(org.Hs.eg.db)
for(i in 0:6){
ind <-which(marker[,"cluster"]==i)
gene <-marker[ind,"gene"]
up <-gene
gs = bitr(up, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
head(gs)
ego.bp = enrichGO(gene=gs$ENTREZID, OrgDb = org.Hs.eg.db,ont= "BP",pAdjustMethod = "BH",pvalueCutoff= 0.05,qvalueCutoff= 0.2,readable= TRUE)
head(ego.bp)
write.csv(ego.bp, file =paste(i,"_go.csv"))
pdf(paste(i,"_go.pdf"))
print(dotplot(ego.bp, showCategory=10,title="up gene GoTerm"))
dev.off()
kk <- enrichKEGG(gene= gs$ENTREZID, organism= 'hsa',pvalueCutoff = 0.05)
write.csv(kk,file =paste(i,"_kegg.csv"))
pdf(paste(i,"_kegg.pdf"))
print(dotplot(kk, showCategory=10,title="KEGG_Up_biological"))
dev.off()
}
###
marker <-read.table("markers_r0.1.txt")
marker1 <- marker[-grep("^RPS|RPL|MT-",marker$gene),]
gene1 <-NULL
for(i in c("0","1","2","3","4","5","6")){
index <- which(marker1$cluster==i)
gene2 <-marker1[index,"gene"][1:10]
gene1 <-c(gene1,gene2)
}
gene <-unique(gene1)
##gene <-c(gene[1:24],"MKI67",gene[26:70])
sub_OS_data1 <-sub_OS[["RNA"]]@data[gene,]
sub_OS_data <-t(scale(t(as.matrix(sub_OS_data1))))
type <-unique(sub_OS@meta.data[,"seurat_clusters"])
bb <-NULL
for(i in 1: length(type)){
ind <-which(sub_OS@meta.data[,"seurat_clusters"]==type[i])
aa <-rowMeans(as.matrix(sub_OS_data[,ind]))
bb <-cbind(bb,aa)
}
colnames(bb) <-type
library(reshape2)
library(ggplot2)
bb[bb>3.8] <-3.8
data_melt<-melt(bb)
data_melt$Var2 <-factor(data_melt$Var2,levels=c("0","1","2","3","4","5","6"))
p <- ggplot(data = data_melt) + geom_tile(aes(x=Var2,y=Var1, fill = value)) +theme_classic() +  theme(axis.ticks = element_blank(), axis.line = element_blank()) + xlab('row name') + ylab('column name')
p1 <- p+scale_fill_gradient2('z-score',low = 'blue', high = 'red', mid = 'white') 
pdf("heatmap_gene_celltype.pdf",width=4,height=8)			
p1
dev.off()
##############
library(presto)
library(msigdbr)
library(fgsea)
library(dplyr)
library(ggplot2)
library(tibble)
# select Hallmarks
m_df<- msigdbr(species = "Homo sapiens", category = "H") 
fgsea_sets = m_df %>% split(x = .$gene_symbol, f = .$gs_name)
summary(fgsea_sets)
pbmc.genes <- wilcoxauc(sub_OS, 'seurat_clusters')
dplyr::count(pbmc.genes, group)
for(i in 0:6){
cluster.genes<- pbmc.genes %>% dplyr::filter(group == i) %>% arrange(desc(auc)) %>% dplyr::select(feature, auc) 
ranks<- deframe(cluster.genes) 
head(ranks)
fgseaRes<- fgsea(fgsea_sets, stats = ranks, nperm = 1000)
setwd("./Analysis/Osteoblasts/GSEA")
fgseaResTidy <- fgseaRes %>% as_tibble() %>% arrange(desc(NES)) 
fgseaResTidy %>% dplyr::select(-leadingEdge, -ES, -nMoreExtreme) %>% arrange(padj) %>% head()
saveRDS(fgseaResTidy,paste(i,'_GSEA-fgsea.rds'))
p <-ggplot(fgseaRes %>% as_tibble() %>% arrange(desc(NES)) %>% filter(pval < 0.05), aes(reorder(pathway, NES), NES)) +geom_col(aes(fill= NES)) +coord_flip() +labs(x="KEGG", y="Normalized Enrichment Score",title="KEGG gene sets NES from GSEA") ##输出差异排秩前20的条目
pdf(paste(i,'_GSEA-fgsea.pdf'))
print(p)
dev.off()
}
#
aa <-readRDS("0 _GSEA-fgsea.rds")
aa1 <-readRDS("1 _GSEA-fgsea.rds")
aa2 <-readRDS("2 _GSEA-fgsea.rds")
aa3 <-readRDS("3 _GSEA-fgsea.rds")
aa4 <-readRDS("4 _GSEA-fgsea.rds")
aa5 <-readRDS("5 _GSEA-fgsea.rds")
aa6 <-readRDS("6 _GSEA-fgsea.rds")
bb <-aa %>% as_tibble() %>% arrange(desc(NES)) %>% filter(pval < 0.05)
cc <-cbind(bb$pathway,bb$NES,"C0")
bb1 <-aa1 %>% as_tibble() %>% arrange(desc(NES)) %>% filter(pval < 0.05)
cc1 <-cbind(bb1$pathway,bb1$NES,"C1")
bb2 <-aa2 %>% as_tibble() %>% arrange(desc(NES)) %>% filter(pval < 0.05)
cc2 <-cbind(bb2$pathway,bb2$NES,"C2")
bb3 <-aa3 %>% as_tibble() %>% arrange(desc(NES)) %>% filter(pval < 0.05)
cc3 <-cbind(bb3$pathway,bb3$NES,"C3")
bb4 <-aa4 %>% as_tibble() %>% arrange(desc(NES)) %>% filter(pval < 0.05)
cc4 <-cbind(bb4$pathway,bb4$NES,"C4")
bb5 <-aa5 %>% as_tibble() %>% arrange(desc(NES)) %>% filter(pval < 0.05)
cc5 <-cbind(bb5$pathway,bb5$NES,"C5")
bb6 <-aa6 %>% as_tibble() %>% arrange(desc(NES)) %>% filter(pval < 0.05)
cc6 <-cbind(bb6$pathway,bb6$NES,"C6")
dd <-data.frame(rbind(cc,cc1,cc2,cc3,cc4,cc5,cc6))
library(pheatmap)
mat <-matrix(0,length(unique(dd$X1)),length(unique(dd$X3)))
rownames(mat) <-unique(dd$X1)
colnames(mat) <-unique(dd$X3)
for(i in 1:nrow(dd)){
mat[dd[i,1],dd[i,3]] <-as.numeric(dd[i,2])
}
bk <- c(seq(-12,-0.01,by=0.01),seq(0.01,11.2,by=0.01))
rownames(mat) <-gsub("kegg ","",gsub("_"," ",tolower(rownames(mat))))
pdf("pheatmap_Func_Hallmark.pdf",width=15,height=15)
pheatmap(mat,cluster_rows = T, cluster_cols=F,show_rownames = T,main = "Heatmap",color = c(colorRampPalette(colors = c('#0A7EB5','#F0F8FF'))(120),colorRampPalette(colors = c('#F0F8FF','#EE6A50'))(112)))
dev.off()
###############
###
DefaultAssay(sub_OS) <- "RNA"
sub_OS@active.ident <-factor(as.matrix(sub_OS@meta.data)[,"donor_status2"])
OS.markers <- FindAllMarkers(object = sub_OS, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.table(OS.markers,"markers_donor_status.txt",row.names=TRUE,col.names=TRUE)
###
###
marker <-OS.markers
marker1 <- marker[-grep("^RPS|RPL|MT-",marker$gene),]
gene1 <-NULL
for(i in c("Tumor","Chemo","Re_Me")){
index <- which(marker1$cluster==i)
gene2 <-marker1[index,"gene"][1:10]
gene1 <-c(gene1,gene2)
}
gene <-unique(gene1)
sub_OS_data1 <-sub_OS[["RNA"]]@data[gene,]
sub_OS_data <-t(scale(t(as.matrix(sub_OS_data1))))
type <-unique(sub_OS@meta.data[,"donor_status2"])
bb <-NULL
for(i in 1: length(type)){
ind <-which(sub_OS@meta.data[,"donor_status2"]==type[i])
aa <-rowMeans(as.matrix(sub_OS_data[,ind]))
bb <-cbind(bb,aa)
}
colnames(bb) <-type
library(reshape2)
library(ggplot2)
bb[bb>3.8] <-3.8
data_melt<-melt(bb)
data_melt$Var2 <-factor(data_melt$Var2,levels=c("Tumor","Chemo","Re_Me"))
p <- ggplot(data = data_melt) + geom_tile(aes(x=Var2,y=Var1, fill = value)) +theme_classic() +  theme(axis.ticks = element_blank(), axis.line = element_blank()) + xlab('row name') + ylab('column name')
p1 <- p+scale_fill_gradient2('z-score',low = 'blue', high = 'red', mid = 'white') 
pdf("heatmap_gene_stage.pdf",width=3,height=8)			
p1
dev.off()
###############
library(ggplot2)
library(clusterProfiler)
library(org.Hs.eg.db)
for(i in c("Tumor","Chemo","Re_Me")){
ind <-which(marker[,"cluster"]==i)
gene <-marker[ind,"gene"]
up <-gene
gs = bitr(up, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
head(gs)
ego.bp = enrichGO(gene=gs$ENTREZID, OrgDb = org.Hs.eg.db,ont= "BP",pAdjustMethod = "BH",pvalueCutoff= 0.05,qvalueCutoff= 0.2,readable= TRUE)
head(ego.bp)
write.csv(ego.bp, file =paste(i,"_go.csv"))
pdf(paste(i,"_go.pdf"))
print(dotplot(ego.bp, showCategory=10,title="up gene GoTerm"))
dev.off()
pdf(paste(i,"_go_box.pdf"))
print(barplot(ego.bp, showCategory=10,title="up gene GoTerm"))
dev.off()
kk <- enrichKEGG(gene= gs$ENTREZID, organism= 'hsa',pvalueCutoff = 0.05)
write.csv(kk,file =paste(i,"_kegg.csv"))
pdf(paste(i,"_kegg.pdf"))
print(dotplot(kk, showCategory=10,title="KEGG_Up_biological"))
dev.off()
}
##

#######################
DefaultAssay(sub_OS) <- "RNA"
sub_OS@active.ident <-factor(as.matrix(sub_OS@meta.data)[,"donor_status2"])
cluster.markers <- FindMarkers(sub_OS, ident.1 = c("Chemo"), ident.2 =c("Tumor"), min.pct = 0.1,logfc.threshold=0)
filt_data <-cluster.markers
data<-c()
for(i in 1:nrow(filt_data)){
	if(filt_data[i,1]<0.05){
		if(filt_data[i,2]>0.25){
			data<-rbind(data,c(filt_data[i,2],filt_data[i,1],"up"))
		}
		else if(filt_data[i,2]<(-0.25)){
			data<-rbind(data,c(filt_data[i,2],filt_data[i,1],"down"))
		}
		else{
		    data<-rbind(data,c(filt_data[i,2],filt_data[i,1],"no"))
	    }
	}	
	else{
		data<-rbind(data,c(filt_data[i,2],filt_data[i,1],"no"))
	}
}
table(data[,3])
DEG_num=as.numeric(table(data[,3])[1]+table(data[,3])[3])
sub_title=paste("DEGs:",DEG_num,sep="")
colnames(data)<-c("logFC","pvalue",sub_title)
rownames(data)<-rownames(filt_data)
write.table(data,file="Cluster_volcano_plot_format.txt",quote=F,sep="\t")
data <-read.table("Cluster_volcano_plot_format.txt",header=1,sep='\t')
up_num=length(intersect(which(data[,2]<0.05),which(data[,1]>0.25)))
down_num=length(intersect(which(data[,2]<0.05),which(data[,1]<(-0.25))))
legend1=paste("up regulated:",up_num,sep=" ")
legend2=paste("down regulated:",down_num,sep=" ")
index <-which(data[,3]=="no")
index1 <-which(data[,3]=="up")
index2 <-which(data[,3]=="down")
ind <-c(index1,index2,index)
data <-data[ind,]
library("ggplot2")
r03=ggplot(data,aes(x=logFC,y=-log10(pvalue))) +geom_point()+ theme_bw()+theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())
r03=r03 + geom_point(aes(color =DEGs.1970))
volcano = r03 +scale_color_manual(values = c("#009966","#0066CC","#FF0033"),breaks=c("up","down"),labels=c(legend1,legend2))+theme(legend.position = c(.25, .75))
volcano=volcano+geom_hline(yintercept=-log10(0.05),linetype=2,size=1)+geom_vline(xintercept=c((-0.25),0.25),linetype=2,size=1)
library(psych)
library(Seurat)
library(dplyr)
library(ggrepel)
data$sign <- ifelse(data$logFC >= 0.8 | data$logFC <= -0.8,rownames(data),NA)
#data$sign <- ifelse( -log10(data$pvalue) >= 10,rownames(data),NA)
pdf(file="volcano_Chemo_Tumor.pdf",width=4,height=5)	
volcano +geom_text_repel(label=data$sign,colour="black",size=2,box.padding = unit(0.3, "lines"), point.padding = unit(0.4, "lines"), show.legend = F) 
dev.off()
##
library(clusterProfiler)
library(org.Hs.eg.db)
data <-read.table("Cluster_volcano_plot_format.txt",header=1,sep='\t')
up <-rownames(data[intersect(which(data[,2]<0.05),which(data[,1]>0.25)),])
down <-rownames(data[intersect(which(data[,2]<0.05),which(data[,1]<(-0.25))),])
gs <- up
gs = bitr(gs, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
head(gs)
ego.bp = enrichGO(gene =  gs$ENTREZID,OrgDb = org.Hs.eg.db,ont  = "BP",pAdjustMethod = "BH",pvalueCutoff  = 0.05,qvalueCutoff  = 0.2,readable = TRUE)					
head(ego.bp)
# 
nam = paste("Chemo_Tumor_UP","Go.csv",sep ="")
write.csv(ego.bp, file =nam)
nam = paste("Chemo_Tumor_UP","Go.pdf",sep ="")
pdf(nam,width=8,height=6)
print(barplot(ego.bp, showCategory=9,title="GO_biological",drop=T))
dev.off()
#
gs <- down
gs = bitr(gs, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
head(gs)
ego.bp = enrichGO(gene =  gs$ENTREZID,OrgDb = org.Hs.eg.db,ont  = "BP",pAdjustMethod = "BH",pvalueCutoff  = 0.05,qvalueCutoff  = 0.2,readable = TRUE)					
# 
nam = paste("Chemo_Tumor_DOWN","Go.csv",sep ="")
write.csv(ego.bp, file =nam)
nam = paste("Chemo_Tumor_DOWN","Go.pdf",sep ="")
pdf(nam,width=12,height=6)
print(barplot(ego.bp, showCategory=10,title="GO_biological",drop=T))
dev.off()

bp1 <-read.csv("Chemo_Tumor_UPGo.csv")
bp2 <-read.csv("Chemo_Tumor_DOWNGo.csv")
bp1$celltype <-"UP"
bp2$celltype <-"DOWN"
bp <-rbind(bp1[c(2,3,4,6,7),],bp2[c(1,3,12,14,15),])
levels <-rev(bp$Description)
pdf("Chemo_Tumor_UP_DOWN_GO_function.pdf",width=10,height=11)
ggplot(bp,aes(x=factor(Description,levels=levels),y=Count,fill=pvalue)) + 
  geom_bar(stat='identity',color='black',width = 0.65) +
  coord_flip() +    
  scale_fill_gradient(low='#FFFFF0', high='darkgoldenrod1')
dev.off()

######
cluster.markers <- FindMarkers(sub_OS, ident.1 = c("Re_Me"), ident.2 =c("Chemo"), min.pct = 0.1,logfc.threshold=0)
filt_data <-cluster.markers
data<-c()
for(i in 1:nrow(filt_data)){
	if(filt_data[i,1]<0.05){
		if(filt_data[i,2]>0.25){
			data<-rbind(data,c(filt_data[i,2],filt_data[i,1],"up"))
		}
		else if(filt_data[i,2]<(-0.25)){
			data<-rbind(data,c(filt_data[i,2],filt_data[i,1],"down"))
		}
		else{
		    data<-rbind(data,c(filt_data[i,2],filt_data[i,1],"no"))
	    }
	}	
	else{
		data<-rbind(data,c(filt_data[i,2],filt_data[i,1],"no"))
	}
}
table(data[,3])
DEG_num=as.numeric(table(data[,3])[1]+table(data[,3])[3])
sub_title=paste("DEGs:",DEG_num,sep="")
colnames(data)<-c("logFC","pvalue",sub_title)
rownames(data)<-rownames(filt_data)
write.table(data,file="Cluster_volcano_plot_format1.txt",quote=F,sep="\t")
data <-read.table("Cluster_volcano_plot_format1.txt",header=1,sep='\t')
up_num=length(intersect(which(data[,2]<0.05),which(data[,1]>0.25)))
down_num=length(intersect(which(data[,2]<0.05),which(data[,1]<(-0.25))))
legend1=paste("up regulated:",up_num,sep=" ")
legend2=paste("down regulated:",down_num,sep=" ")
index <-which(data[,3]=="no")
index1 <-which(data[,3]=="up")
index2 <-which(data[,3]=="down")
ind <-c(index1,index2,index)
data <-data[ind,]
library("ggplot2")
r03=ggplot(data,aes(x=logFC,y=-log10(pvalue))) +geom_point()+ theme_bw()+theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())
r03=r03 + geom_point(aes(color =DEGs.2025))
volcano = r03 +scale_color_manual(values = c("#009966","#0066CC","#FF0033"),breaks=c("up","down"),labels=c(legend1,legend2))+theme(legend.position = c(.25, .75))
volcano=volcano+geom_hline(yintercept=-log10(0.05),linetype=2,size=1)+geom_vline(xintercept=c((-0.25),0.25),linetype=2,size=1)
library(psych)
library(Seurat)
library(dplyr)
library(ggrepel)
data$sign <- ifelse(data$logFC >= 0.8 | data$logFC <= -0.8,rownames(data),NA)
#data$sign <- ifelse( -log10(data$pvalue) >= 10,rownames(data),NA)
pdf(file="volcano_Re_Chemo.pdf",width=4,height=5)	
volcano +geom_text_repel(label=data$sign,colour="black",size=2,box.padding = unit(0.3, "lines"), point.padding = unit(0.4, "lines"), show.legend = F) 
dev.off()
##
library(clusterProfiler)
library(org.Hs.eg.db)
data <-read.table("Cluster_volcano_plot_format1.txt",header=1,sep='\t')
up <-rownames(data[intersect(which(data[,2]<0.05),which(data[,1]>0.25)),])
down <-rownames(data[intersect(which(data[,2]<0.05),which(data[,1]<(-0.25))),])
gs <- up
gs = bitr(gs, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
head(gs)
ego.bp = enrichGO(gene =  gs$ENTREZID,OrgDb = org.Hs.eg.db,ont  = "BP",pAdjustMethod = "BH",pvalueCutoff  = 0.05,qvalueCutoff  = 0.2,readable = TRUE)					
head(ego.bp)
# 
nam = paste("Re_Chemo_UP","Go.csv",sep ="")
write.csv(ego.bp, file =nam)
nam = paste("Re_Chemo_UP","Go.pdf",sep ="")
pdf(nam,width=8,height=6)
print(barplot(ego.bp, showCategory=9,title="GO_biological",drop=T))
dev.off()
#
gs <- down
gs = bitr(gs, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
head(gs)
ego.bp = enrichGO(gene =  gs$ENTREZID,OrgDb = org.Hs.eg.db,ont  = "BP",pAdjustMethod = "BH",pvalueCutoff  = 0.05,qvalueCutoff  = 0.2,readable = TRUE)					
# 
nam = paste("Re_Chemo_DOWN","Go.csv",sep ="")
write.csv(ego.bp, file =nam)
nam = paste("Re_Chemo_DOWN","Go.pdf",sep ="")
pdf(nam,width=12,height=6)
print(barplot(ego.bp, showCategory=10,title="GO_biological",drop=T))
dev.off()

bp1 <-read.csv("Re_Chemo_UPGo.csv")
bp2 <-read.csv("Re_Chemo_DOWNGo.csv")
bp1$celltype <-"UP"
bp2$celltype <-"DOWN"
bp <-rbind(bp1[c(2,3,17,120,128),],bp2[c(1,2,8,14,17),])
levels <-rev(bp$Description)
pdf("Re_Chemo_UP_DOWN_GO_function.pdf",width=10,height=11)
ggplot(bp,aes(x=factor(Description,levels=levels),y=Count,fill=pvalue)) + 
  geom_bar(stat='identity',color='black',width = 0.65) +
  coord_flip() +    
  scale_fill_gradient(low='#FFFFF0', high='darkgoldenrod1')
dev.off()
####
###### copykat
library(Seurat)
library(copykat)
library(tidyverse)
rm(list = ls())
setwd("./Analysis/Osteoblasts/copykat")
OS@active.ident <-factor(as.matrix(OS@meta.data)[,"donor_ID"])
sub1 <-subset(OS,idents=c("OS3_"))
setwd("./Analysis/Osteoblasts/copykat/OS3")
scRNA <-sub1
counts <- as.matrix(scRNA@assays$RNA@counts)
cnv <- copykat(rawmat=counts, ngene.chr=10, sam.name="Epi", n.cores=8)

saveRDS(cnv, "cnv.rds")
cnv <-readRDS("./Osteoblastic/copykat/OS3/cnv.rds")
pred.test <- data.frame(cnv$prediction)
pred.test <- pred.test[-which(pred.test$copykat.pred == "not.defined"), ]
scRNA1 <- subset(scRNA, cells =  pred.test[,1] )
scRNA1@meta.data$copykat.pred <- pred.test$copykat.pred
scRNA1@meta.data$copykat.tumor.pred <- pred.test$copykat.pred
p1 <- DimPlot(scRNA1, label = T,reduction = "umap",group.by ="celltype")
p2 <- DimPlot(scRNA1, group.by = "copykat.pred",reduction = "umap")
pc <-p1 + p2 
ggsave("pred_mallignant.pdf", pc, width = 12, height = 5)

##
setwd("./Analysis/Osteoblasts/copykat/BC5")
OS@active.ident <-factor(as.matrix(OS@meta.data)[,"donor_ID"])
sub1 <-subset(OS,idents=c("BC5_"))
scRNA <-sub1
cnv <-readRDS("./Osteoblastic/copykat/BC5/cnv.rds")
pred.test <- data.frame(cnv$prediction)
pred.test <- pred.test[-which(pred.test$copykat.pred == "not.defined"), ]
scRNA1 <- subset(scRNA, cells =  pred.test[,1] )
scRNA1@meta.data$copykat.pred <- pred.test$copykat.pred
scRNA1@meta.data$copykat.tumor.pred <- pred.test$copykat.pred
p1 <- DimPlot(scRNA1, label = T,reduction = "umap",group.by ="celltype")
p2 <- DimPlot(scRNA1, group.by = "copykat.pred",reduction = "umap")
p1 + p2 
ggsave("pred_mallignant.pdf", pc, width = 12, height = 5)
##
setwd("./Analysis/Osteoblasts/copykat/BC22")
sub1 <-subset(OS,idents=c("BC22"))
scRNA <-sub1
cnv <-readRDS("./Osteoblastic/copykat/BC22/cnv.rds")
pred.test <- data.frame(cnv$prediction)
pred.test <- pred.test[-which(pred.test$copykat.pred == "not.defined"), ]
scRNA1 <- subset(scRNA, cells =  pred.test[,1] )
scRNA1@meta.data$copykat.pred <- pred.test$copykat.pred
scRNA1@meta.data$copykat.tumor.pred <- pred.test$copykat.pred
p1 <- DimPlot(scRNA1, label = T,reduction = "umap",group.by ="celltype")
p2 <- DimPlot(scRNA1, group.by = "copykat.pred",reduction = "umap")
pc <-p1 + p2 
ggsave("pred_mallignant.pdf", pc, width = 12, height = 5)
#
setwd("./Analysis/Osteoblasts/copykat/BC11")
sub1 <-subset(OS,idents=c("BC11"))
scRNA <-sub1
cnv <-readRDS("./Osteoblastic/copykat/BC11/cnv.rds")
pred.test <- data.frame(cnv$prediction)
pred.test <- pred.test[-which(pred.test$copykat.pred == "not.defined"), ]
scRNA1 <- subset(scRNA, cells =  pred.test[,1] )
scRNA1@meta.data$copykat.pred <- pred.test$copykat.pred
scRNA1@meta.data$copykat.tumor.pred <- pred.test$copykat.pred
p1 <- DimPlot(scRNA1, label = T,reduction = "umap",group.by ="celltype")
p2 <- DimPlot(scRNA1, group.by = "copykat.pred",reduction = "umap")
pc <-p1 + p2 
ggsave("pred_mallignant.pdf", pc, width = 12, height = 5)
#













$######
###########
######
######
########
############
##############
###############Osteoclasts
setwd("./Analysis/Osteoclasts")
OS@active.ident <-factor(as.matrix(OS@meta.data)[,"seurat_clusters"])
sub_OS <-subset(OS,idents= c("5"))
DefaultAssay(sub_OS) <- "integrated"
# Run the standard workflow for visualization and clustering
#sub_OS <- ScaleData(sub_OS, verbose = FALSE)
sub_OS <- RunPCA(sub_OS, npcs = 30, verbose = FALSE)
sub_OS <- RunUMAP(sub_OS, reduction = "pca", dims = 1:30)
sub_OS <- FindNeighbors(sub_OS, reduction = "pca", dims = 1:30)
sub_OS <- FindClusters(sub_OS, resolution = 0.2)
p0 <- DimPlot(object = sub_OS, reduction = "umap",label=TRUE)
p1<- DimPlot(object = sub_OS, reduction = "umap",group.by ="donor_ID",label=TRUE)
p2 <- DimPlot(object = sub_OS, reduction = "umap",group.by ="donor_status",label=TRUE)
pdf("Umap_sub_OS.pdf",width=22,height=20)
CombinePlots(plots  = list(p0, p1, p2))
dev.off()
saveRDS(sub_OS, file="data_rPCA_sub_OS.rds")
#
sub_OS <-readRDS("data_rPCA_sub_OS.rds")
UmapPlot<- function(obj) {
  library(ggplot2)
  library(dplyr)
  UMAP = obj@reductions$umap@cell.embeddings %>% as.data.frame() %>% cbind(cell_type = obj@meta.data$seurat_clusters) 
  allcolour=c("#DC143C","#0000FF","#20B2AA","#FFA500","#9370DB","#98FB98","#F08080","#1E90FF","#7CFC00","#FFFF00","#808000","#FF00FF","#FA8072","#7B68EE","#9400D3","#800080","#A0522D","#D2B48C","#D2691E","#87CEEB","#40E0D0","#5F9EA0","#FF1493","#0000CD","#008B8B","#FFE4B5","#8A2BE2","#228B22","#E9967A","#4682B4","#32CD32","#F0E68C","#FFFFE0","#EE82EE","#FF6347","#6A5ACD","#9932CC","#8B008B","#8B4513","#DEB887")
  my36colors <-c('#E5D2DD', '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87','#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#8C549C', '#585658', '#9FA3A8', '#E0D4CA', '#5F3D69', '#C5DEBA', '#58A4C3', '#E4C755', '#F7F398','#AA9A59', '#E63863', '#E39A35', '#C1E6F3', '#6778AE', '#91D0BE', '#B53E2B','#712820', '#DCC1DD', '#CCE0F5',  '#CCC9E6', '#625D9E', '#68A180', '#3A6963','#968175')
  p <- ggplot(UMAP,aes(x= UMAP_1, y = UMAP_2 ,color = cell_type)) +  geom_point(size=0.2, alpha =1 )  +  scale_color_manual(values = my36colors)
  p2 <- p  +theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),  panel.border = element_blank(),  axis.title = element_blank(),  axis.text = element_blank(), axis.ticks = element_blank(),panel.background = element_rect(fill = 'white'), plot.background=element_rect(fill="white"))
  p3 <- p2 + theme(legend.title = element_blank(),legend.key=element_rect(fill='white'), legend.text = element_text(size=20),legend.key.size=unit(1,'cm') ) + guides(color = guide_legend(override.aes = list(size=5))) 
#  p4 <- p3 + geom_segment(aes(x = min(TSNE$tSNE_1) , y = min(TSNE$tSNE_2) ,xend = min(TSNE$tSNE_1) +3, yend = min(TSNE$tSNE_2) ),colour = "black", size=1,arrow = arrow(length = unit(0.3,"cm")))+ geom_segment(aes(x = min(TSNE$tSNE_1), y = min(TSNE$tSNE_2),xend = min(TSNE$tSNE_1) , yend = min(TSNE$tSNE_2) + 3),colour = "black", size=1,arrow = arrow(length = unit(0.3,"cm"))) + annotate("text", x = min(TSNE$tSNE_1) +1.5, y = min(TSNE$tSNE_2) -1, label = "tSNE_1",color="black",size = 3, fontface="bold" ) + annotate("text", x = min(TSNE$tSNE_1) -1, y = min(TSNE$tSNE_2) + 1.5, label = "tSNE_2",color="black",size = 3, fontface="bold" ,angle=90) 
  cell_type_med <- UMAP %>%group_by(cell_type) %>%summarise(UMAP_1 = median(UMAP_1),UMAP_2 = median(UMAP_2))
  library(ggrepel)
  p5 <- p3 + geom_label_repel(aes(label=cell_type), fontface="bold",data = cell_type_med, point.padding=unit(1.5, "lines")) + theme(legend.position = "none")                  
  return(p5)
}
pdf("UMAP_Seurat.pdf",width=4,height=4)
UmapPlot(sub_OS)
dev.off()

UmapPlot<- function(obj) {
  library(ggplot2)
  library(dplyr)
  UMAP = obj@reductions$umap@cell.embeddings %>% as.data.frame() %>% cbind(cell_type = obj@meta.data$donor_ID) 
  allcolour=c("#DC143C","#0000FF","#20B2AA","#FFA500","#9370DB","#98FB98","#F08080","#1E90FF","#7CFC00","#FFFF00","#808000","#FF00FF","#FA8072","#7B68EE","#9400D3","#800080","#A0522D","#D2B48C","#D2691E","#87CEEB","#40E0D0","#5F9EA0","#FF1493","#0000CD","#008B8B","#FFE4B5","#8A2BE2","#228B22","#E9967A","#4682B4","#32CD32","#F0E68C","#FFFFE0","#EE82EE","#FF6347","#6A5ACD","#9932CC","#8B008B","#8B4513","#DEB887")
  my36colors <-c('#E5D2DD', '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87','#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#8C549C', '#585658', '#9FA3A8', '#E0D4CA', '#5F3D69', '#C5DEBA', '#58A4C3', '#E4C755', '#F7F398','#AA9A59', '#E63863', '#E39A35', '#C1E6F3', '#6778AE', '#91D0BE', '#B53E2B','#712820', '#DCC1DD', '#CCE0F5',  '#CCC9E6', '#625D9E', '#68A180', '#3A6963','#968175') 
  p <- ggplot(UMAP,aes(x= UMAP_1, y = UMAP_2 ,color = cell_type)) +  geom_point(size=0.2, alpha =1 )  +  scale_color_manual(values = my36colors)
  p2 <- p  +theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),  panel.border = element_blank(),  axis.title = element_blank(),  axis.text = element_blank(), axis.ticks = element_blank(),panel.background = element_rect(fill = 'white'), plot.background=element_rect(fill="white"))
  p3 <- p2 + theme(legend.title = element_blank(),legend.key=element_rect(fill='white'), legend.text = element_text(size=20),legend.key.size=unit(1,'cm') ) + guides(color = guide_legend(override.aes = list(size=5))) 
#  p4 <- p3 + geom_segment(aes(x = min(TSNE$tSNE_1) , y = min(TSNE$tSNE_2) ,xend = min(TSNE$tSNE_1) +3, yend = min(TSNE$tSNE_2) ),colour = "black", size=1,arrow = arrow(length = unit(0.3,"cm")))+ geom_segment(aes(x = min(TSNE$tSNE_1), y = min(TSNE$tSNE_2),xend = min(TSNE$tSNE_1) , yend = min(TSNE$tSNE_2) + 3),colour = "black", size=1,arrow = arrow(length = unit(0.3,"cm"))) + annotate("text", x = min(TSNE$tSNE_1) +1.5, y = min(TSNE$tSNE_2) -1, label = "tSNE_1",color="black",size = 3, fontface="bold" ) + annotate("text", x = min(TSNE$tSNE_1) -1, y = min(TSNE$tSNE_2) + 1.5, label = "tSNE_2",color="black",size = 3, fontface="bold" ,angle=90) 
  cell_type_med <- UMAP %>%group_by(cell_type) %>%summarise(UMAP_1 = median(UMAP_1),UMAP_2 = median(UMAP_2))
  library(ggrepel)
  p5 <- p3 + geom_label_repel(aes(label=cell_type), fontface="bold",data = cell_type_med, point.padding=unit(1.5, "lines")) + theme(legend.position = "none")                  
  return(p5)
}
pdf("UMAP_ID.pdf",width=5,height=5)
UmapPlot(sub_OS)
dev.off()

UmapPlot<- function(obj) {
  library(ggplot2)
  library(dplyr)
  UMAP = obj@reductions$umap@cell.embeddings %>% as.data.frame() %>% cbind(cell_type = obj@meta.data$donor_status1) 
  allcolour=c("#DC143C","#0000FF","#20B2AA","#FFA500","#9370DB","#98FB98","#F08080","#1E90FF","#7CFC00","#FFFF00","#808000","#FF00FF","#FA8072","#7B68EE","#9400D3","#800080","#A0522D","#D2B48C","#D2691E","#87CEEB","#40E0D0","#5F9EA0","#FF1493","#0000CD","#008B8B","#FFE4B5","#8A2BE2","#228B22","#E9967A","#4682B4","#32CD32","#F0E68C","#FFFFE0","#EE82EE","#FF6347","#6A5ACD","#9932CC","#8B008B","#8B4513","#DEB887")
  my36colors <-c('#E5D2DD', '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87','#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#8C549C', '#585658', '#9FA3A8', '#E0D4CA', '#5F3D69', '#C5DEBA', '#58A4C3', '#E4C755', '#F7F398','#AA9A59', '#E63863', '#E39A35', '#C1E6F3', '#6778AE', '#91D0BE', '#B53E2B','#712820', '#DCC1DD', '#CCE0F5',  '#CCC9E6', '#625D9E', '#68A180', '#3A6963','#968175') 
  p <- ggplot(UMAP,aes(x= UMAP_1, y = UMAP_2 ,color = cell_type)) +  geom_point(size=0.2, alpha =1 )  +  scale_color_manual(values = my36colors)
  p2 <- p  +theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),  panel.border = element_blank(),  axis.title = element_blank(),  axis.text = element_blank(), axis.ticks = element_blank(),panel.background = element_rect(fill = 'white'), plot.background=element_rect(fill="white"))
  p3 <- p2 + theme(legend.title = element_blank(),legend.key=element_rect(fill='white'), legend.text = element_text(size=20),legend.key.size=unit(1,'cm') ) + guides(color = guide_legend(override.aes = list(size=5))) 
#  p4 <- p3 + geom_segment(aes(x = min(TSNE$tSNE_1) , y = min(TSNE$tSNE_2) ,xend = min(TSNE$tSNE_1) +3, yend = min(TSNE$tSNE_2) ),colour = "black", size=1,arrow = arrow(length = unit(0.3,"cm")))+ geom_segment(aes(x = min(TSNE$tSNE_1), y = min(TSNE$tSNE_2),xend = min(TSNE$tSNE_1) , yend = min(TSNE$tSNE_2) + 3),colour = "black", size=1,arrow = arrow(length = unit(0.3,"cm"))) + annotate("text", x = min(TSNE$tSNE_1) +1.5, y = min(TSNE$tSNE_2) -1, label = "tSNE_1",color="black",size = 3, fontface="bold" ) + annotate("text", x = min(TSNE$tSNE_1) -1, y = min(TSNE$tSNE_2) + 1.5, label = "tSNE_2",color="black",size = 3, fontface="bold" ,angle=90) 
  cell_type_med <- UMAP %>%group_by(cell_type) %>%summarise(UMAP_1 = median(UMAP_1),UMAP_2 = median(UMAP_2))
  library(ggrepel)
  p5 <- p3 + geom_label_repel(aes(label=cell_type), fontface="bold",data = cell_type_med, point.padding=unit(1.5, "lines")) + theme(legend.position = "none")                  
  return(p5)
}
pdf("UMAP_stage.pdf",width=5,height=5)
UmapPlot(sub_OS)
dev.off()

######
DefaultAssay(sub_OS) <- "RNA"
OS.markers <- FindAllMarkers(object = sub_OS, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.table(OS.markers,"markers_r0.2.txt",row.names=TRUE,col.names=TRUE)
DefaultAssay(sub_OS) <- "RNA"
sub_OS@active.ident <-factor(as.matrix(sub_OS@meta.data)[,"donor_status2"])
OS.markers <- FindAllMarkers(object = sub_OS, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.table(OS.markers,"markers_r0.2_donor_status.txt",row.names=TRUE,col.names=TRUE)

#########
library(ggplot2)
library(patchwork)
marker <-read.table("./Analysis/Osteoclasts/markers_r0.2.txt")
marker1 <- marker[-grep("^RPS|RPL|MT-",marker$gene),]
gene1 <-NULL
for(i in c("0","1","2","3","4")){
index <- which(marker1$cluster==i)
gene2 <-marker1[index,"gene"][1:10]
gene1 <-c(gene1,gene2)
}
gene <-unique(gene1)
sub_OS_data1 <-sub_OS[["RNA"]]@data[gene,]
sub_OS_data <-t(scale(t(as.matrix(sub_OS_data1))))
type <-unique(sub_OS@meta.data[,"seurat_clusters"])
bb <-NULL
for(i in 1: length(type)){
ind <-which(sub_OS@meta.data[,"seurat_clusters"]==type[i])
aa <-rowMeans(as.matrix(sub_OS_data[,ind]))
bb <-cbind(bb,aa)
}
colnames(bb) <-type
library(reshape2)
library(ggplot2)
bb[bb>3.8] <-3.8
data_melt<-melt(bb)
data_melt$Var2 <-factor(data_melt$Var2,levels=c("0","1","2","3","4"))
p <- ggplot(data = data_melt) + geom_tile(aes(x=Var2,y=Var1, fill = value)) +theme_classic() +  theme(axis.ticks = element_blank(), axis.line = element_blank()) + xlab('row name') + ylab('column name')
p1 <- p+scale_fill_gradient2('z-score',low = 'blue', high = 'red', mid = 'white') 
pdf("heatmap_gene_celltype.pdf",width=4,height=8)			
p1
dev.off()

setwd("./Analysis/Osteoclasts/gene")
gene <-c("CTSK","ACP5","CD74", "CD14","TOP2A")
library(Seurat)
DefaultAssay(sub_OS) <- "RNA"
plots <- VlnPlot(sub_OS, features =gene, group.by = "seurat_clusters", pt.size = 0, combine = FALSE)
pdf("Gene.pdf",width=7,height=10)
wrap_plots(plots = plots, ncol = 1)
dev.off() 
##
pdf("Tsne_CTSK.pdf",width=20,height=10)
print(FeaturePlot(sub_OS, features =c("CTSK","ACP5"), slot = "data",cols = c("grey", "red"),blend.threshold = 1,reduction = "umap",pt.size=3,label.size = 10))
dev.off()
pdf("Tsne_CD74.pdf",width=20,height=10)
print(FeaturePlot(sub_OS, features =c("CD74","CD14"), slot = "data",cols = c("grey", "red"),blend.threshold = 1,reduction = "umap",pt.size=3,))
dev.off()
##
sub_OS@meta.data$donor_status2 <-factor(sub_OS@meta.data$donor_status2,levels=c("Tumor","Chemo","Re_Me"))
pdf("Vlo_CD74.pdf",width=10,height=5)
print(VlnPlot(sub_OS, features =c("CD74"), split.by = "donor_status2", group.by = "seurat_clusters", pt.size = 0, combine = FALSE))
dev.off()
pdf("Vlo_CD14.pdf",width=10,height=5)
print(VlnPlot(sub_OS, features =c("CD14"), split.by = "donor_status2", group.by = "seurat_clusters", pt.size = 0, combine = FALSE))
dev.off()
####
marker <-read.table("./Analysis/Osteoclasts/markers_r0.2.txt")
library(ggplot2)
library(clusterProfiler)
library(org.Hs.eg.db)
for(i in 0:4){
ind <-which(marker[,"cluster"]==i)
gene <-marker[ind,"gene"]
up <-gene
gs = bitr(up, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
head(gs)
ego.bp = enrichGO(gene=gs$ENTREZID, OrgDb = org.Hs.eg.db,ont= "BP",pAdjustMethod = "BH",pvalueCutoff= 0.05,qvalueCutoff= 0.2,readable= TRUE)
head(ego.bp)
write.csv(ego.bp, file =paste(i,"_go.csv"))
pdf(paste(i,"_go.pdf"))
print(dotplot(ego.bp, showCategory=10,title="up gene GoTerm"))
dev.off()
pdf(paste(i,"_go_box.pdf"))
print(barplot(ego.bp, showCategory=10,title="up gene GoTerm"))
dev.off()
kk <- enrichKEGG(gene= gs$ENTREZID, organism= 'hsa',pvalueCutoff = 0.05)
write.csv(kk,file =paste(i,"_kegg.csv"))
pdf(paste(i,"_kegg.pdf"))
print(dotplot(kk, showCategory=10,title="KEGG_Up_biological"))
dev.off()
}
##
data <-as.data.frame(read.csv("0 _go.csv",sep=",",header=T))
data <-data[c(48,50,52,54,58),]
p <-ggplot(data,aes(x=Description,y=Count,fill=pvalue)) + 
  geom_bar(stat='identity',color='black',width = 0.65) +
  coord_flip() +    
  scale_fill_gradient(low='#FFFFF0', high='darkgoldenrod1')  
mytheme <- theme(axis.title=element_text(face="bold", size=10,colour = 'gray25'), 
                 axis.text=element_text(face="bold", size=10,colour = 'gray25'), 
                 axis.line = element_line(size=0.5, colour = 'black'), 
                 axis.line.y = element_blank(),  
                 axis.ticks.y = element_blank(), 
                 panel.background = element_rect(fill="white"), 
                 panel.grid.major.y=element_blank(), 
                 panel.grid.minor.y=element_blank(), 
                 panel.grid.minor.x=element_blank()) 
pdf("c0_barplot.pdf",width=15,height=5)
p + labs(x='',y='Counts',fill='p value') 
dev.off()
#
data <-as.data.frame(read.csv("1 _go.csv",sep=",",header=T))
data <-data[c(1,3,4,5,6),]
p <-ggplot(data,aes(x=Description,y=Count,fill=pvalue)) + 
  geom_bar(stat='identity',color='black',width = 0.65) +
  coord_flip() +    
  scale_fill_gradient(low='#FFFFF0', high='darkgoldenrod1')  
mytheme <- theme(axis.title=element_text(face="bold", size=10,colour = 'gray25'), 
                 axis.text=element_text(face="bold", size=10,colour = 'gray25'), 
                 axis.line = element_line(size=0.5, colour = 'black'), 
                 axis.line.y = element_blank(),  
                 axis.ticks.y = element_blank(), 
                 panel.background = element_rect(fill="white"), 
                 panel.grid.major.y=element_blank(), 
                 panel.grid.minor.y=element_blank(), 
                 panel.grid.minor.x=element_blank()) 
pdf("c1_barplot.pdf",width=15,height=5)
p + labs(x='',y='Counts',fill='p value') 
dev.off()
#
data <-as.data.frame(read.csv("2 _go.csv",sep=",",header=T))
data <-data[c(2,3,4,5,7),]
p <-ggplot(data,aes(x=Description,y=Count,fill=pvalue)) + 
  geom_bar(stat='identity',color='black',width = 0.65) +
  coord_flip() +    
  scale_fill_gradient(low='#FFFFF0', high='darkgoldenrod1')  
mytheme <- theme(axis.title=element_text(face="bold", size=10,colour = 'gray25'), 
                 axis.text=element_text(face="bold", size=10,colour = 'gray25'), 
                 axis.line = element_line(size=0.5, colour = 'black'), 
                 axis.line.y = element_blank(),  
                 axis.ticks.y = element_blank(), 
                 panel.background = element_rect(fill="white"), 
                 panel.grid.major.y=element_blank(), 
                 panel.grid.minor.y=element_blank(), 
                 panel.grid.minor.x=element_blank()) 
pdf("c2_barplot.pdf",width=15,height=5)
p + labs(x='',y='Counts',fill='p value') 
dev.off()
#
data <-as.data.frame(read.csv("3 _go.csv",sep=",",header=T))
data <-data[c(3,10,12,21,27),]
p <-ggplot(data,aes(x=Description,y=Count,fill=pvalue)) + 
  geom_bar(stat='identity',color='black',width = 0.65) +
  coord_flip() +    
  scale_fill_gradient(low='#FFFFF0', high='darkgoldenrod1')  
mytheme <- theme(axis.title=element_text(face="bold", size=10,colour = 'gray25'), 
                 axis.text=element_text(face="bold", size=10,colour = 'gray25'), 
                 axis.line = element_line(size=0.5, colour = 'black'), 
                 axis.line.y = element_blank(),  
                 axis.ticks.y = element_blank(), 
                 panel.background = element_rect(fill="white"), 
                 panel.grid.major.y=element_blank(), 
                 panel.grid.minor.y=element_blank(), 
                 panel.grid.minor.x=element_blank()) 
pdf("c3_barplot.pdf",width=15,height=5)
p + labs(x='',y='Counts',fill='p value') 
dev.off()
#
data <-as.data.frame(read.csv("4 _go.csv",sep=",",header=T))
data <-data[c(1,3,4,15,18),]
p <-ggplot(data,aes(x=Description,y=Count,fill=pvalue)) + 
  geom_bar(stat='identity',color='black',width = 0.65) +
  coord_flip() +    
  scale_fill_gradient(low='#FFFFF0', high='darkgoldenrod1')  
mytheme <- theme(axis.title=element_text(face="bold", size=10,colour = 'gray25'), 
                 axis.text=element_text(face="bold", size=10,colour = 'gray25'), 
                 axis.line = element_line(size=0.5, colour = 'black'), 
                 axis.line.y = element_blank(),  
                 axis.ticks.y = element_blank(), 
                 panel.background = element_rect(fill="white"), 
                 panel.grid.major.y=element_blank(), 
                 panel.grid.minor.y=element_blank(), 
                 panel.grid.minor.x=element_blank()) 
pdf("c4_barplot.pdf",width=15,height=5)
p + labs(x='',y='Counts',fill='p value') 
dev.off()
##########比例
bb <-NULL
sample <- unique(sub_OS@meta.data[,"donor_ID"])
for(i in 1:length(sample)){
aa <-table(sub_OS@meta.data[which(sub_OS@meta.data[,"donor_ID"]==sample[i]),"seurat_clusters"])
bb <-cbind(bb,aa)
}
colnames(bb) <-sample
cc <-NULL
for(i in 1:5){
aa <-bb[i,]/sum(bb[i,])
cc <-rbind(cc,aa)
}
rownames(cc) <-rownames(bb)
library(reshape2)
library(ggplot2)
library(RColorBrewer)
colorr <- brewer.pal(4,"Set2")
df <-melt(cc)
df$Var2 <-factor(df$Var2,levels=colnames(cc))
df$Var1 <-factor(df$Var1,levels=c("0","1","2","3","4"))
pdf(file="box_precent_Donor.pdf",width=20,height=10)
ggplot(df,aes(x=Var1,y=value,fill=Var2))+scale_fill_manual(values=my36colors[1:15])+geom_bar(stat="identity")+theme(legend.title=element_blank()) 
dev.off()
##############
##
bb <-NULL
sample <- unique(sub_OS@meta.data[,"donor_status1"])
for(i in 1:length(sample)){
aa <-table(sub_OS@meta.data[which(sub_OS@meta.data[,"donor_status1"]==sample[i]),"seurat_clusters"])
bb <-cbind(bb,aa)
}
colnames(bb) <-sample
cc <-NULL
for(i in 1:5){
aa <-bb[i,]/sum(bb[i,])
cc <-rbind(cc,aa)
}
rownames(cc) <-rownames(bb)
library(reshape2)
library(ggplot2)
library(RColorBrewer)
colorr <- brewer.pal(4,"Set2")
df <-melt(cc)
df$Var2 <-factor(df$Var2,levels=colnames(cc))
df$Var1 <-factor(df$Var1,levels=c("0","1","2","3","4"))
pdf(file="box_precent_Stage.pdf",width=20,height=10)
ggplot(df,aes(x=Var1,y=value,fill=Var2))+scale_fill_manual(values=my36colors[1:3])+geom_bar(stat="identity")+theme(legend.title=element_blank()) 
dev.off()
##
bb <-NULL
sample <- unique(sub_OS@meta.data[,"donor_status2"])
for(i in 1:length(sample)){
aa <-table(sub_OS@meta.data[which(sub_OS@meta.data[,"donor_status2"]==sample[i]),"seurat_clusters"])
bb <-cbind(bb,aa)
}
colnames(bb) <-sample
cc <-NULL
for(i in 1:5){
aa <-bb[i,]/sum(bb[i,])
cc <-rbind(cc,aa)
}
rownames(cc) <-rownames(bb)
library(reshape2)
library(ggplot2)
library(RColorBrewer)
colorr <- brewer.pal(4,"Set2")
df <-melt(cc)
df$Var2 <-factor(df$Var2,levels=colnames(cc))
df$Var1 <-factor(df$Var1,levels=c("0","1","2","3","4"))
pdf(file="box_precent_Stage1.pdf",width=20,height=10)
ggplot(df,aes(x=Var1,y=value,fill=Var2))+scale_fill_manual(values=my36colors[1:4])+geom_bar(stat="identity")+theme(legend.title=element_blank()) 
dev.off()
###
library(reshape2)
library(ggplot2)
library(RColorBrewer)
dat <- sub_OS@meta.data[,c("donor_status2","seurat_clusters","donor_ID")]
dat1 <- dat[which(dat$donor_status2 =="Tumor"),]
bb <-NULL
sample <- unique(dat1[,"donor_ID"])
for(i in 1:length(sample)){
aa <-table(dat1[which(dat1[,"donor_ID"]==sample[i]),"seurat_clusters"])
bb <-cbind(bb,aa)}
colnames(bb) <-sample
cc <-t(t(bb)/as.vector(table(OS$donor_ID)[sample]))
dd1 <-melt(cc)
dd1$status <-"Tumor"
dat1 <- dat[which(dat$donor_status2 =="Chemo"),]
bb <-NULL
sample <- unique(dat1[,"donor_ID"])
for(i in 1:length(sample)){
aa <-table(dat1[which(dat1[,"donor_ID"]==sample[i]),"seurat_clusters"])
bb <-cbind(bb,aa)}
colnames(bb) <-sample
cc <-t(t(bb)/as.vector(table(OS$donor_ID)[sample]))
dd2 <-melt(cc)
dd2$status <-"Chemo"
dat1 <- dat[which(dat$donor_status2 =="Re_Me"),]
bb <-NULL
sample <- unique(dat1[,"donor_ID"])
for(i in 1:length(sample)){
aa <-table(dat1[which(dat1[,"donor_ID"]==sample[i]),"seurat_clusters"])
bb <-cbind(bb,aa)}
colnames(bb) <-sample
cc <-t(t(bb)/as.vector(table(OS$donor_ID)[sample]))
dd3 <-melt(cc)
dd3$status <-"Re_Me"

#
df <-rbind(dd1,dd2,dd3)
levels <-c("0","1","2","3","4")
df$status = factor(df$status, levels=c("Tumor","Chemo","Re_Me"))
dodge <- position_dodge(width = 0.8)
pdf(file="box_precent_dot_stage2.pdf",width=10,height=4)
ggplot(df, aes(x=factor(Var1,levels=levels), y=value,fill = status))  + stat_boxplot(geom="errorbar",width=0.3,position = dodge,size=2)+geom_boxplot(width=0.8,outlier.colour = NA)+
ylim(0,0.15)+theme(legend.title=element_blank())+theme(axis.text.x = element_text(size = 12,angle = 45, hjust = 1, vjust = 1, face = "bold"))+labs(x="Cluster", y = "overlap")
dev.off()
#############t.test/wilcox
bb <-NULL
for(i in c(0:4)){
aa <- c(wilcox.test(df[which(df[,1]==i),][which(df[which(df[,1]==i),][,4]=="Tumor"),3],df[which(df[,1]==i),][which(df[which(df[,1]==i),][,4]=="Chemo"),3])$p.value,
wilcox.test(df[which(df[,1]==i),][which(df[which(df[,1]==i),][,4]=="Tumor"),3],df[which(df[,1]==i),][which(df[which(df[,1]==i),][,4]=="Re_Me"),3])$p.value,
wilcox.test(df[which(df[,1]==i),][which(df[which(df[,1]==i),][,4]=="Chemo"),3],df[which(df[,1]==i),][which(df[which(df[,1]==i),][,4]=="Re_Me"),3])$p.value)
bb <-rbind(bb,aa)
}
bb
aa 0.93722944 1.0000000 0.8571429
aa 0.06493506 0.6428571 0.1428571
aa 0.24025974 0.6428571 0.2857143
aa 0.93504305 0.8643948 0.8643948
aa 0.29710698 1.0000000 0.6428571
#######
sub_OS@active.ident <-factor(as.matrix(sub_OS@meta.data)[,"donor_status2"])
DefaultAssay(sub_OS) <- "RNA"
cluster.markers <- FindMarkers(sub_OS, ident.1 = c("Chemo"), ident.2 =c("Tumor"), min.pct = 0.1,logfc.threshold=0)
filt_data <-cluster.markers
data<-c()
for(i in 1:nrow(filt_data)){
	if(filt_data[i,1]<0.05){
		if(filt_data[i,2]>0.25){
			data<-rbind(data,c(filt_data[i,2],filt_data[i,1],"up"))
		}
		else if(filt_data[i,2]<(-0.25)){
			data<-rbind(data,c(filt_data[i,2],filt_data[i,1],"down"))
		}
		else{
		    data<-rbind(data,c(filt_data[i,2],filt_data[i,1],"no"))
	    }
	}	
	else{
		data<-rbind(data,c(filt_data[i,2],filt_data[i,1],"no"))
	}
}
table(data[,3])
DEG_num=as.numeric(table(data[,3])[1]+table(data[,3])[3])
sub_title=paste("DEGs:",DEG_num,sep="")
colnames(data)<-c("logFC","pvalue",sub_title)
rownames(data)<-rownames(filt_data)
write.table(data,file="Cluster_volcano_plot_format.txt",quote=F,sep="\t")
data <-read.table("Cluster_volcano_plot_format.txt",header=1,sep='\t')
up_num=length(intersect(which(data[,2]<0.05),which(data[,1]>0.25)))
down_num=length(intersect(which(data[,2]<0.05),which(data[,1]<(-0.25))))
legend1=paste("up regulated:",up_num,sep=" ")
legend2=paste("down regulated:",down_num,sep=" ")
index <-which(data[,3]=="no")
index1 <-which(data[,3]=="up")
index2 <-which(data[,3]=="down")
ind <-c(index1,index2,index)
data <-data[ind,]
library("ggplot2")
r03=ggplot(data,aes(x=logFC,y=-log10(pvalue))) +geom_point()+ theme_bw()+theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())
r03=r03 + geom_point(aes(color =DEGs.694))
volcano = r03 +scale_color_manual(values = c("#009966","#0066CC","#FF0033"),breaks=c("up","down"),labels=c(legend1,legend2))+theme(legend.position = c(.25, .75))
volcano=volcano+geom_hline(yintercept=-log10(0.05),linetype=2,size=1)+geom_vline(xintercept=c((-0.25),0.25),linetype=2,size=1)
library(psych)
library(Seurat)
library(dplyr)
library(ggrepel)
data$sign <- ifelse(data$logFC >= 1 | data$logFC <= -1,rownames(data),NA)
#data$sign <- ifelse( -log10(data$pvalue) >= 10,rownames(data),NA)
pdf(file="volcano_Chemo_Tumor.pdf",width=7,height=5)	
volcano +geom_text_repel(label=data$sign,colour="black",size=2,box.padding = unit(0.3, "lines"), point.padding = unit(0.4, "lines"), show.legend = F) 
dev.off()
####
library(clusterProfiler)
library(org.Hs.eg.db)
data <-read.table("Cluster_volcano_plot_format.txt",header=1,sep='\t')
up <-rownames(data[intersect(which(data[,2]<0.05),which(data[,1]>0.25)),])
down <-rownames(data[intersect(which(data[,2]<0.05),which(data[,1]<(-0.25))),])
gs <- up
gs = bitr(gs, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
head(gs)
ego.bp = enrichGO(gene =  gs$ENTREZID,OrgDb = org.Hs.eg.db,ont  = "BP",pAdjustMethod = "BH",pvalueCutoff  = 0.05,qvalueCutoff  = 0.2,readable = TRUE)					
head(ego.bp)
# 
nam = paste("Chemo_Tumor_UP","Go.csv",sep ="")
write.csv(ego.bp, file =nam)
nam = paste("Chemo_Tumor_UP","Go.pdf",sep ="")
pdf(nam,width=8,height=6)
print(barplot(ego.bp, showCategory=9,title="GO_biological",drop=T))
dev.off()
#
gs <- down
gs = bitr(gs, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
head(gs)
ego.bp = enrichGO(gene =  gs$ENTREZID,OrgDb = org.Hs.eg.db,ont  = "BP",pAdjustMethod = "BH",pvalueCutoff  = 0.05,qvalueCutoff  = 0.2,readable = TRUE)					
# 
nam = paste("Chemo_Tumor_DOWN","Go.csv",sep ="")
write.csv(ego.bp, file =nam)
nam = paste("Chemo_Tumor_DOWN","Go.pdf",sep ="")
pdf(nam,width=12,height=6)
print(barplot(ego.bp, showCategory=10,title="GO_biological",drop=T))
dev.off()
##
bp1 <-read.csv("Chemo_Tumor_UPGo.csv")
bp2 <-read.csv("Chemo_Tumor_DOWNGo.csv")
bp1$celltype <-"UP"
bp2$celltype <-"DOWN"
bp <-rbind(bp1[c(1,2,3,4,5),],bp2[c(3,4,6,10,18),])
levels <-rev(bp$Description)
pdf("UP_DOWN_GO_function.pdf",width=10,height=11)
ggplot(bp,aes(x=factor(Description,levels=levels),y=Count,fill=pvalue)) + 
  geom_bar(stat='identity',color='black',width = 0.65) +
  coord_flip() +    
  scale_fill_gradient(low='#FFFFF0', high='darkgoldenrod1')
dev.off()
##
nodeid.tbl_tree <- utils::getFromNamespace("nodeid.tbl_tree", "tidytree")
rootnode.tbl_tree <- utils::getFromNamespace("rootnode.tbl_tree", "tidytree")
offspring.tbl_tree <- utils::getFromNamespace("offspring.tbl_tree", "tidytree")
child.tbl_tree <- utils::getFromNamespace("child.tbl_tree", "tidytree")
parent.tbl_tree <- utils::getFromNamespace("parent.tbl_tree", "tidytree")
offspring.tbl_tree_item <- function(.data, .node, tiponly = FALSE, self_include = FALSE, ...) {
    x <- child.tbl_tree(.data, .node)

    ## https://github.com/GuangchuangYu/ggtree/issues/239
    rn <- rootnode.tbl_tree(.data)$node
    x <- x[x$node != rn, ]

    if (nrow(x) == 0) {
        if (self_include) {
            x <- .data[.data$node == .node, ]
        } 

        return(x)
    }

    ## id <- x$node
    ## i <- 1
    ## while(i <= length(id)) {
    ##     id <- c(id, child(.data, id[i])$node)
    ##     i <- i + 1
    ## }
    ## filter_(.data, ~ node %in% id)

    parent <- .data$parent
    children <- .data$node
    ## n <- length(parent)
    n <- max(parent)

    kids <- vector("list", n)
    for (i in seq_along(parent)) {
        kids[[parent[i]]] <-c(kids[[parent[i]]], children[i])
    }

    id <- x$node
    i <- 1
    while(i <= length(id)) {
        id <- c(id, kids[[id[i]]])
        i <- i + 1
    }

    if (self_include) {
        id <- c(.node, id)
    }

    sp <- .data[children %in% id,]
    if (tiponly) {
        return(sp[sp$node < rn,])
    }
    return(sp)
}
library(ggplot2)
library(clusterProfiler)
library(org.Hs.eg.db)
gs = bitr(up, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
ego.bp = enrichGO(gene=gs$ENTREZID, OrgDb = org.Hs.eg.db,ont= "BP",pAdjustMethod = "BH",pvalueCutoff= 0.05,qvalueCutoff= 0.2,readable= TRUE)
library(clusterProfiler)
library(DOSE)
library(enrichplot)
library(tidytree)
edox2 <- pairwise_termsim(ego.bp)
p1 <- treeplot(edox2)
p2 <- treeplot(edox2, hclust_method = "average")
ggsave(p1,filename =paste("cnetplot_Chemo_Tumor_Up","json","_tree1.pdf"), width =12,height =10)
ggsave(p2,filename =paste("cnetplot_Chemo_Tumor_Up","json","_tree2.pdf"), width =15,height =10)
#####
gs = bitr(down, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
ego.bp = enrichGO(gene=gs$ENTREZID, OrgDb = org.Hs.eg.db,ont= "BP",pAdjustMethod = "BH",pvalueCutoff= 0.05,qvalueCutoff= 0.2,readable= TRUE)
library(clusterProfiler)
library(DOSE)
library(enrichplot)
edox2 <- pairwise_termsim(ego.bp)
p1 <- treeplot(edox2)
p2 <- treeplot(edox2, hclust_method = "average")
ggsave(p1,filename =paste("cnetplot_Chemo_Tumor_Down","json","_tree1.pdf"), width =12,height =10)
ggsave(p2,filename =paste("cnetplot_Chemo_Tumor_Down","json","_tree2.pdf"), width =15,height =10)

## Monocle analysis
setwd("./Analysis/Osteoclasts/Monocle")
library(monocle3)
scRNA <-sub_OS
data<-GetAssayData(scRNA,assay ='RNA',slot ='counts')
cell_metadata <-scRNA@meta.data
gene_annotation <-data.frame(gene_short_name =rownames(data))
rownames(gene_annotation)<-rownames(data)
cds <-new_cell_data_set(data,cell_metadata =cell_metadata,gene_metadata =gene_annotation)
cds <- preprocess_cds(cds, num_dim = 30)
cds <- align_cds(cds, alignment_group = "donor_ID")
cds <- reduce_dimension(cds)
cds <- cluster_cells(cds)
cds <- learn_graph(cds)
pdf("clu_pseu.pdf",width=7,height=5)
plot_cells(cds, color_cells_by = "seurat_clusters",label_cell_groups=FALSE,label_leaves=FALSE,label_branch_points=FALSE,graph_label_size=1.5)
dev.off()
#
cds <- order_cells(cds)
pdf("clu_pseu1",width=7,height=5)
plot_cells(cds,color_cells_by = "pseudotime",label_cell_groups=FALSE,label_leaves=FALSE,label_branch_points=FALSE,graph_label_size=1.5)
dev.off()







#####
###########
############
#########
###########Immune ITL
######
setwd("./Analysis/Tcell")
sub_OS <-subset(OS,idents= c("3","19","18","21"))
DefaultAssay(sub_OS) <- "integrated"
# Run the standard workflow for visualization and clustering
#sub_OS <- ScaleData(sub_OS, verbose = FALSE)
sub_OS <- RunPCA(sub_OS, npcs = 30, verbose = FALSE)
sub_OS <- RunUMAP(sub_OS, reduction = "pca", dims = 1:30)
sub_OS <- FindNeighbors(sub_OS, reduction = "pca", dims = 1:30)
sub_OS <- FindClusters(sub_OS, resolution = 0.15)
p0 <- DimPlot(object = sub_OS, reduction = "umap",label=TRUE)
p1<- DimPlot(object = sub_OS, reduction = "umap",group.by ="donor_ID",label=TRUE)
p2 <- DimPlot(object = sub_OS, reduction = "umap",group.by ="donor_status",label=TRUE)
pdf("Umap__OS.pdf",width=22,height=20)
CombinePlots(plots  = list(p0, p1, p2))
dev.off()
saveRDS(sub_OS, file="data_rPCA_sub_OS.rds")
##sub_OS <-readRDS("data_rPCA_sub_OS.rds")
DefaultAssay(sub_OS) <- "RNA"
OS.markers <- FindAllMarkers(object = sub_OS, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.table(OS.markers,"markers_r0.15.txt",row.names=TRUE,col.names=TRUE)
##
gene <-c("CD3D","CD3E","CD8A","CD8B","GZMK","GZMA","NKG7","GNLY","CD4","TIGIT","CD79A","CD79B","IGKC","MS4A1","TPSB2","CPA3","TCF4","IRF8","MKI67","TOP2A") 
library(Seurat)
DefaultAssay(sub_OS) <- "RNA"
sub_OS$seurat_clusters <- factor(sub_OS$seurat_clusters, levels = c("6","1","2","8","3","4","5","9","7","0"))
plots <- VlnPlot(sub_OS, features =gene,  group.by = "seurat_clusters", pt.size = 0, combine = FALSE)
pdf("Gene.pdf",width=10,height=30)
wrap_plots(plots = plots, ncol = 1)
dev.off() 


#######
UmapPlot<- function(obj) {
  library(ggplot2)
  library(dplyr)
  UMAP = obj@reductions$umap@cell.embeddings %>% as.data.frame() %>% cbind(cell_type = obj@meta.data$seurat_clusters) 
  allcolour=c("#DC143C","#0000FF","#20B2AA","#FFA500","#9370DB","#98FB98","#F08080","#1E90FF","#7CFC00","#FFFF00","#808000","#FF00FF","#FA8072","#7B68EE","#9400D3","#800080","#A0522D","#D2B48C","#D2691E","#87CEEB","#40E0D0","#5F9EA0","#FF1493","#0000CD","#008B8B","#FFE4B5","#8A2BE2","#228B22","#E9967A","#4682B4","#32CD32","#F0E68C","#FFFFE0","#EE82EE","#FF6347","#6A5ACD","#9932CC","#8B008B","#8B4513","#DEB887")
  my36colors <-c('#E5D2DD', '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87','#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#8C549C', '#585658', '#9FA3A8', '#E0D4CA', '#5F3D69', '#C5DEBA', '#58A4C3', '#E4C755', '#F7F398','#AA9A59', '#E63863', '#E39A35', '#C1E6F3', '#6778AE', '#91D0BE', '#B53E2B','#712820', '#DCC1DD', '#CCE0F5',  '#CCC9E6', '#625D9E', '#68A180', '#3A6963','#968175')
  p <- ggplot(UMAP,aes(x= UMAP_1, y = UMAP_2 ,color = cell_type)) +  geom_point(size=0.2, alpha =1 )  +  scale_color_manual(values = my36colors)
  p2 <- p  +theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),  panel.border = element_blank(),  axis.title = element_blank(),  axis.text = element_blank(), axis.ticks = element_blank(),panel.background = element_rect(fill = 'white'), plot.background=element_rect(fill="white"))
  p3 <- p2 + theme(legend.title = element_blank(),legend.key=element_rect(fill='white'), legend.text = element_text(size=20),legend.key.size=unit(1,'cm') ) + guides(color = guide_legend(override.aes = list(size=5))) 
#  p4 <- p3 + geom_segment(aes(x = min(TSNE$tSNE_1) , y = min(TSNE$tSNE_2) ,xend = min(TSNE$tSNE_1) +3, yend = min(TSNE$tSNE_2) ),colour = "black", size=1,arrow = arrow(length = unit(0.3,"cm")))+ geom_segment(aes(x = min(TSNE$tSNE_1), y = min(TSNE$tSNE_2),xend = min(TSNE$tSNE_1) , yend = min(TSNE$tSNE_2) + 3),colour = "black", size=1,arrow = arrow(length = unit(0.3,"cm"))) + annotate("text", x = min(TSNE$tSNE_1) +1.5, y = min(TSNE$tSNE_2) -1, label = "tSNE_1",color="black",size = 3, fontface="bold" ) + annotate("text", x = min(TSNE$tSNE_1) -1, y = min(TSNE$tSNE_2) + 1.5, label = "tSNE_2",color="black",size = 3, fontface="bold" ,angle=90) 
  cell_type_med <- UMAP %>%group_by(cell_type) %>%summarise(UMAP_1 = median(UMAP_1),UMAP_2 = median(UMAP_2))
  library(ggrepel)
  p5 <- p3 + geom_label_repel(aes(label=cell_type), fontface="bold",data = cell_type_med, point.padding=unit(1.5, "lines")) + theme(legend.position = "none")                  
  return(p5)
}
pdf("UMAP_Clusters.pdf",width=4,height=4)
UmapPlot(sub_OS)
dev.off()

UmapPlot<- function(obj) {
  library(ggplot2)
  library(dplyr)
  UMAP = obj@reductions$umap@cell.embeddings %>% as.data.frame() %>% cbind(cell_type = obj@meta.data$donor_ID) 
  allcolour=c("#DC143C","#0000FF","#20B2AA","#FFA500","#9370DB","#98FB98","#F08080","#1E90FF","#7CFC00","#FFFF00","#808000","#FF00FF","#FA8072","#7B68EE","#9400D3","#800080","#A0522D","#D2B48C","#D2691E","#87CEEB","#40E0D0","#5F9EA0","#FF1493","#0000CD","#008B8B","#FFE4B5","#8A2BE2","#228B22","#E9967A","#4682B4","#32CD32","#F0E68C","#FFFFE0","#EE82EE","#FF6347","#6A5ACD","#9932CC","#8B008B","#8B4513","#DEB887")
  my36colors <-c('#E5D2DD', '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87','#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#8C549C', '#585658', '#9FA3A8', '#E0D4CA', '#5F3D69', '#C5DEBA', '#58A4C3', '#E4C755', '#F7F398','#AA9A59', '#E63863', '#E39A35', '#C1E6F3', '#6778AE', '#91D0BE', '#B53E2B','#712820', '#DCC1DD', '#CCE0F5',  '#CCC9E6', '#625D9E', '#68A180', '#3A6963','#968175')
  p <- ggplot(UMAP,aes(x= UMAP_1, y = UMAP_2 ,color = cell_type)) +  geom_point(size=0.2, alpha =1 )  +  scale_color_manual(values = my36colors)
  p2 <- p  +theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),  panel.border = element_blank(),  axis.title = element_blank(),  axis.text = element_blank(), axis.ticks = element_blank(),panel.background = element_rect(fill = 'white'), plot.background=element_rect(fill="white"))
  p3 <- p2 + theme(legend.title = element_blank(),legend.key=element_rect(fill='white'), legend.text = element_text(size=20),legend.key.size=unit(1,'cm') ) + guides(color = guide_legend(override.aes = list(size=5))) 
#  p4 <- p3 + geom_segment(aes(x = min(TSNE$tSNE_1) , y = min(TSNE$tSNE_2) ,xend = min(TSNE$tSNE_1) +3, yend = min(TSNE$tSNE_2) ),colour = "black", size=1,arrow = arrow(length = unit(0.3,"cm")))+ geom_segment(aes(x = min(TSNE$tSNE_1), y = min(TSNE$tSNE_2),xend = min(TSNE$tSNE_1) , yend = min(TSNE$tSNE_2) + 3),colour = "black", size=1,arrow = arrow(length = unit(0.3,"cm"))) + annotate("text", x = min(TSNE$tSNE_1) +1.5, y = min(TSNE$tSNE_2) -1, label = "tSNE_1",color="black",size = 3, fontface="bold" ) + annotate("text", x = min(TSNE$tSNE_1) -1, y = min(TSNE$tSNE_2) + 1.5, label = "tSNE_2",color="black",size = 3, fontface="bold" ,angle=90) 
  cell_type_med <- UMAP %>%group_by(cell_type) %>%summarise(UMAP_1 = median(UMAP_1),UMAP_2 = median(UMAP_2))
  library(ggrepel)
  p5 <- p3 + geom_label_repel(aes(label=cell_type), fontface="bold",data = cell_type_med, point.padding=unit(1.5, "lines")) + theme(legend.position = "none")                  
  return(p5)
}
pdf("UMAP_Donor.pdf",width=10,height=10)
UmapPlot(sub_OS)
dev.off()

UmapPlot<- function(obj) {
  library(ggplot2)
  library(dplyr)
  UMAP = obj@reductions$umap@cell.embeddings %>% as.data.frame() %>% cbind(cell_type = obj@meta.data$donor_status1) 
  allcolour=c("#DC143C","#0000FF","#20B2AA","#FFA500","#9370DB","#98FB98","#F08080","#1E90FF","#7CFC00","#FFFF00","#808000","#FF00FF","#FA8072","#7B68EE","#9400D3","#800080","#A0522D","#D2B48C","#D2691E","#87CEEB","#40E0D0","#5F9EA0","#FF1493","#0000CD","#008B8B","#FFE4B5","#8A2BE2","#228B22","#E9967A","#4682B4","#32CD32","#F0E68C","#FFFFE0","#EE82EE","#FF6347","#6A5ACD","#9932CC","#8B008B","#8B4513","#DEB887")
  my36colors <-c('#E5D2DD', '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87','#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#8C549C', '#585658', '#9FA3A8', '#E0D4CA', '#5F3D69', '#C5DEBA', '#58A4C3', '#E4C755', '#F7F398','#AA9A59', '#E63863', '#E39A35', '#C1E6F3', '#6778AE', '#91D0BE', '#B53E2B','#712820', '#DCC1DD', '#CCE0F5',  '#CCC9E6', '#625D9E', '#68A180', '#3A6963','#968175')
  p <- ggplot(UMAP,aes(x= UMAP_1, y = UMAP_2 ,color = cell_type)) +  geom_point(size=0.2, alpha =1 )  +  scale_color_manual(values = my36colors)
  p2 <- p  +theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),  panel.border = element_blank(),  axis.title = element_blank(),  axis.text = element_blank(), axis.ticks = element_blank(),panel.background = element_rect(fill = 'white'), plot.background=element_rect(fill="white"))
  p3 <- p2 + theme(legend.title = element_blank(),legend.key=element_rect(fill='white'), legend.text = element_text(size=20),legend.key.size=unit(1,'cm') ) + guides(color = guide_legend(override.aes = list(size=5))) 
#  p4 <- p3 + geom_segment(aes(x = min(TSNE$tSNE_1) , y = min(TSNE$tSNE_2) ,xend = min(TSNE$tSNE_1) +3, yend = min(TSNE$tSNE_2) ),colour = "black", size=1,arrow = arrow(length = unit(0.3,"cm")))+ geom_segment(aes(x = min(TSNE$tSNE_1), y = min(TSNE$tSNE_2),xend = min(TSNE$tSNE_1) , yend = min(TSNE$tSNE_2) + 3),colour = "black", size=1,arrow = arrow(length = unit(0.3,"cm"))) + annotate("text", x = min(TSNE$tSNE_1) +1.5, y = min(TSNE$tSNE_2) -1, label = "tSNE_1",color="black",size = 3, fontface="bold" ) + annotate("text", x = min(TSNE$tSNE_1) -1, y = min(TSNE$tSNE_2) + 1.5, label = "tSNE_2",color="black",size = 3, fontface="bold" ,angle=90) 
  cell_type_med <- UMAP %>%group_by(cell_type) %>%summarise(UMAP_1 = median(UMAP_1),UMAP_2 = median(UMAP_2))
  library(ggrepel)
  p5 <- p3 + geom_label_repel(aes(label=cell_type), fontface="bold",data = cell_type_med, point.padding=unit(1.5, "lines")) + theme(legend.position = "none")                  
  return(p5)
}
pdf("UMAP_status.pdf",width=4,height=4)
UmapPlot(sub_OS)
dev.off()

####比例
######
bb <-NULL
sample <- unique(sub_OS@meta.data[,"donor_ID"])
for(i in 1:length(sample)){
aa <-table(sub_OS@meta.data[which(sub_OS@meta.data[,"donor_ID"]==sample[i]),"seurat_clusters"])
bb <-cbind(bb,aa)
}
colnames(bb) <-sample
cc <-NULL
for(i in 1:10){
aa <-bb[i,]/sum(bb[i,])
cc <-rbind(cc,aa)
}
rownames(cc) <-rownames(bb)
library(reshape2)
library(ggplot2)
library(RColorBrewer)
colorr <- brewer.pal(4,"Set2")
df <-melt(cc)
df$Var2 <-factor(df$Var2,levels=colnames(cc))
df$Var1 <-factor(df$Var1,levels=c("6","1","2","8","3","4","5","9","7","0"))
pdf(file="box_precent_donor.pdf",width=20,height=10)
ggplot(df,aes(x=Var1,y=value,fill=Var2))+scale_fill_manual(values=my36colors[1:17])+geom_bar(stat="identity")+theme(legend.title=element_blank()) 
dev.off()
##
bb <-NULL
sample <- unique(sub_OS@meta.data[,"donor_status1"])
for(i in 1:length(sample)){
aa <-table(sub_OS@meta.data[which(sub_OS@meta.data[,"donor_status1"]==sample[i]),"seurat_clusters"])
bb <-cbind(bb,aa)
}
colnames(bb) <-sample
cc <-NULL
for(i in 1:10){
aa <-bb[i,]/sum(bb[i,])
cc <-rbind(cc,aa)
}
rownames(cc) <-rownames(bb)
library(reshape2)
library(ggplot2)
library(RColorBrewer)
colorr <- brewer.pal(4,"Set2")
df <-melt(cc)
df$Var2 <-factor(df$Var2,levels=colnames(cc))
df$Var1 <-factor(df$Var1,levels=c("6","1","2","8","3","4","5","9","7","0"))
pdf(file="box_precent_status.pdf",width=20,height=10)
ggplot(df,aes(x=Var1,y=value,fill=Var2))+scale_fill_manual(values=my36colors[1:4])+geom_bar(stat="identity")+theme(legend.title=element_blank()) 
dev.off()
##
bb <-NULL
sample <- unique(sub_OS@meta.data[,"donor_status2"])
for(i in 1:length(sample)){
aa <-table(sub_OS@meta.data[which(sub_OS@meta.data[,"donor_status2"]==sample[i]),"seurat_clusters"])
bb <-cbind(bb,aa)
}
colnames(bb) <-sample
cc <-NULL
for(i in 1:10){
aa <-bb[i,]/sum(bb[i,])
cc <-rbind(cc,aa)
}
rownames(cc) <-rownames(bb)
library(reshape2)
library(ggplot2)
library(RColorBrewer)
colorr <- brewer.pal(4,"Set2")
df <-melt(cc)
df$Var2 <-factor(df$Var2,levels=colnames(cc))
df$Var1 <-factor(df$Var1,levels=c("6","1","2","8","3","4","5","9","7","0"))
pdf(file="box_precent_status1.pdf",width=20,height=10)
ggplot(df,aes(x=Var1,y=value,fill=Var2))+scale_fill_manual(values=my36colors[1:3])+geom_bar(stat="identity")+theme(legend.title=element_blank()) 
dev.off()

###
library(reshape2)
library(ggplot2)
library(RColorBrewer)
dat <- sub_OS@meta.data[,c("donor_status2","seurat_clusters","donor_ID")]
dat1 <- dat[which(dat$donor_status2 =="Tumor"),]
bb <-NULL
sample <- unique(dat1[,"donor_ID"])
for(i in 1:length(sample)){
aa <-table(dat1[which(dat1[,"donor_ID"]==sample[i]),"seurat_clusters"])
bb <-cbind(bb,aa)}
colnames(bb) <-sample
cc <-t(t(bb)/as.vector(table(OS$donor_ID)[sample]))
dd1 <-melt(cc)
dd1$status <-"Tumor"
dat1 <- dat[which(dat$donor_status2 =="Chemo"),]
bb <-NULL
sample <- unique(dat1[,"donor_ID"])
for(i in 1:length(sample)){
aa <-table(dat1[which(dat1[,"donor_ID"]==sample[i]),"seurat_clusters"])
bb <-cbind(bb,aa)}
colnames(bb) <-sample
cc <-t(t(bb)/as.vector(table(OS$donor_ID)[sample]))
dd2 <-melt(cc)
dd2$status <-"Chemo"
#
dat1 <- dat[which(dat$donor_status2 =="Re_Me"),]
bb <-NULL
sample <- unique(dat1[,"donor_ID"])
for(i in 1:length(sample)){
aa <-table(dat1[which(dat1[,"donor_ID"]==sample[i]),"seurat_clusters"])
bb <-cbind(bb,aa)}
colnames(bb) <-sample
cc <-t(t(bb)/as.vector(table(OS$donor_ID)[sample]))
dd3 <-melt(cc)
dd3$status <-"Re_Me"

df <-rbind(dd1,dd2,dd3)
levels <- c("6","1","2","8","3","4","5","9","7","0")
df$status = factor(df$status, levels=c("Tumor","Chemo","Re_Me"))
dodge <- position_dodge(width = 0.8)
pdf(file="box_precent_dot_stage.pdf",width=10,height=4)
ggplot(df, aes(x=factor(Var1,levels=levels), y=value,fill = status))  + stat_boxplot(geom="errorbar",width=0.3,position = dodge,size=2)+geom_boxplot(width=0.8,outlier.colour = NA)+scale_fill_manual(values=my36colors[1:15])+
ylim(0,0.2)+theme(legend.title=element_blank())+theme(axis.text.x = element_text(size = 12,angle = 45, hjust = 1, vjust = 1, face = "bold"))+labs(x="Cluster", y = "overlap")
dev.off()
#############t.test/wilcox
bb <-NULL
for(i in c("6","1","2","8","3","4","5","9","7","0")){
aa <- c(wilcox.test(df[which(df[,1]==i),][which(df[which(df[,1]==i),][,4]=="Tumor"),3],df[which(df[,1]==i),][which(df[which(df[,1]==i),][,4]=="Chemo"),3])$p.value,
wilcox.test(df[which(df[,1]==i),][which(df[which(df[,1]==i),][,4]=="Tumor"),3],df[which(df[,1]==i),][which(df[which(df[,1]==i),][,4]=="Re_Me"),3])$p.value,
wilcox.test(df[which(df[,1]==i),][which(df[which(df[,1]==i),][,4]=="Chemo"),3],df[which(df[,1]==i),][which(df[which(df[,1]==i),][,4]=="Re_Me"),3])$p.value)
bb <-rbind(bb,aa)
}
bb
##
aa 0.18065268 0.11428571 1.0000000
aa 0.13752914 0.01904762 0.7878788
aa 0.10139860 0.01904762 0.6484848
aa 0.02474611 0.01141722 0.1861741
aa 0.10139860 0.17142857 0.6484848
aa 0.05247016 0.16452182 0.9210564
aa 0.13308040 0.17142857 1.0000000
aa 1.00000000 0.28495764 0.1848753
aa 0.28331448 0.04219696 0.4984039
aa 0.07342657 0.11428571 0.9272727
#######
##########

sub_OS$newcluster <-paste0(sub_OS@meta.data$seurat_clusters,sub_OS@meta.data$donor_status2)
DefaultAssay(sub_OS) <- "RNA"
sub_OS@active.ident <-factor(sub_OS$newcluster)
#gene
setwd("./Analysis/Tcell/DEG")
marker_6 <-FindMarkers(sub_OS,ident.1=c("6Chemo"),ident.2=c("6Tumor"), min.pct = 0.25, logfc.threshold = 0.25)
marker_6 <-marker_6[which(marker_6[,1]<=0.05),]
up_A <-rownames(marker_6[which(marker_6[,2]>=0),])
down_A <-rownames(marker_6[which(marker_6[,2]<=0),])
length(up_A)
length(down_A)
marker_1 <-FindMarkers(sub_OS,ident.1=c("1Chemo"),ident.2=c("1Tumor"), min.pct = 0.25, logfc.threshold = 0.25)
marker_1 <-marker_1[which(marker_1[,1]<=0.05),]
up_B <-rownames(marker_1[which(marker_1[,2]>=0),])
down_B <-rownames(marker_1[which(marker_1[,2]<=0),])
length(up_B)
length(down_B)
marker_2 <-FindMarkers(sub_OS,ident.1=c("2Chemo"),ident.2=c("2Tumor"), min.pct = 0.25, logfc.threshold = 0.25)
marker_2 <-marker_2[which(marker_2[,1]<=0.05),]
up_C <-rownames(marker_2[which(marker_2[,2]>=0),])
down_C <-rownames(marker_2[which(marker_2[,2]<=0),])
length(up_C)
length(down_C)
marker_8 <-FindMarkers(sub_OS,ident.1=c("8Chemo"),ident.2=c("8Tumor"), min.pct = 0.25, logfc.threshold = 0.25)
marker_8 <-marker_8[which(marker_8[,1]<=0.05),]
up_D <-rownames(marker_8[which(marker_8[,2]>=0),])
down_D <-rownames(marker_8[which(marker_8[,2]<=0),])
length(up_D)
length(down_D)
marker_3 <-FindMarkers(sub_OS,ident.1=c("3Chemo"),ident.2=c("3Tumor"), min.pct = 0.25, logfc.threshold = 0.25)
marker_3 <-marker_3[which(marker_3[,1]<=0.05),]
up_E <-rownames(marker_3[which(marker_3[,2]>=0),])
down_E <-rownames(marker_3[which(marker_3[,2]<=0),])
length(up_E)
length(down_E)

up <-intersect(up_A,intersect(up_B,intersect(up_C,intersect(up_D,up_E))))
down <-intersect(down_A,intersect(down_B,intersect(down_C,intersect(down_D,down_E))))
#
##Venn
library(venn)
x <- list(up = up_A, up1 = up_B, up2 = up_C,up3 = up_D,up4 = up_E)
pdf("up.pdf",width=8,height=6)
venn(x,zcolor='style', ellipse = TRUE)
dev.off()
x <- list(down = down_A, down1 = down_B, down2 = down_C,down3 = down_D,down4 = down_E)
pdf("down.pdf",width=8,height=6)
venn(x,zcolor='style', ellipse = TRUE)
dev.off()
saveRDS(up,"Up.rds")
saveRDS(down,"Down.rds")
#####
library(ggplot2)
library(clusterProfiler)
library(org.Hs.eg.db)
#up <-readRDS("Up.rds")
#down <-readRDS("Down.rds")
gs <- up
gs = bitr(up, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
ego.bp = enrichGO(gene=gs$ENTREZID, OrgDb = org.Hs.eg.db,ont= "BP",pAdjustMethod = "BH",pvalueCutoff= 0.05,qvalueCutoff= 0.2,readable= TRUE)								
head(ego.bp)
write.csv(ego.bp, file = paste("UP","Go.csv",sep =""))
pdf(paste("UP","Go.pdf",sep =""),width=8,height=6)
print(barplot(ego.bp, showCategory=10,title="GO_biological",drop=T))
dev.off()
gs <- down
gs = bitr(down, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
ego.bp = enrichGO(gene=gs$ENTREZID, OrgDb = org.Hs.eg.db,ont= "BP",pAdjustMethod = "BH",pvalueCutoff= 0.05,qvalueCutoff= 0.2,readable= TRUE)												
head(ego.bp)
write.csv(ego.bp, file =paste("DOWN","Go.csv",sep =""))
pdf(paste("DOWN","Go.pdf",sep =""),width=12,height=6)
print(barplot(ego.bp, showCategory=10,title="GO_biological",drop=T))
dev.off()
##
gs = bitr(up, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
ego.bp = enrichGO(gene=gs$ENTREZID, OrgDb = org.Hs.eg.db,ont= "BP",pAdjustMethod = "BH",pvalueCutoff= 0.05,qvalueCutoff= 0.2,readable= TRUE)
library(clusterProfiler)
library(DOSE)
library(enrichplot)
library(tidytree)
edox2 <- pairwise_termsim(ego.bp)
p1 <- treeplot(edox2)
p2 <- treeplot(edox2, hclust_method = "average")
ggsave(p1,filename =paste("cnetplot_Up","json","_tree1.pdf"), width =12,height =10)
ggsave(p2,filename =paste("cnetplot_Up","json","_tree2.pdf"), width =15,height =10)
#####
gs = bitr(down, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
ego.bp = enrichGO(gene=gs$ENTREZID, OrgDb = org.Hs.eg.db,ont= "BP",pAdjustMethod = "BH",pvalueCutoff= 0.05,qvalueCutoff= 0.2,readable= TRUE)
library(clusterProfiler)
library(DOSE)
library(enrichplot)
edox2 <- pairwise_termsim(ego.bp)
p1 <- treeplot(edox2)
p2 <- treeplot(edox2, hclust_method = "average")
ggsave(p1,filename =paste("cnetplot_Down","json","_tree1.pdf"), width =12,height =10)
ggsave(p2,filename =paste("cnetplot_Down","json","_tree2.pdf"), width =15,height =10)

###
gs = bitr(down, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
ego.bp = enrichGO(gene=gs$ENTREZID, OrgDb = org.Hs.eg.db,ont= "BP",pAdjustMethod = "BH",pvalueCutoff= 0.05,qvalueCutoff= 0.2,readable= TRUE)
y <-c("T cell receptor signaling pathway","antigen receptor-mediated signaling pathway","alpha-beta T cell activation","regulation of T cell activation","immune response-activating cell surface receptor signaling pathway")
cnetp2 <- cnetplot(ego.bp,  
                   showCategory = y,
                   node_label = 'gene',
                   circular = T, 
                   colorEdge = T)
ggsave(cnetp2,filename = paste("cnetplot_","json","_cir2.pdf"), width =15,height =10)

###
bp1 <-read.csv("UPGo.csv")
bp2 <-read.csv("DOWNGo.csv")
bp1$celltype <-"UP"
bp2$celltype <-"DOWN"
bp <-rbind(bp1[c(1,2,3,4,5),],bp2[c(1,2,3,4,5),])
levels <-rev(bp$Description)
pdf("UP_DOWN_GO_function1.pdf",width=10,height=11)
ggplot(bp,aes(x=factor(Description,levels=levels),y=Count,fill=pvalue)) + 
  geom_bar(stat='identity',color='black',width = 0.65) +
  coord_flip() +    
  scale_fill_gradient(low='#FFFFF0', high='darkgoldenrod1')
dev.off()
#####################
####################

sub_OS$newcluster <-paste0(sub_OS@meta.data$seurat_clusters,sub_OS@meta.data$donor_status2)
DefaultAssay(sub_OS) <- "RNA"
sub_OS@active.ident <-factor(sub_OS$newcluster)
#gene
setwd("./Analysis/Tcell/DEG1")
marker_6 <-FindMarkers(sub_OS,ident.1=c("6Re_Me"),ident.2=c("6Chemo"), min.pct = 0.25, logfc.threshold = 0.25)
marker_6 <-marker_6[which(marker_6[,1]<=0.05),]
up_A <-rownames(marker_6[which(marker_6[,2]>=0),])
down_A <-rownames(marker_6[which(marker_6[,2]<=0),])
length(up_A)
length(down_A)
marker_1 <-FindMarkers(sub_OS,ident.1=c("1Re_Me"),ident.2=c("1Chemo"), min.pct = 0.25, logfc.threshold = 0.25)
marker_1 <-marker_1[which(marker_1[,1]<=0.05),]
up_B <-rownames(marker_1[which(marker_1[,2]>=0),])
down_B <-rownames(marker_1[which(marker_1[,2]<=0),])
length(up_B)
length(down_B)
marker_2 <-FindMarkers(sub_OS,ident.1=c("2Re_Me"),ident.2=c("2Chemo"), min.pct = 0.25, logfc.threshold = 0.25)
marker_2 <-marker_2[which(marker_2[,1]<=0.05),]
up_C <-rownames(marker_2[which(marker_2[,2]>=0),])
down_C <-rownames(marker_2[which(marker_2[,2]<=0),])
length(up_C)
length(down_C)
marker_8 <-FindMarkers(sub_OS,ident.1=c("8Re_Me"),ident.2=c("8Chemo"), min.pct = 0.25, logfc.threshold = 0.25)
marker_8 <-marker_8[which(marker_8[,1]<=0.05),]
up_D <-rownames(marker_8[which(marker_8[,2]>=0),])
down_D <-rownames(marker_8[which(marker_8[,2]<=0),])
length(up_D)
length(down_D)
marker_3 <-FindMarkers(sub_OS,ident.1=c("3Re_Me"),ident.2=c("3Chemo"), min.pct = 0.25, logfc.threshold = 0.25)
marker_3 <-marker_3[which(marker_3[,1]<=0.05),]
up_E <-rownames(marker_3[which(marker_3[,2]>=0),])
down_E <-rownames(marker_3[which(marker_3[,2]<=0),])
length(up_E)
length(down_E)

up <-intersect(up_A,intersect(up_B,intersect(up_C,up_E)))
down <-intersect(down_A,intersect(down_B,intersect(down_C,down_E)))
#
##Venn
library(venn)
x <- list(up = up_A, up1 = up_B, up2 = up_C,up4 = up_E)
pdf("up.pdf",width=8,height=6)
venn(x,zcolor='style', ellipse = TRUE)
dev.off()
x <- list(down = down_A, down1 = down_B, down2 = down_C,down4 = down_E)
pdf("down.pdf",width=8,height=6)
venn(x,zcolor='style', ellipse = TRUE)
dev.off()
saveRDS(up,"Up.rds")
saveRDS(down,"Down.rds")
#####
library(ggplot2)
library(clusterProfiler)
library(org.Hs.eg.db)
#up <-readRDS("Up.rds")
#down <-readRDS("Down.rds")
gs <- up
gs = bitr(up, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
ego.bp = enrichGO(gene=gs$ENTREZID, OrgDb = org.Hs.eg.db,ont= "BP",pAdjustMethod = "BH",pvalueCutoff= 0.05,qvalueCutoff= 0.2,readable= TRUE)								
head(ego.bp)
write.csv(ego.bp, file = paste("UP","Go.csv",sep =""))
pdf(paste("UP","Go.pdf",sep =""),width=8,height=6)
print(barplot(ego.bp, showCategory=10,title="GO_biological",drop=T))
dev.off()
gs <- down
gs = bitr(down, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
ego.bp = enrichGO(gene=gs$ENTREZID, OrgDb = org.Hs.eg.db,ont= "BP",pAdjustMethod = "BH",pvalueCutoff= 0.05,qvalueCutoff= 0.2,readable= TRUE)												
head(ego.bp)
write.csv(ego.bp, file =paste("DOWN","Go.csv",sep =""))
pdf(paste("DOWN","Go.pdf",sep =""),width=12,height=6)
print(barplot(ego.bp, showCategory=10,title="GO_biological",drop=T))
dev.off()
##
gs = bitr(up, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
ego.bp = enrichGO(gene=gs$ENTREZID, OrgDb = org.Hs.eg.db,ont= "BP",pAdjustMethod = "BH",pvalueCutoff= 0.05,qvalueCutoff= 0.2,readable= TRUE)
library(clusterProfiler)
library(DOSE)
library(enrichplot)
library(tidytree)
edox2 <- pairwise_termsim(ego.bp)
p1 <- treeplot(edox2)
p2 <- treeplot(edox2, hclust_method = "average")
ggsave(p1,filename =paste("cnetplot_Up","json","_tree1.pdf"), width =12,height =10)
ggsave(p2,filename =paste("cnetplot_Up","json","_tree2.pdf"), width =15,height =10)
#####
gs = bitr(down, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
ego.bp = enrichGO(gene=gs$ENTREZID, OrgDb = org.Hs.eg.db,ont= "BP",pAdjustMethod = "BH",pvalueCutoff= 0.05,qvalueCutoff= 0.2,readable= TRUE)
library(clusterProfiler)
library(DOSE)
library(enrichplot)
edox2 <- pairwise_termsim(ego.bp)
p1 <- treeplot(edox2)
p2 <- treeplot(edox2, hclust_method = "average")
ggsave(p1,filename =paste("cnetplot_Down","json","_tree1.pdf"), width =12,height =10)
ggsave(p2,filename =paste("cnetplot_Down","json","_tree2.pdf"), width =15,height =10)

###
gs = bitr(down, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
ego.bp = enrichGO(gene=gs$ENTREZID, OrgDb = org.Hs.eg.db,ont= "BP",pAdjustMethod = "BH",pvalueCutoff= 0.05,qvalueCutoff= 0.2,readable= TRUE)
y <-c("T cell receptor signaling pathway","antigen receptor-mediated signaling pathway","alpha-beta T cell activation","regulation of T cell activation","immune response-activating cell surface receptor signaling pathway")
cnetp2 <- cnetplot(ego.bp,  
                   showCategory = y,
                   node_label = 'gene',
                   circular = T, 
                   colorEdge = T)
ggsave(cnetp2,filename = paste("cnetplot_","json","_cir2.pdf"), width =15,height =10)

###
bp1 <-read.csv("UPGo.csv")
bp2 <-read.csv("DOWNGo.csv")
bp1$celltype <-"UP"
bp2$celltype <-"DOWN"
bp <-rbind(bp1[c(1,2,3,4,5),],bp2[c(1,2,3,4,5),])
levels <-rev(bp$Description)
pdf("UP_DOWN_GO_function1.pdf",width=10,height=11)
ggplot(bp,aes(x=factor(Description,levels=levels),y=Count,fill=pvalue)) + 
  geom_bar(stat='identity',color='black',width = 0.65) +
  coord_flip() +    
  scale_fill_gradient(low='#FFFFF0', high='darkgoldenrod1')
dev.off()



####survival analysis
data4 <-readRDS("./Bulk-seq/sur_data.rds")
rownames(data4) <- gsub("-", ".", rownames(data4)) 
########
library(ggpubr)
library(survival)
library(survminer)
# add.all = TRUE  

gene <-gsub("-", ".",c("event","OS",down))
gene1 <-intersect(gene,rownames(data4))
dat <-as.data.frame(apply(t(data4[gene1,]),2,as.numeric))
covariates <- colnames(dat)[3:ncol(dat)]
univ_formulas <- sapply(covariates,function(x) as.formula(paste('Surv(OS, event)~', x)))
univ_models <- lapply( univ_formulas, function(x){coxph(x, dat)})
# Extract data
univ_results <- lapply(univ_models,function(x){
x <- summary(x)
p.value<-signif(x$wald["pvalue"], digits=2)
wald.test<-signif(x$wald["test"], digits=2)
beta<-signif(x$coef[1], digits=2)
HR <-signif(x$coef[2], digits=2);#exp(beta)
HR.confint.lower <- signif(x$conf.int[,"lower .95"], 2)
HR.confint.upper <- signif(x$conf.int[,"upper .95"],2)
HR <- paste0(HR, " (",HR.confint.lower, "-", HR.confint.upper, ")")
res<-c(beta, HR, wald.test, p.value)
names(res)<-c("beta", "HR (95% CI for HR)", "wald.test",
"p.value")
return(res)
})
res <- as.data.frame(t(as.data.frame(univ_results, check.names = FALSE)))
res%>% arrange(p.value)
           beta HR (95% CI for HR) wald.test p.value
TRAC       -0.34   0.71 (0.58-0.89)       9.4  0.0022
SH3BGRL3   -0.52    0.6 (0.42-0.84)       8.9  0.0029
ICAM3      -0.35    0.7 (0.54-0.91)       7.3  0.0069
PRKCH       -0.4    0.67 (0.5-0.91)       6.7  0.0094
ARL4C      -0.43   0.65 (0.47-0.91)       6.4   0.012
IL2RG       -0.3   0.74 (0.59-0.94)       6.3   0.012
BRK1       -0.53   0.59 (0.39-0.89)       6.3   0.012
RAB5IF      -0.6   0.55 (0.34-0.88)       6.2   0.013
AKNA       -0.41   0.66 (0.48-0.92)         6   0.015
CORO1A     -0.31   0.73 (0.57-0.94)       5.8   0.016
HCLS1      -0.29   0.75 (0.59-0.95)       5.7   0.017
VPS37B      -0.5    0.61 (0.4-0.92)       5.6   0.018
RHOG       -0.44   0.65 (0.45-0.93)       5.4    0.02
TNFRSF1B   -0.33   0.72 (0.54-0.95)       5.3   0.022
METRNL     -0.35   0.71 (0.52-0.95)       5.2   0.023
GTF2B      -0.54   0.58 (0.36-0.93)       5.2   0.023
LIMD2      -0.42   0.66 (0.45-0.95)       5.1   0.024
ALOX5AP    -0.25   0.78 (0.62-0.97)       5.1   0.024
SRGN       -0.26   0.77 (0.61-0.97)       4.9   0.026
CD3E        -0.2   0.82 (0.69-0.98)       4.9   0.027
CREM       -0.47   0.63 (0.41-0.95)       4.8   0.028
ORAI1      -0.47   0.62 (0.41-0.95)       4.8   0.028
LCK        -0.26   0.77 (0.62-0.98)       4.7    0.03
PTPRC      -0.26   0.77 (0.61-0.98)       4.6   0.031
GSPT1      -0.49   0.61 (0.39-0.96)       4.6   0.032
ARHGAP9    -0.25   0.78 (0.62-0.98)       4.5   0.033
SAMSN1      -0.3   0.74 (0.56-0.98)       4.4   0.035
HLA.A      -0.31   0.73 (0.55-0.99)       4.2   0.041
OAZ1       -0.38   0.68 (0.47-0.99)       4.1   0.043
HNRNPA0    -0.54   0.58 (0.34-0.98)       4.1   0.043
APOBEC3G   -0.33   0.72 (0.53-0.99)       4.1   0.043
TRBC2      -0.17      0.84 (0.71-1)         4   0.046

Gene <-rownames(res%>% arrange(p.value)%>% filter(p.value < 0.05))
Gene <-"TRAC"
for(i in 1:length(Gene)){
gene <-Gene[i]
index <-which(rownames(data4)== gene)
dat <-as.data.frame(apply(t(data4[c(1,2,index),]),2,as.numeric))
dat$Gene=ifelse(dat[,gene] > median(dat[,gene]),'high','low')
fit <-survfit(Surv(OS, event) ~ Gene, data =dat)
p <-ggsurvplot(fit,pval = TRUE, conf.int = TRUE,
       risk.table = FALSE, # Add risk table
       risk.table.col = "strata", # Change risk table color by groups
       linetype = "strata", # Change line type by groups
       surv.median.line = "hv", ggtheme = theme_bw(), 
       main = "Survival curves",palette = c("#E7B800", "#2E9FDF"),xlim = c(0, max(dat[,2])),break.time.by = 365,axes.offset = T,xlab = "Overall Survival (months)",  #x轴的label
       ylab = "Survival Proportion")    
pdf(paste(gene,"_sur.pdf"),width=6,height=6)
print(p)
dev.off()
}



##################
#######
#####Immune Myeloid cells 
####
setwd("./Analysis/Monocell")
sub_OS <-subset(OS,idents= c("0","12","22","11"))
DefaultAssay(sub_OS) <- "integrated"
# Run the standard workflow for visualization and clustering
#sub_OS <- ScaleData(sub_OS, verbose = FALSE)
sub_OS <- RunPCA(sub_OS, npcs = 30, verbose = FALSE)
sub_OS <- RunUMAP(sub_OS, reduction = "pca", dims = 1:30)
sub_OS <- FindNeighbors(sub_OS, reduction = "pca", dims = 1:30)
sub_OS <- FindClusters(sub_OS, resolution = 0.1)
p0 <- DimPlot(object = sub_OS, reduction = "umap",label=TRUE)
p1<- DimPlot(object = sub_OS, reduction = "umap",group.by ="donor_ID",label=TRUE)
p2 <- DimPlot(object = sub_OS, reduction = "umap",group.by ="donor_status",label=TRUE)
pdf("Umap_sub.pdf",width=22,height=20)
CombinePlots(plots  = list(p0, p1, p2))
dev.off()
saveRDS(sub_OS, file="data_rPCA_sub_OS.rds")
DefaultAssay(sub_OS) <- "RNA"
OS.markers <- FindAllMarkers(object = sub_OS, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.table(OS.markers,"markers_r0.1.txt",row.names=TRUE,col.names=TRUE)
##
gene <-c("CD14","FCGR3A","C1QC","CD68","CD163","MRC1","MAF","CD1C","FCER1A","S100A8","S100A9","MKI67","TOP2A") 
library(Seurat)
DefaultAssay(sub_OS) <- "RNA"
sub_OS$seurat_clusters <- factor(sub_OS$seurat_clusters, levels = c("0","1","2","3","4","5"))
plots <- VlnPlot(sub_OS, features =gene,  group.by = "seurat_clusters", pt.size = 0, combine = FALSE)
pdf("Gene.pdf",width=6,height=21)
wrap_plots(plots = plots, ncol = 1)
dev.off() 
##
sub_OS@meta.data$donor_status2 <-factor(sub_OS@meta.data$donor_status2,levels=c("Tumor","Chemo","Re_Me"))
sub_OS@active.ident <-factor(sub_OS$seurat_clusters)
sub_OS_C1 <-subset(sub_OS,idents= c("0"))
DefaultAssay(sub_OS_C1) <- "RNA"
pdf("Gene_86.pdf",width=6,height=3)
print(VlnPlot(sub_OS_C1, features =c("CD86"), split.by = "donor_status2", group.by = "seurat_clusters", pt.size = 0, combine = FALSE))
dev.off() 
pdf("Gene_163.pdf",width=6,height=3)
print(VlnPlot(sub_OS_C1, features =c("CD163"), split.by = "donor_status2", group.by = "seurat_clusters", pt.size = 0, combine = FALSE))
dev.off() 


##
sub_OS$newcluster <-paste0(sub_OS@meta.data$seurat_clusters,sub_OS@meta.data$donor_status2)
sub_OS@active.ident <-factor(sub_OS$newcluster)
gene <-intersect(c("CXXL5","CXCL1","CCL4","CXCL3","CXCL2","CCL20","CCL3L1","CCL3","CCL2","CXCL8","IL1A","IL1B","IL6","TNF"),rownames(sub_OS))
OS_data1 <-sub_OS[["RNA"]]@data[gene,]
OS_data <-t(scale(t(as.matrix(OS_data1))))
type <-unique(sub_OS@meta.data[,"newcluster"])
bb <-NULL
for(i in 1: length(type)){
ind <-which(sub_OS@meta.data[,"newcluster"]==type[i])
aa <-rowMeans(as.matrix(OS_data[,ind]))
bb <-cbind(bb,aa)
}
colnames(bb) <-type
library(reshape2)
library(ggplot2)
bb[bb>4] <-3.8
data_melt<-melt(bb)
data_melt$Var2 <-factor(data_melt$Var2,levels=c("0Tumor","0Chemo","0Re_Me","1Tumor","1Chemo","1Re_Me","2Tumor","2Chemo","2Re_Me","3Tumor","3Chemo","3Re_Me","4Tumor","4Chemo","4Re_Me","5Re_Me"))
p <- ggplot(data = data_melt) + geom_tile(aes(x=Var2,y=Var1, fill = value)) +theme_classic() +  theme(axis.ticks = element_blank(), axis.line = element_blank()) + xlab('row name') + ylab('column name')
p1 <- p+scale_fill_gradient2('z-score',low = 'blue', high = 'red', mid = 'white') 
pdf("heatmap_gene.pdf",width=10,height=15)	
p1
dev.off()
##
setwd("./Analysis/Monocell/function")
marker <-read.table("./Analysis/Monocell/markers_r0.1.txt")
library(ggplot2)
library(clusterProfiler)
library(org.Hs.eg.db)
for(i in 0:5){
ind <-which(marker[,"cluster"]==i)
gene <-marker[ind,"gene"]
up <-gene
gs = bitr(up, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
head(gs)
ego.bp = enrichGO(gene=gs$ENTREZID, OrgDb = org.Hs.eg.db,ont= "BP",pAdjustMethod = "BH",pvalueCutoff= 0.05,qvalueCutoff= 0.2,readable= TRUE)
head(ego.bp)
write.csv(ego.bp, file =paste(i,"_go.csv"))
pdf(paste(i,"_go.pdf"))
print(dotplot(ego.bp, showCategory=10,title="up gene GoTerm"))
dev.off()
kk <- enrichKEGG(gene= gs$ENTREZID, organism= 'hsa',pvalueCutoff = 0.05)
write.csv(kk,file =paste(i,"_kegg.csv"))
pdf(paste(i,"_kegg.pdf"))
print(dotplot(kk, showCategory=10,title="KEGG_Up_biological"))
dev.off()
}
###C0
library(psych)
library(Seurat)
lisub_OS@active.ident <-factor(sub_OS$newcluster)
DefaultAssay(sub_OS) <- "RNA"
cluster.markers <- FindMarkers(sub_OS, ident.1 = c("0Chemo"), ident.2 =c("0Tumor"), min.pct = 0.1,logfc.threshold=0)
filt_data <-cluster.markers
data<-c()
for(i in 1:nrow(filt_data)){
	if(filt_data[i,1]<0.05){
		if(filt_data[i,2]>0.25){
			data<-rbind(data,c(filt_data[i,2],filt_data[i,1],"up"))
		}
		else if(filt_data[i,2]<(-0.25)){
			data<-rbind(data,c(filt_data[i,2],filt_data[i,1],"down"))
		}
		else{
		    data<-rbind(data,c(filt_data[i,2],filt_data[i,1],"no"))
	    }
	}	
	else{
		data<-rbind(data,c(filt_data[i,2],filt_data[i,1],"no"))
	}
}
table(data[,3])
DEG_num=as.numeric(table(data[,3])[1]+table(data[,3])[3])
sub_title=paste("DEGs:",DEG_num,sep="")
colnames(data)<-c("logFC","pvalue",sub_title)
rownames(data)<-rownames(filt_data)
write.table(data,file="Cluster_volcano_plot_format.txt",quote=F,sep="\t")
data <-read.table("Cluster_volcano_plot_format.txt",header=1,sep='\t')
up_num=length(intersect(which(data[,2]<0.05),which(data[,1]>0.25)))
down_num=length(intersect(which(data[,2]<0.05),which(data[,1]<(-0.25))))
legend1=paste("up regulated:",up_num,sep=" ")
legend2=paste("down regulated:",down_num,sep=" ")
index <-which(data[,3]=="no")
index1 <-which(data[,3]=="up")
index2 <-which(data[,3]=="down")
ind <-c(index1,index2,index)
data <-data[ind,]
library("ggplot2")
r03=ggplot(data,aes(x=logFC,y=-log10(pvalue))) +geom_point()+ theme_bw()+theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())
r03=r03 + geom_point(aes(color =DEGs.1323))
volcano = r03 +scale_color_manual(values = c("#009966","#0066CC","#FF0033"),breaks=c("up","down"),labels=c(legend1,legend2))+theme(legend.position = c(.25, .75))
volcano=volcano+geom_hline(yintercept=-log10(0.05),linetype=2,size=1)+geom_vline(xintercept=c((-0.25),0.25),linetype=2,size=1)
library(dplyr)
library(ggrepel)
data$sign <- ifelse(data$logFC >= 1 | data$logFC <= -1,rownames(data),NA)
#data$sign <- ifelse( -log10(data$pvalue) >= 10,rownames(data),NA)
pdf(file="volcano_Chemo_Tumor.pdf",width=7,height=5)	
volcano +geom_text_repel(label=data$sign,colour="black",size=2,box.padding = unit(0.3, "lines"), point.padding = unit(0.4, "lines"), show.legend = F) 
dev.off()
#
library(clusterProfiler)
library(org.Hs.eg.db)
data <-read.table("Cluster_volcano_plot_format.txt",header=1,sep='\t')
up <-rownames(data[intersect(which(data[,2]<0.05),which(data[,1]>0.25)),])
down <-rownames(data[intersect(which(data[,2]<0.05),which(data[,1]<(-0.25))),])
gs <- up
gs = bitr(gs, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
head(gs)
ego.bp = enrichGO(gene =  gs$ENTREZID,OrgDb = org.Hs.eg.db,ont  = "BP",pAdjustMethod = "BH",pvalueCutoff  = 0.05,qvalueCutoff  = 0.2,readable = TRUE)					
head(ego.bp)
# 
nam = paste("Chemo_Tumor_UP","Go.csv",sep ="")
write.csv(ego.bp, file =nam)
nam = paste("Chemo_Tumor_UP","Go.pdf",sep ="")
pdf(nam,width=8,height=6)
print(barplot(ego.bp, showCategory=9,title="GO_biological",drop=T))
dev.off()
#
gs <- down
gs = bitr(gs, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
head(gs)
ego.bp = enrichGO(gene =  gs$ENTREZID,OrgDb = org.Hs.eg.db,ont  = "BP",pAdjustMethod = "BH",pvalueCutoff  = 0.05,qvalueCutoff  = 0.2,readable = TRUE)					
# 
nam = paste("Chemo_Tumor_DOWN","Go.csv",sep ="")
write.csv(ego.bp, file =nam)
nam = paste("Chemo_Tumor_DOWN","Go.pdf",sep ="")
pdf(nam,width=12,height=6)
print(barplot(ego.bp, showCategory=10,title="GO_biological",drop=T))
dev.off()
##
bp1 <-read.csv("Chemo_Tumor_UPGo.csv")
bp2 <-read.csv("Chemo_Tumor_DOWNGo.csv")
bp1$celltype <-"UP"
bp2$celltype <-"DOWN"
bp <-rbind(bp1[c(1,2,3,4,5),],bp2[c(1,2,3,4,142),])
levels <-rev(bp$Description)
pdf("UP_DOWN_GO_function.pdf",width=10,height=11)
ggplot(bp,aes(x=factor(Description,levels=levels),y=Count,fill=pvalue)) + 
  geom_bar(stat='identity',color='black',width = 0.65) +
  coord_flip() +    
  scale_fill_gradient(low='#FFFFF0', high='darkgoldenrod1')
dev.off()


######
bb <-NULL
sample <- unique(sub_OS@meta.data[,"donor_ID"])
for(i in 1:length(sample)){
aa <-table(sub_OS@meta.data[which(sub_OS@meta.data[,"donor_ID"]==sample[i]),"seurat_clusters"])
bb <-cbind(bb,aa)
}
colnames(bb) <-sample
bb <-t(t(bb)/as.vector(table(OS$donor_ID)[sample]))
cc <-NULL
for(i in 1:6){
aa <-bb[i,]/sum(bb[i,])
cc <-rbind(cc,aa)
}
rownames(cc) <-rownames(bb)
library(reshape2)
library(ggplot2)
library(RColorBrewer)
colorr <- brewer.pal(4,"Set2")
df <-melt(cc)
df$Var2 <-factor(df$Var2,levels=colnames(cc))
df$Var1 <-factor(df$Var1,levels=c("0","1","2","3","4","5"))
pdf(file="box_precent_donor.pdf",width=20,height=10)
ggplot(df,aes(x=Var1,y=value,fill=Var2))+scale_fill_manual(values=my36colors[1:17])+geom_bar(stat="identity")+theme(legend.title=element_blank()) 
dev.off()
##
bb <-NULL
sample <- unique(sub_OS@meta.data[,"donor_status2"])
for(i in 1:length(sample)){
aa <-table(sub_OS@meta.data[which(sub_OS@meta.data[,"donor_status2"]==sample[i]),"seurat_clusters"])
bb <-cbind(bb,aa)
}
colnames(bb) <-sample
cc <-NULL
for(i in 1:6){
aa <-bb[i,]/sum(bb[i,])
cc <-rbind(cc,aa)
}
rownames(cc) <-rownames(bb)
library(reshape2)
library(ggplot2)
library(RColorBrewer)
colorr <- brewer.pal(4,"Set2")
df <-melt(cc)
df$Var2 <-factor(df$Var2,levels=colnames(cc))
df$Var1 <-factor(df$Var1,levels=c("0","1","2","3","4","5"))
pdf(file="box_precent_status.pdf",width=20,height=10)
ggplot(df,aes(x=Var1,y=value,fill=Var2))+scale_fill_manual(values=my36colors[1:3])+geom_bar(stat="identity")+theme(legend.title=element_blank()) 
dev.off()

###
library(reshape2)
library(ggplot2)
library(RColorBrewer)
dat <- sub_OS@meta.data[,c("donor_status2","seurat_clusters","donor_ID")]
dat1 <- dat[which(dat$donor_status2 =="Tumor"),]
bb <-NULL
sample <- unique(dat1[,"donor_ID"])
for(i in 1:length(sample)){
aa <-table(dat1[which(dat1[,"donor_ID"]==sample[i]),"seurat_clusters"])
bb <-cbind(bb,aa)}
colnames(bb) <-sample
cc <-t(t(bb)/as.vector(table(OS$donor_ID)[sample]))
dd1 <-melt(cc)
dd1$status <-"Tumor"
dat1 <- dat[which(dat$donor_status2 =="Chemo"),]
bb <-NULL
sample <- unique(dat1[,"donor_ID"])
for(i in 1:length(sample)){
aa <-table(dat1[which(dat1[,"donor_ID"]==sample[i]),"seurat_clusters"])
bb <-cbind(bb,aa)}
colnames(bb) <-sample
cc <-t(t(bb)/as.vector(table(OS$donor_ID)[sample]))
dd2 <-melt(cc)
dd2$status <-"Chemo"
#
dat1 <- dat[which(dat$donor_status2 =="Re_Me"),]
bb <-NULL
sample <- unique(dat1[,"donor_ID"])
for(i in 1:length(sample)){
aa <-table(dat1[which(dat1[,"donor_ID"]==sample[i]),"seurat_clusters"])
bb <-cbind(bb,aa)}
colnames(bb) <-sample
cc <-t(t(bb)/as.vector(table(OS$donor_ID)[sample]))
dd3 <-melt(cc)
dd3$status <-"Re_Me"

df <-rbind(dd1,dd2,dd3)
levels <-c("0","1","2","3","4","5")
df$status = factor(df$status, levels=c("Tumor","Chemo","Re_Me"))
dodge <- position_dodge(width = 0.8)
pdf(file="box_precent_dot_stage.pdf",width=10,height=4)
ggplot(df, aes(x=factor(Var1,levels=levels), y=value,fill = status))  + stat_boxplot(geom="errorbar",width=0.3,position = dodge,size=2)+geom_boxplot(width=0.8,outlier.colour = NA)+scale_fill_manual(values=my36colors[1:15])+
ylim(0,0.6)+theme(legend.title=element_blank())+theme(axis.text.x = element_text(size = 12,angle = 45, hjust = 1, vjust = 1, face = "bold"))+labs(x="Cluster", y = "overlap")
dev.off()
#############t.test/wilcox
bb <-NULL
for(i in c("0","1","2","3","4","5")){
aa <- c(wilcox.test(df[which(df[,1]==i),][which(df[which(df[,1]==i),][,4]=="Tumor"),3],df[which(df[,1]==i),][which(df[which(df[,1]==i),][,4]=="Chemo"),3])$p.value,
wilcox.test(df[which(df[,1]==i),][which(df[which(df[,1]==i),][,4]=="Tumor"),3],df[which(df[,1]==i),][which(df[which(df[,1]==i),][,4]=="Re_Me"),3])$p.value,
wilcox.test(df[which(df[,1]==i),][which(df[which(df[,1]==i),][,4]=="Chemo"),3],df[which(df[,1]==i),][which(df[which(df[,1]==i),][,4]=="Re_Me"),3])$p.value)
bb <-rbind(bb,aa)
}
bb
##
      [,1]       [,2]       [,3]
aa 0.1013986 0.03809524 1.00000000
aa 0.9452214 0.25714286 0.41212121
aa 0.1013986 0.11428571 0.78787879
aa 0.7307692 0.11428571 0.04242424
aa 1.0000000 0.60952381 0.41212121
aa       NaN 0.30743417 0.25683926

#######
UmapPlot<- function(obj) {
  library(ggplot2)
  library(dplyr)
  UMAP = obj@reductions$umap@cell.embeddings %>% as.data.frame() %>% cbind(cell_type = obj@meta.data$seurat_clusters) 
  allcolour=c("#DC143C","#0000FF","#20B2AA","#FFA500","#9370DB","#98FB98","#F08080","#1E90FF","#7CFC00","#FFFF00","#808000","#FF00FF","#FA8072","#7B68EE","#9400D3","#800080","#A0522D","#D2B48C","#D2691E","#87CEEB","#40E0D0","#5F9EA0","#FF1493","#0000CD","#008B8B","#FFE4B5","#8A2BE2","#228B22","#E9967A","#4682B4","#32CD32","#F0E68C","#FFFFE0","#EE82EE","#FF6347","#6A5ACD","#9932CC","#8B008B","#8B4513","#DEB887")
  my36colors <-c('#E5D2DD', '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87','#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#8C549C', '#585658', '#9FA3A8', '#E0D4CA', '#5F3D69', '#C5DEBA', '#58A4C3', '#E4C755', '#F7F398','#AA9A59', '#E63863', '#E39A35', '#C1E6F3', '#6778AE', '#91D0BE', '#B53E2B','#712820', '#DCC1DD', '#CCE0F5',  '#CCC9E6', '#625D9E', '#68A180', '#3A6963','#968175')
  p <- ggplot(UMAP,aes(x= UMAP_1, y = UMAP_2 ,color = cell_type)) +  geom_point(size=0.2, alpha =1 )  +  scale_color_manual(values = my36colors)
  p2 <- p  +theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),  panel.border = element_blank(),  axis.title = element_blank(),  axis.text = element_blank(), axis.ticks = element_blank(),panel.background = element_rect(fill = 'white'), plot.background=element_rect(fill="white"))
  p3 <- p2 + theme(legend.title = element_blank(),legend.key=element_rect(fill='white'), legend.text = element_text(size=20),legend.key.size=unit(1,'cm') ) + guides(color = guide_legend(override.aes = list(size=5))) 
#  p4 <- p3 + geom_segment(aes(x = min(TSNE$tSNE_1) , y = min(TSNE$tSNE_2) ,xend = min(TSNE$tSNE_1) +3, yend = min(TSNE$tSNE_2) ),colour = "black", size=1,arrow = arrow(length = unit(0.3,"cm")))+ geom_segment(aes(x = min(TSNE$tSNE_1), y = min(TSNE$tSNE_2),xend = min(TSNE$tSNE_1) , yend = min(TSNE$tSNE_2) + 3),colour = "black", size=1,arrow = arrow(length = unit(0.3,"cm"))) + annotate("text", x = min(TSNE$tSNE_1) +1.5, y = min(TSNE$tSNE_2) -1, label = "tSNE_1",color="black",size = 3, fontface="bold" ) + annotate("text", x = min(TSNE$tSNE_1) -1, y = min(TSNE$tSNE_2) + 1.5, label = "tSNE_2",color="black",size = 3, fontface="bold" ,angle=90) 
  cell_type_med <- UMAP %>%group_by(cell_type) %>%summarise(UMAP_1 = median(UMAP_1),UMAP_2 = median(UMAP_2))
  library(ggrepel)
  p5 <- p3 + geom_label_repel(aes(label=cell_type), fontface="bold",data = cell_type_med, point.padding=unit(1.5, "lines")) + theme(legend.position = "none")                  
  return(p5)
}
pdf("UMAP_Cluster.pdf",width=4,height=4)
UmapPlot(sub_OS)
dev.off()

UmapPlot<- function(obj) {
  library(ggplot2)
  library(dplyr)
  UMAP = obj@reductions$umap@cell.embeddings %>% as.data.frame() %>% cbind(cell_type = obj@meta.data$donor_ID) 
  allcolour=c("#DC143C","#0000FF","#20B2AA","#FFA500","#9370DB","#98FB98","#F08080","#1E90FF","#7CFC00","#FFFF00","#808000","#FF00FF","#FA8072","#7B68EE","#9400D3","#800080","#A0522D","#D2B48C","#D2691E","#87CEEB","#40E0D0","#5F9EA0","#FF1493","#0000CD","#008B8B","#FFE4B5","#8A2BE2","#228B22","#E9967A","#4682B4","#32CD32","#F0E68C","#FFFFE0","#EE82EE","#FF6347","#6A5ACD","#9932CC","#8B008B","#8B4513","#DEB887")
  my36colors <-c('#E5D2DD', '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87','#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#8C549C', '#585658', '#9FA3A8', '#E0D4CA', '#5F3D69', '#C5DEBA', '#58A4C3', '#E4C755', '#F7F398','#AA9A59', '#E63863', '#E39A35', '#C1E6F3', '#6778AE', '#91D0BE', '#B53E2B','#712820', '#DCC1DD', '#CCE0F5',  '#CCC9E6', '#625D9E', '#68A180', '#3A6963','#968175')
   p <- ggplot(UMAP,aes(x= UMAP_1, y = UMAP_2 ,color = cell_type)) +  geom_point(size=0.2, alpha =1 )  +  scale_color_manual(values = my36colors)
  p2 <- p  +theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),  panel.border = element_blank(),  axis.title = element_blank(),  axis.text = element_blank(), axis.ticks = element_blank(),panel.background = element_rect(fill = 'white'), plot.background=element_rect(fill="white"))
  p3 <- p2 + theme(legend.title = element_blank(),legend.key=element_rect(fill='white'), legend.text = element_text(size=20),legend.key.size=unit(1,'cm') ) + guides(color = guide_legend(override.aes = list(size=5))) 
#  p4 <- p3 + geom_segment(aes(x = min(TSNE$tSNE_1) , y = min(TSNE$tSNE_2) ,xend = min(TSNE$tSNE_1) +3, yend = min(TSNE$tSNE_2) ),colour = "black", size=1,arrow = arrow(length = unit(0.3,"cm")))+ geom_segment(aes(x = min(TSNE$tSNE_1), y = min(TSNE$tSNE_2),xend = min(TSNE$tSNE_1) , yend = min(TSNE$tSNE_2) + 3),colour = "black", size=1,arrow = arrow(length = unit(0.3,"cm"))) + annotate("text", x = min(TSNE$tSNE_1) +1.5, y = min(TSNE$tSNE_2) -1, label = "tSNE_1",color="black",size = 3, fontface="bold" ) + annotate("text", x = min(TSNE$tSNE_1) -1, y = min(TSNE$tSNE_2) + 1.5, label = "tSNE_2",color="black",size = 3, fontface="bold" ,angle=90) 
  cell_type_med <- UMAP %>%group_by(cell_type) %>%summarise(UMAP_1 = median(UMAP_1),UMAP_2 = median(UMAP_2))
  library(ggrepel)
  p5 <- p3 + geom_label_repel(aes(label=cell_type), fontface="bold",data = cell_type_med, point.padding=unit(1.5, "lines")) + theme(legend.position = "none")                  
  return(p5)
}
pdf("UMAP_donor.pdf",width=10,height=10)
UmapPlot(sub_OS)
dev.off()

UmapPlot<- function(obj) {
  library(ggplot2)
  library(dplyr)
  UMAP = obj@reductions$umap@cell.embeddings %>% as.data.frame() %>% cbind(cell_type = obj@meta.data$donor_status1) 
  allcolour=c("#DC143C","#0000FF","#20B2AA","#FFA500","#9370DB","#98FB98","#F08080","#1E90FF","#7CFC00","#FFFF00","#808000","#FF00FF","#FA8072","#7B68EE","#9400D3","#800080","#A0522D","#D2B48C","#D2691E","#87CEEB","#40E0D0","#5F9EA0","#FF1493","#0000CD","#008B8B","#FFE4B5","#8A2BE2","#228B22","#E9967A","#4682B4","#32CD32","#F0E68C","#FFFFE0","#EE82EE","#FF6347","#6A5ACD","#9932CC","#8B008B","#8B4513","#DEB887")
  my36colors <-c('#E5D2DD', '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87','#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#8C549C', '#585658', '#9FA3A8', '#E0D4CA', '#5F3D69', '#C5DEBA', '#58A4C3', '#E4C755', '#F7F398','#AA9A59', '#E63863', '#E39A35', '#C1E6F3', '#6778AE', '#91D0BE', '#B53E2B','#712820', '#DCC1DD', '#CCE0F5',  '#CCC9E6', '#625D9E', '#68A180', '#3A6963','#968175')
  p <- ggplot(UMAP,aes(x= UMAP_1, y = UMAP_2 ,color = cell_type)) +  geom_point(size=0.2, alpha =1 )  +  scale_color_manual(values = my36colors)
  p2 <- p  +theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),  panel.border = element_blank(),  axis.title = element_blank(),  axis.text = element_blank(), axis.ticks = element_blank(),panel.background = element_rect(fill = 'white'), plot.background=element_rect(fill="white"))
  p3 <- p2 + theme(legend.title = element_blank(),legend.key=element_rect(fill='white'), legend.text = element_text(size=20),legend.key.size=unit(1,'cm') ) + guides(color = guide_legend(override.aes = list(size=5))) 
#  p4 <- p3 + geom_segment(aes(x = min(TSNE$tSNE_1) , y = min(TSNE$tSNE_2) ,xend = min(TSNE$tSNE_1) +3, yend = min(TSNE$tSNE_2) ),colour = "black", size=1,arrow = arrow(length = unit(0.3,"cm")))+ geom_segment(aes(x = min(TSNE$tSNE_1), y = min(TSNE$tSNE_2),xend = min(TSNE$tSNE_1) , yend = min(TSNE$tSNE_2) + 3),colour = "black", size=1,arrow = arrow(length = unit(0.3,"cm"))) + annotate("text", x = min(TSNE$tSNE_1) +1.5, y = min(TSNE$tSNE_2) -1, label = "tSNE_1",color="black",size = 3, fontface="bold" ) + annotate("text", x = min(TSNE$tSNE_1) -1, y = min(TSNE$tSNE_2) + 1.5, label = "tSNE_2",color="black",size = 3, fontface="bold" ,angle=90) 
  cell_type_med <- UMAP %>%group_by(cell_type) %>%summarise(UMAP_1 = median(UMAP_1),UMAP_2 = median(UMAP_2))
  library(ggrepel)
  p5 <- p3 + geom_label_repel(aes(label=cell_type), fontface="bold",data = cell_type_med, point.padding=unit(1.5, "lines")) + theme(legend.position = "none")                  
  return(p5)
}
pdf("UMAP_stage.pdf",width=4,height=4)
UmapPlot(sub_OS)
dev.off()


########################
####GSEA
library(msigdbr) 
library(fgsea)
library(dplyr)
library(tibble)
library(Seurat)
library(ggplot2)
library(clusterProfiler)
library(org.Hs.eg.db) 
sub_OS$newcluster <-paste0(sub_OS@meta.data$seurat_clusters,sub_OS@meta.data$donor_status)
DefaultAssay(sub_OS) <- "RNA"
sub_OS@active.ident <-factor(sub_OS$newcluster)
markers <-FindMarkers(sub_OS,ident.1=c("0Chemo"),ident.2=c("0Tumor"), min.pct = 0.25, logfc.threshold = 0)
markers$genes = rownames(markers)
cluster.genes<- markers %>% arrange(desc(avg_log2FC)) %>% dplyr::select(genes,avg_log2FC) 
ranks<- deframe(cluster.genes)
m_df <- msigdbr(species = "Homo sapiens", category = "H")
fgsea_sets = m_df %>% split(x = .$gene_symbol, f = .$gs_name)
length(fgsea_sets)
fgseaRes<- fgsea(fgsea_sets, stats = ranks, nperm = 1000) 
setwd("./Analysis/Monocell/GSVA")
saveRDS(fgseaRes,"C0.rds")
p <-ggplot(fgseaRes %>% as_tibble() %>% arrange(desc(NES)) %>% filter(pval < 0.05), aes(reorder(pathway, NES), NES)) +geom_col(aes(fill= NES)) +coord_flip() +labs(x="KEGG", y="Normalized Enrichment Score",title="KEGG gene sets NES from GSEA") ##输出差异排秩前20的条目
pdf('C0_GSEA-fgsea.pdf')
print(p)
dev.off()
#
markers <-FindMarkers(sub_OS,ident.1=c("1Chemo"),ident.2=c("1Tumor"), min.pct = 0.25, logfc.threshold = 0)
markers$genes = rownames(markers)
cluster.genes<- markers %>% arrange(desc(avg_log2FC)) %>% dplyr::select(genes,avg_log2FC) 
ranks<- deframe(cluster.genes)
m_df <- msigdbr(species = "Homo sapiens", category = "H")
fgsea_sets = m_df %>% split(x = .$gene_symbol, f = .$gs_name)
length(fgsea_sets)
fgseaRes<- fgsea(fgsea_sets, stats = ranks, nperm = 1000) 
saveRDS(fgseaRes,"C1.rds")
p <-ggplot(fgseaRes %>% as_tibble() %>% arrange(desc(NES)) %>% filter(pval < 0.05), aes(reorder(pathway, NES), NES)) +geom_col(aes(fill= NES)) +coord_flip() +labs(x="KEGG", y="Normalized Enrichment Score",title="KEGG gene sets NES from GSEA") ##输出差异排秩前20的条目
pdf('C1_GSEA-fgsea.pdf')
print(p)
dev.off()
#
markers <-FindMarkers(sub_OS,ident.1=c("2Chemo"),ident.2=c("2Tumor"), min.pct = 0.25, logfc.threshold = 0)
markers$genes = rownames(markers)
cluster.genes<- markers %>% arrange(desc(avg_log2FC)) %>% dplyr::select(genes,avg_log2FC) 
ranks<- deframe(cluster.genes)
m_df <- msigdbr(species = "Homo sapiens", category = "H")
fgsea_sets = m_df %>% split(x = .$gene_symbol, f = .$gs_name)
length(fgsea_sets)
fgseaRes<- fgsea(fgsea_sets, stats = ranks, nperm = 1000) 
saveRDS(fgseaRes,"C2.rds")
p <-ggplot(fgseaRes %>% as_tibble() %>% arrange(desc(NES)) %>% filter(pval < 0.05), aes(reorder(pathway, NES), NES)) +geom_col(aes(fill= NES)) +coord_flip() +labs(x="KEGG", y="Normalized Enrichment Score",title="KEGG gene sets NES from GSEA") ##输出差异排秩前20的条目
pdf('C2_GSEA-fgsea.pdf')
print(p)
dev.off()
#
markers <-FindMarkers(sub_OS,ident.1=c("3Chemo"),ident.2=c("3Tumor"), min.pct = 0.25, logfc.threshold = 0)
markers$genes = rownames(markers)
cluster.genes<- markers %>% arrange(desc(avg_log2FC)) %>% dplyr::select(genes,avg_log2FC) 
ranks<- deframe(cluster.genes)
m_df <- msigdbr(species = "Homo sapiens", category = "H")
fgsea_sets = m_df %>% split(x = .$gene_symbol, f = .$gs_name)
length(fgsea_sets)
fgseaRes<- fgsea(fgsea_sets, stats = ranks, nperm = 1000) 
saveRDS(fgseaRes,"C3.rds")
p <-ggplot(fgseaRes %>% as_tibble() %>% arrange(desc(NES)) %>% filter(pval < 0.05), aes(reorder(pathway, NES), NES)) +geom_col(aes(fill= NES)) +coord_flip() +labs(x="KEGG", y="Normalized Enrichment Score",title="KEGG gene sets NES from GSEA") ##输出差异排秩前20的条目
pdf('C3_GSEA-fgsea.pdf')
print(p)
dev.off()
#
markers <-FindMarkers(sub_OS,ident.1=c("4Chemo"),ident.2=c("4Tumor"), min.pct = 0.25, logfc.threshold = 0)
markers$genes = rownames(markers)
cluster.genes<- markers %>% arrange(desc(avg_log2FC)) %>% dplyr::select(genes,avg_log2FC) 
ranks<- deframe(cluster.genes)
m_df <- msigdbr(species = "Homo sapiens", category = "H")
fgsea_sets = m_df %>% split(x = .$gene_symbol, f = .$gs_name)
length(fgsea_sets)
fgseaRes<- fgsea(fgsea_sets, stats = ranks, nperm = 1000) 
saveRDS(fgseaRes,"C4.rds")
p <-ggplot(fgseaRes %>% as_tibble() %>% arrange(desc(NES)) %>% filter(pval < 0.05), aes(reorder(pathway, NES), NES)) +geom_col(aes(fill= NES)) +coord_flip() +labs(x="KEGG", y="Normalized Enrichment Score",title="KEGG gene sets NES from GSEA") ##输出差异排秩前20的条目
pdf('C4_GSEA-fgsea.pdf')
print(p)
dev.off()

#
aa <-readRDS("C0.rds")
aa1 <-readRDS("C1.rds")
aa2 <-readRDS("C2.rds")
aa3 <-readRDS("C3.rds")
aa4 <-readRDS("C4.rds")
bb <-aa %>% as_tibble() %>% arrange(desc(NES)) %>% filter(pval < 0.05)
cc <-cbind(bb$pathway,bb$NES,"C0")
bb1 <-aa1 %>% as_tibble() %>% arrange(desc(NES)) %>% filter(pval < 0.05)
cc1 <-cbind(bb1$pathway,bb1$NES,"C1")
bb2 <-aa2 %>% as_tibble() %>% arrange(desc(NES)) %>% filter(pval < 0.05)
cc2 <-cbind(bb2$pathway,bb2$NES,"C2")
bb3 <-aa3 %>% as_tibble() %>% arrange(desc(NES)) %>% filter(pval < 0.05)
cc3 <-cbind(bb3$pathway,bb3$NES,"C3")
bb4 <-aa4 %>% as_tibble() %>% arrange(desc(NES)) %>% filter(pval < 0.05)
cc4 <-cbind(bb4$pathway,bb4$NES,"C4")

dd <-data.frame(rbind(cc,cc1,cc2,cc3,cc4))
library(pheatmap)
mat <-matrix(0,length(unique(dd$X1)),length(unique(dd$X3)))
rownames(mat) <-unique(dd$X1)
colnames(mat) <-unique(dd$X3)
for(i in 1:nrow(dd)){
mat[dd[i,1],dd[i,3]] <-as.numeric(dd[i,2])
}
min(mat)
max(mat)
bk <- c(seq(-3.3,-0.01,by=0.01),seq(0.01,3.6,by=0.01))
rownames(mat) <-gsub("kegg ","",gsub("_"," ",tolower(rownames(mat))))
pdf("pheatmap_Func_H.pdf",width=15,height=15)
pheatmap(mat[,c("C0","C1","C2","C3","C4")],cluster_rows = T, cluster_cols=F,show_rownames = T,main = "Heatmap",color = c(colorRampPalette(colors = c('#0A7EB5','#F0F8FF'))(33),colorRampPalette(colors = c('#F0F8FF','#EE6A50'))(36)))
dev.off()
##gene
setwd("./Analysis/Monocell/DEG")
marker_0 <-FindMarkers(sub_OS,ident.1=c("0Chemo"),ident.2=c("0Tumor"), min.pct = 0.25, logfc.threshold = 0.25)
marker_0 <-marker_0[which(marker_0[,1]<=0.05),]
up_A <-rownames(marker_0[which(marker_0[,2]>=0),])
down_A <-rownames(marker_0[which(marker_0[,2]<=0),])
length(up_A)
length(down_A)
marker_1 <-FindMarkers(sub_OS,ident.1=c("1Chemo"),ident.2=c("1Tumor"), min.pct = 0.25, logfc.threshold = 0.25)
marker_1 <-marker_1[which(marker_1[,1]<=0.05),]
up_B <-rownames(marker_1[which(marker_1[,2]>=0),])
down_B <-rownames(marker_1[which(marker_1[,2]<=0),])
length(up_B)
length(down_B)
marker_2 <-FindMarkers(sub_OS,ident.1=c("2Chemo"),ident.2=c("2Tumor"), min.pct = 0.25, logfc.threshold = 0.25)
marker_2 <-marker_2[which(marker_2[,1]<=0.05),]
up_C <-rownames(marker_2[which(marker_2[,2]>=0),])
down_C <-rownames(marker_2[which(marker_2[,2]<=0),])
length(up_C)
length(down_C)
marker_3 <-FindMarkers(sub_OS,ident.1=c("3Chemo"),ident.2=c("3Tumor"), min.pct = 0.25, logfc.threshold = 0.25)
marker_3 <-marker_3[which(marker_3[,1]<=0.05),]
up_D <-rownames(marker_3[which(marker_3[,2]>=0),])
down_D <-rownames(marker_3[which(marker_3[,2]<=0),])
length(up_D)
length(down_D)
marker_4 <-FindMarkers(sub_OS,ident.1=c("4Chemo"),ident.2=c("4Tumor"), min.pct = 0.25, logfc.threshold = 0.25)
marker_4 <-marker_4[which(marker_3[,1]<=0.05),]
up_E <-rownames(marker_4[which(marker_3[,2]>=0),])
down_E <-rownames(marker_4[which(marker_3[,2]<=0),])
length(up_E)
length(down_E)

up <-intersect(up_A,intersect(up_B,intersect(up_C,intersect(up_D,up_E))))
down <-intersect(down_A,intersect(down_B,intersect(down_C,intersect(down_D,down_E))))
#
##Venn
library(venn)
x <- list(up = up_A, up1 = up_B, up2 = up_C,up3 = up_D,up4 = up_E)
pdf("up.pdf",width=8,height=6)
venn(x,zcolor='style', ellipse = TRUE)
dev.off()
x <- list(down = down_A, down1 = down_B, down2 = down_C,down3 = down_D,down4 = down_E)
pdf("down.pdf",width=8,height=6)
venn(x,zcolor='style', ellipse = TRUE)
dev.off()
saveRDS(up,"Up.rds")
saveRDS(down,"Down.rds")

#####
library(ggplot2)
library(clusterProfiler)
library(org.Hs.eg.db)
#up <-readRDS("Up.rds")
#down <-readRDS("Down.rds")
gs <- up
gs = bitr(up, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
ego.bp = enrichGO(gene=gs$ENTREZID, OrgDb = org.Hs.eg.db,ont= "BP",pAdjustMethod = "BH",pvalueCutoff= 0.05,qvalueCutoff= 0.2,readable= TRUE)								
head(ego.bp)
write.csv(ego.bp, file = paste("UP","Go.csv",sep =""))
pdf(paste("UP","Go.pdf",sep =""),width=8,height=6)
print(barplot(ego.bp, showCategory=9,title="GO_biological",drop=T))
dev.off()
kk <- enrichKEGG(gene= gs$ENTREZID, organism   = 'mmu',pvalueCutoff = 0.05)
write.csv(kk,file ="UP_kegg.csv")
pdf("UP_kegg.pdf")
print(dotplot(kk, showCategory=10,title="KEGG_Up_biological"))
dev.off()
gs <- down
gs = bitr(down, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
ego.bp = enrichGO(gene=gs$ENTREZID, OrgDb = org.Hs.eg.db,ont= "BP",pAdjustMethod = "BH",pvalueCutoff= 0.05,qvalueCutoff= 0.2,readable= TRUE)												
head(ego.bp)
write.csv(ego.bp, file =paste("DOWN","Go.csv",sep =""))
pdf(paste("DOWN","Go.pdf",sep =""),width=12,height=6)
print(barplot(ego.bp, showCategory=10,title="GO_biological",drop=T))
dev.off()
kk <- enrichKEGG(gene= gs$ENTREZID, organism   = 'mmu',pvalueCutoff = 0.05)
write.csv(kk,file ="Down_kegg.csv")
pdf("Down_kegg.pdf")
print(dotplot(kk, showCategory=10,title="KEGG_biological"))
dev.off()
###
bp1 <-read.csv("UPGo.csv")
bp2 <-read.csv("DOWNGo.csv")
bp1$celltype <-"UP"
bp2$celltype <-"DOWN"
bp <-rbind(bp1[c(1,2,3,4,5),],bp2[c(1,5,7,8,26),])
levels <-rev(bp$Description)
pdf("UP_DOWN_GO_function2.pdf",width=10,height=11)
ggplot(bp,aes(x=factor(Description,levels=levels),y=Count,fill=pvalue)) + 
  geom_bar(stat='identity',color='black',width = 0.65) +
  coord_flip() +    
  scale_fill_gradient(low='#FFFFF0', high='darkgoldenrod1')
dev.off()

##
####生存分析
data4 <-readRDS("./Bulk-seq/sur_data.rds")
rownames(data4) <- gsub("-", ".", rownames(data4)) 
########
library(ggpubr)
library(survival)
library(survminer)
# add.all = TRUE  

gene <-gsub("-", ".",c("event","OS",down))
gene1 <-intersect(gene,rownames(data4))
dat <-as.data.frame(apply(t(data4[gene1,]),2,as.numeric))
covariates <- colnames(dat)[3:ncol(dat)]
univ_formulas <- sapply(covariates,function(x) as.formula(paste('Surv(OS, event)~', x)))
univ_models <- lapply( univ_formulas, function(x){coxph(x, dat)})
# Extract data
univ_results <- lapply(univ_models,function(x){
x <- summary(x)
p.value<-signif(x$wald["pvalue"], digits=2)
wald.test<-signif(x$wald["test"], digits=2)
beta<-signif(x$coef[1], digits=2)
HR <-signif(x$coef[2], digits=2);#exp(beta)
HR.confint.lower <- signif(x$conf.int[,"lower .95"], 2)
HR.confint.upper <- signif(x$conf.int[,"upper .95"],2)
HR <- paste0(HR, " (",HR.confint.lower, "-", HR.confint.upper, ")")
res<-c(beta, HR, wald.test, p.value)
names(res)<-c("beta", "HR (95% CI for HR)", "wald.test",
"p.value")
return(res)
})
res <- as.data.frame(t(as.data.frame(univ_results, check.names = FALSE)))
res%>% arrange(p.value) 
VASP        -0.7    0.5 (0.32-0.77)       9.6  0.0019
MSN        -0.64   0.53 (0.35-0.79)       9.4  0.0021
TNFAIP8    -0.63    0.53 (0.36-0.8)       9.4  0.0022
ETV6       -0.61   0.55 (0.37-0.81)         9  0.0027
GNAI2      -0.58   0.56 (0.37-0.84)       7.9  0.0049
SLC43A2    -0.43   0.65 (0.48-0.88)       7.7  0.0056
RBM47      -0.38    0.69 (0.52-0.9)       7.3  0.0071
SH2B3      -0.51    0.6 (0.41-0.87)       7.1  0.0076
RAB5IF      -0.6   0.55 (0.34-0.88)       6.2   0.013
TGIF1      -0.52    0.59 (0.39-0.9)         6   0.014
ANXA5      -0.53    0.59 (0.38-0.9)         6   0.014
FTH1       -0.43   0.65 (0.46-0.93)       5.6   0.018
TNFRSF1B   -0.33   0.72 (0.54-0.95)       5.3   0.022
METRNL     -0.35   0.71 (0.52-0.95)       5.2   0.023
SRGN       -0.26   0.77 (0.61-0.97)       4.9   0.026
CREM       -0.47   0.63 (0.41-0.95)       4.8   0.028
ABL2       -0.46   0.63 (0.42-0.96)       4.7   0.029
RNF149     -0.45   0.64 (0.42-0.97)       4.3   0.037
PLEK       -0.25   0.78 (0.62-0.99)       4.3   0.037
LMNA        -0.4   0.67 (0.46-0.98)       4.3   0.039
FAM107B    -0.39   0.68 (0.46-0.99)       4.1   0.043
WNK1       -0.41   0.66 (0.44-0.99)         4   0.045
ZNF267      -0.4      0.67 (0.45-1)       3.9   0.048
CXCL16     -0.33      0.72 (0.52-1)       3.8    0.05

gene <- rownames(res%>% arrange(p.value)%>% filter(p.value < 0.05))
for(i in 1:length(gene)) {
index <-which(rownames(data4)== gene[i])
dat <-as.data.frame(apply(t(data4[c(1,2,index),]),2,as.numeric))
dat$Gene=ifelse(dat[,3] > median(dat[,3]),'high','low')
fit <-survfit(Surv(OS, event) ~ Gene, data =dat)
p <-ggsurvplot(fit,pval = TRUE, conf.int = TRUE,
       risk.table = FALSE, # Add risk table
       risk.table.col = "strata", # Change risk table color by groups
       linetype = "strata", # Change line type by groups
       surv.median.line = "hv", ggtheme = theme_bw(), 
       main = "Survival curves",palette = c("#E7B800", "#2E9FDF"),xlim = c(0, max(dat[,2])),break.time.by = 365,axes.offset = T,xlab = "Overall Survival (months)",  #x轴的label
       ylab = "Survival Proportion")    
pdf(paste(gene[i],"_sur.pdf"),width=6,height=6)
print(p)
dev.off()
}

####
#####
####









####################
###############
########### communication 
#######cellchat
library(CellChat)
library(ggplot2)
library(ggalluvial)
library(Seurat)
setwd("./Analysis/Cellchat")
sub_OS_TIL <-readRDS("./Analysis/Tcell/data_rPCA_sub_OS.rds")
current.cluster.ids <- c("6","1","2","8","3","4","5","9","7","0")
new.cluster.ids <- c("Pro-T","CD8+T","NK","CD4+T","Treg","B","B","Mast","pDC","Other")
sub_OS_TIL@meta.data$group <- plyr::mapvalues(x = sub_OS_TIL@meta.data[,"seurat_clusters"], from = current.cluster.ids, to = new.cluster.ids)

sub_OS_Mono <-readRDS("./Analysis/Monocell/data_rPCA_sub_OS.rds")
current.cluster.ids <- c(0,1,2,3,4,5)
new.cluster.ids <- c("Mac","Mono","Mono","Pro-MAC/Mono","DC","Neutrophil")
sub_OS_Mono@meta.data$group <- plyr::mapvalues(x = sub_OS_Mono@meta.data[,"seurat_clusters"], from = current.cluster.ids, to = new.cluster.ids)

OS@active.ident <-factor(OS$celltype)
sub_OS_Osteoblastic <-subset(OS,ident=c("Osteoblasts"))
sub_OS_Osteoblastic@meta.data$group <-c("Osteoblasts")
sub_OS_Chondroblastic <-subset(OS,ident=c("Chondroblastic"))
sub_OS_Chondroblastic@meta.data$group <-c("Chondroblastic")
sub_OS_Osteoclast <-subset(OS,ident=c("Osteoclast"))
sub_OS_Osteoclast@meta.data$group <-c("Osteoclast")
sub_OS_MSC <-subset(OS,ident=c("MSC"))
sub_OS_MSC@meta.data$group <-c("MSC")
sub_OS_CAF <-subset(OS,ident=c("Fibroblast"))
sub_OS_CAF@meta.data$group <-c("CAF")
sub_OS_Pericyte <-subset(OS,ident=c("Pericyte"))
sub_OS_Pericyte@meta.data$group <-c("Pericyte")

sub_OS <-merge(sub_OS_Osteoblastic,y=c(sub_OS_Chondroblastic,sub_OS_Osteoclast,sub_OS_MSC,sub_OS_CAF,sub_OS_Pericyte,sub_OS_TIL,sub_OS_Mono))
sub_OS@active.ident <-factor(sub_OS$group)
sub_OS_M <-subset(sub_OS,ident=c("CD8+T","Treg","NK","B","Mac","Mono","DC","Osteoblasts","Chondroblastic","Osteoclast","MSC","CAF","Pericyte"))

saveRDS(sub_OS_M,"sub_OS_M.rds")

#sub_OS_M <-readRDS("./Analysis/Cellchat/sub_OS_M.rds")
current.cluster.ids <- c("B","Treg","CD8+T","NK","Mac","Mono","DC","Osteoblasts","Osteoclast","Chondroblastic","MSC","CAF","Pericyte")
new.cluster.ids <- c("AA_B","AB_Treg","AC_CD8+T","AD_NK","AE_Mac","AF_Mono","AG_DC","AH_Osteoblastic","AI_Osteoclast","AJ_Chondroblastic","AK_MSC","AM_CAF","AN_Pericyte")
sub_OS_M@meta.data$group1 <- plyr::mapvalues(x = sub_OS_M@meta.data[,"group"], from = current.cluster.ids, to = new.cluster.ids)
table(sub_OS_M@meta.data$group1,sub_OS_M@meta.data$donor_status)    
                 Metastasis Chemo rechem  umor
  AA_B                      73      84     10  1137
  AB_Treg                  246     169     39   479
  AC_CD8+T                  36     722     38  1794
  AD_NK                    268     528     55  1124
  AE_Mac                    95    4206    231  8506
  AF_Mono                 1203    9975   1192  5174
  AG_DC                    157    1177     56   498
  AH_Osteoblastic        10614   16628  10259  8887
  AI_Osteoclast            423    7170      0  1791
  AJ_Chondroblastic        110    1939   4251   398
  AK_MSC                  1244    2551     50   426
  AM_CAF                   259   10699     90   752
  AN_Pericyte              459    2893    802   835
  

sub_OS_M@active.ident <-factor(sub_OS_M$donor_status2)
sub_OS1 <-subset(sub_OS_M,ident=c("Tumor"))
sub_OS2 <-subset(sub_OS_M,ident=c("Chemo"))
sub_OS3 <-subset(sub_OS_M,ident=c("Re_Me"))
Tumor <- sub_OS1
Chem <- sub_OS2
rechem <-sub_OS3

cco.Young <- createCellChat(Tumor$RNA@data) 
identity = data.frame(group =as.vector(Tumor$group1), row.names = Cells(Tumor)) 
unique(identity$group) # check the cell labels
cco.Young <- addMeta(cco.Young, meta = identity, meta.name = "labels") 
cco.Young <- setIdent(cco.Young, ident.use = "labels") # set "labels" as default cell identity
levels(cco.Young@idents) # show factor levels of the cell labels
cco.Old <- createCellChat(Chem$RNA@data) 
identity = data.frame(group =as.vector(Chem$group1), row.names = Cells(Chem)) 
unique(identity$group) # check the cell labels
cco.Old <- addMeta(cco.Old, meta = identity, meta.name = "labels") 
cco.Old <- setIdent(cco.Old, ident.use = "labels") # set "labels" as default cell identity
levels(cco.Old@idents) # show factor levels of the cell labels
cco.rechem <- createCellChat(rechem$RNA@data) 
identity = data.frame(group =as.vector(rechem$group1), row.names = Cells(rechem)) 
unique(identity$group) # check the cell labels
cco.rechem <- addMeta(cco.rechem, meta = identity, meta.name = "labels") 
cco.rechem <- setIdent(cco.rechem, ident.use = "labels") # set "labels" as default cell identity
levels(cco.rechem@idents) # show factor levels of the cell labels
save(cco.Young, file = "cco_Tumor.rda")
save(cco.Old, file = "cco_Chem.rda")
save(cco.rechem, file = "cco_Rechem.rda")


setwd("./Analysis/Cellchat/Tumor")
load("./Analysis/Cellchat/cco_Tumor.rda")
cellchat <- cco.Young 
cellchat@DB <- CellChatDB.human # use CellChatDB.human if running on human data
##future::plan("multisession", workers = 1) # do parallel
cellchat <- subsetData(cellchat) 
cellchat <- identifyOverExpressedGenes(cellchat) 
cellchat <- identifyOverExpressedInteractions(cellchat) 
cellchat <- projectData(cellchat, PPI.human)
cellchat <- computeCommunProb(cellchat, raw.use = TRUE, population.size = TRUE) 
cellchat <- filterCommunication(cellchat, min.cells = 10) 
cellchat <- computeCommunProbPathway(cellchat) 
cellchat <- aggregateNet(cellchat) 
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP") 
cellchat <- computeNetSimilarity(cellchat, type = "functional") 
#library(uwot)
#cellchat <- netEmbedding(cellchat, umap.method = 'uwot',type = "functional")
cellchat <- netEmbedding(cellchat,umap.method = 'uwot', type = "functional") 
cellchat <- netClustering1(cellchat, type = "functional") 
cellchat <- computeNetSimilarity(cellchat, type = "structural") 
cellchat <- netEmbedding(cellchat, umap.method = 'uwot', type = "structural") 
cellchat <- netClustering1(cellchat, type = "structural") 
cco.Young <- cellchat 
saveRDS(cco.Young, "cco_Tumor1.rds")
#
setwd("./Analysis/Cellchat/Chemo")
load("./Analysis/Cellchat/cco_Chem.rda")
cellchat <- cco.Old 
cellchat@DB <- CellChatDB.human 
cellchat <- subsetData(cellchat) 
cellchat <- identifyOverExpressedGenes(cellchat) 
cellchat <- identifyOverExpressedInteractions(cellchat) 
cellchat <- projectData(cellchat, PPI.human) 
cellchat <- computeCommunProb(cellchat, raw.use = TRUE, population.size = TRUE) 
cellchat <- filterCommunication(cellchat, min.cells =10) 
cellchat <- computeCommunProbPathway(cellchat) 
cellchat <- aggregateNet(cellchat) 
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP") 
cellchat <- computeNetSimilarity(cellchat, type = "functional") 
cellchat <- netEmbedding(cellchat, umap.method = 'uwot', type = "functional") 
cellchat <- netClustering1(cellchat, type = "functional") 
cellchat <- computeNetSimilarity(cellchat, type = "structural") 
cellchat <- netEmbedding(cellchat, umap.method = 'uwot', type = "structural") 
cellchat <- netClustering1(cellchat, type = "structural") 
cco.Old <- cellchat 
saveRDS(cco.Old, "cco_Chem1.rds")
#
setwd("./Analysis/Cellchat/Re_Me")
load("./Analysis/Cellchat/cco_Rechem.rda")
cellchat <- cco.rechem 
cellchat@DB <- CellChatDB.human 
cellchat <- subsetData(cellchat) 
cellchat <- identifyOverExpressedGenes(cellchat) 
cellchat <- identifyOverExpressedInteractions(cellchat) 
cellchat <- projectData(cellchat, PPI.human) 
cellchat <- computeCommunProb(cellchat, raw.use = TRUE, population.size = TRUE) 
cellchat <- filterCommunication(cellchat, min.cells =10) 
cellchat <- computeCommunProbPathway(cellchat) 
cellchat <- aggregateNet(cellchat) 
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP") 
cellchat <- computeNetSimilarity(cellchat, type = "functional") 
cellchat <- netEmbedding(cellchat, umap.method = 'uwot', type = "functional") 
cellchat <- netClustering1(cellchat, type = "functional") 
cellchat <- computeNetSimilarity(cellchat, type = "structural") 
cellchat <- netEmbedding(cellchat, umap.method = 'uwot', type = "structural") 
cellchat <- netClustering1(cellchat, type = "structural") 
cco.rechem <- cellchat 
saveRDS(cco.rechem, "cco_Rechem1.rds")
###

setwd("./Analysis/Cellchat/Comp")
cco.Young <-readRDS("./Analysis/Cellchat/Tumor/cco_Tumor1.rds")
cco.Old <-readRDS("./Analysis/Cellchat/Chemo/cco_Chem1.rds")
cco.Old1 <-readRDS("./Analysis/Cellchat/Re_Me/cco_rechem1.rds")
object.list <- list(Tumor = cco.Young,Chem = cco.Old,Rechem = cco.Old1) 
cellchat <- mergeCellChat(object.list, add.names = names(object.list),cell.prefix = TRUE)
gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2,3), measure = "count") 
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2,3), measure = "weight") 
p <- gg1 + gg2 
ggsave("Overview_number_strength.pdf", p, width = 6, height = 4)
#
pdf(file = "Counts_Compare_net.pdf",width=15,height=8)
par(mfrow = c(1,2)) 
weight.max <- getMaxWeight(object.list, attribute = c("idents","count")) 
i=1
netVisual_circle(object.list[[i]]@net$count, weight.scale = T, label.edge= F, edge.weight.max = weight.max[2], edge.width.max = 12, 
title.name = paste0("Number of interactions - ", names(object.list)[i]))
i=2
netVisual_circle(object.list[[i]]@net$count, weight.scale = T, label.edge= F, edge.weight.max = weight.max[2], edge.width.max = 12, 
title.name = paste0("Number of interactions - ", names(object.list)[i]))
i=3
netVisual_circle(object.list[[i]]@net$count, weight.scale = T, label.edge= F, edge.weight.max = weight.max[2], edge.width.max = 12, 
title.name = paste0("Number of interactions - ", names(object.list)[i]))
dev.off()
##
pathway.union <- unique(c(object.list[[1]]@netP$pathways, object.list[[2]]@netP$pathways, object.list[[3]]@netP$pathways))
ht1 = netAnalysis_signalingRole_heatmap(object.list[[1]], pattern = "all", signaling = pathway.union, title = names(object.list)[1], width = 8, height = 16, color.heatmap = "OrRd") 
ht2 = netAnalysis_signalingRole_heatmap(object.list[[2]], pattern = "all", signaling = pathway.union, title = names(object.list)[2], width = 8, height = 16, color.heatmap = "OrRd") 
ht3 = netAnalysis_signalingRole_heatmap(object.list[[3]], pattern = "all", signaling = pathway.union, title = names(object.list)[3], width = 8, height = 16, color.heatmap = "OrRd") 
pdf(file = "Compare_signal_pattern_all2.pdf", width = 12, height = 10)
ht1 + ht2 +ht3
dev.off()
###

levels(cellchat@idents$joint) 
p <- netVisual_bubble(cellchat, sources.use = c(1,2,4,5,6,7), targets.use = c(3), comparison = c(1, 2,3), angle.x = 45) 
ggsave("Compare_LR_bubble.pdf", p, width =7, height =14)

p <- netVisual_bubble(cellchat, sources.use = c(5), targets.use = c(8,9,10,11,12,13), comparison = c(1, 2,3), angle.x = 45) 
ggsave("Compare_LR_bubble2.pdf", p, width = 7, height =14)

p <- netVisual_bubble(cellchat, sources.use = c(8), targets.use = c(1,2,3,4,5,6,7), comparison = c(1, 2,3), angle.x = 45) 
ggsave("Compare_LR_bubble12.pdf", p, width = 10, height =14)

###
pathways.show <- c("MHC-I")
weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show)
par(mfrow = c(1,3), xpd = TRUE)
for(i in 1:length(object.list)){
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1],
                      edge.width.max = 10, signaling.name = paste(pathways.show, names(object.list)[i]))
}
######
setwd("./Analysis/Cellchat/Tumor")
pathways.show <- cco.Young@netP$pathways
for (i in 1:length(pathways.show)) 
{
pdf(file = paste0(pathways.show[i], "_gene.pdf"))
print(plotGeneExpression(cco.Young , signaling =pathways.show[i]))
dev.off()
}
setwd("./pathway")
cellchat <-cco.Young
pathways.show <- cco.Young@netP$pathways
groupSize <- as.numeric(table(cellchat@idents)) # number of cells in each cell group
for (i in 1:length(pathways.show)) {
  gg <- netAnalysis_contribution(cellchat, signaling = pathways.show[i])
  ggsave(filename=paste0(pathways.show[i], "_L-R_contribution.pdf"), plot=gg, width = 3, height = 2, units = 'in', dpi = 300)
  # Visualize signaling roles of cell groups
  pdf(file = paste0(pathways.show[i], "_signalRole.pdf"))
  netAnalysis_signalingRole_network(cellchat, signaling = pathways.show[i], width = 8, height = 2.5, font.size = 10)
  dev.off()
}
netVisual_aggregate(cellchat, signaling = "MHC-II", layout = "circle")
for (i in 1:length(pathways.show)) 
{
netVisual(cellchat, signaling = pathways.show[i], layout = "circle", vertex.weight = groupSize,pt.title=20,vertex.label.cex = 1.7)
}

##
setwd("./Analysis/Cellchat/Chemo/pathway")
pathways.show <- cco.Old@netP$pathways
for (i in 1:length(pathways.show)) 
{
netVisual(cco.Old, signaling = pathways.show[i], layout = "circle", vertex.weight = groupSize,pt.title=20,vertex.label.cex = 1.7)
}
setwd("./Analysis/Cellchat/Tumor")
pathways.show <- cco.Old1@netP$pathways
for (i in 1:length(pathways.show)) 
{
netVisual(cco.Old1, signaling = pathways.show[i], layout = "circle", vertex.weight = groupSize,pt.title=20,vertex.label.cex = 1.7)
}


#netVisual_bubble(cellchat, signaling = c("MIF","ANNEXIN") )
pdf(file = "Gene_Compare_net.pdf",width=15,height=8)
par(mfrow = c(1,2), xpd=TRUE)
netVisual_chord_gene(object.list[[1]], sources.use = 5, targets.use = c(8,9,10), slot.name = 'net', lab.cex = 0.8, small.gap = 3.5, title.name = paste0("Up-regulated signaling in ", names(object.list)[2]))
#> Note: The first link end is drawn out of sector 'MIF'.
netVisual_chord_gene(object.list[[2]], sources.use = 5, targets.use = c(8,9,10), slot.name = 'net', lab.cex = 0.8, small.gap = 3.5, title.name = paste0("Down-regulated signaling in ", names(object.list)[2]))
netVisual_chord_gene(object.list[[3]], sources.use = 5, targets.use = c(8,9,10), slot.name = 'net', lab.cex = 0.8, small.gap = 3.5, title.name = paste0("Down-regulated signaling in ", names(object.list)[2]))
dev.off()

##
setwd("./Analysis/Cellchat/Comp")
object.list <- list(Tumor = cco.Young,Chem = cco.Old) 
cellchat <- mergeCellChat(object.list, add.names = names(object.list),cell.prefix = TRUE)
gg1 <- rankNet(cellchat, mode = "comparison", stacked = T, do.stat = TRUE) 
gg2 <- rankNet(cellchat, mode = "comparison", stacked = F, do.stat = TRUE) 
p <- gg1 + gg2 
ggsave("Compare_pathway_strengh.pdf", p, width = 10, height = 8)
##
object.list <- list(Chem = cco.Old,Rechem = cco.Old1) 
cellchat <- mergeCellChat(object.list, add.names = names(object.list),cell.prefix = TRUE)
gg1 <- rankNet(cellchat, mode = "comparison", stacked = T, do.stat = TRUE) 
gg2 <- rankNet(cellchat, mode = "comparison", stacked = F, do.stat = TRUE) 
p <- gg1 + gg2 
ggsave("Compare_pathway_strengh1.pdf", p, width = 10, height = 8)

#####
current.cluster.ids <- c("Tumor","Chemo","Re_Me")
new.cluster.ids <- c("AA_Tumor","AB_Chemo","AC_Re_Me")
sub_OS_M@meta.data$donor_status4 <- plyr::mapvalues(x = sub_OS_M@meta.data[,"donor_status2"], from = current.cluster.ids, to = new.cluster.ids)
#sub_OS_M@active.ident <-factor(sub_OS_M$donor_status4)
#sub_OS_M1 <-subset(sub_OS_M,ident=c("AA_Tumor","AB_Chemo"))
plots <- VlnPlot(sub_OS_M, features = c("HLA-A", "HLA-B","HLA-C", "HLA-E", "HLA-F","HLA-DPA1","HLA-DPA1","HLA-DPB1","HLA-DQA1","HLA-DMA","HLA-DMB","HLA-DQA2","HLA-DOA","HLA-DQB1"), split.by = "donor_status4", group.by = "group1", pt.size = 0, combine = FALSE)
pdf("Gene2.pdf",width=10,height=46)
wrap_plots(plots = plots, ncol = 1)
dev.off()  


#########
library(survival)
library(survminer)
sur <- read.csv("./Bulk-seq/TARGET-OS.clinical.tsv",sep="\t")
sur[,1]= gsub("-", ".", sur[,1])
sur <-sur[grep("01A",sur[,1]),]
data <- read.csv("./Bulk-seq/TARGET-OS.htseq_counts.tsv",row.names=1,sep="\t")
index <-NULL
for(i in 1:88) {
index1 <- which(sur[,1] ==colnames(data)[i])
index <-c(index,index1)
}
sur1 <-sur[index,]
sur1 <-sur1[,c(1,6,7,9,10,14,26,30)]
###0 liver 1 Dead 
sur2 <-sur1[grep("Alive|Dead",sur1[,4]),]
sur3 <-sur2[which(sur2[,5]!=0),]
event <- c("Alive","Dead")
event1 <- c("0","1")
sur3$event <- plyr::mapvalues(x = sur3[,"Vital.Status"], from = event, to = event1)
#sur3$OS <- sur3$Overall.Survival.Time.in.Days/30
sur3$OS <- sur3$Overall.Survival.Time.in.Days
rownames(sur3) <-sur3[,1]
#####

data1 <-data[,sur3[,1]][c(1:60483),]
data1$ensemble <-rownames(data1)
data1$ensemble <- substr(data1$ensemble,1,15)
library(clusterProfiler)
library(org.Hs.eg.db)
# Ensembl_ID
Ensembl_ID <- data1$ensemble
gene_symbol <- bitr(Ensembl_ID, fromType="ENSEMBL", toType=c("SYMBOL", "ENTREZID"), OrgDb="org.Hs.eg.db")
head(gene_symbol)
data2=data.frame(gene_symbol,data1[match(gene_symbol$ENSEMBL,data1$ensemble),]
data2 <-data2[,c(-1,-3,-88)]
data3 <- aggregate(data2[,c(2:85)], by=list(type=data2$SYMBOL),mean)
rownames(data3)=data3$type
data3 <-data3[,-1]
data4 <-rbind(t(sur3[,c(9,10)]),data3)
setwd("./sur")
saveRDS(data4,"sur_data.rds")
data4 <-rbind(t(sur3),data3)
setwd("./sur")
saveRDS(data4,"sur_data1.rds")
########
data4 <-readRDS("./sur/sur_data.rds")  ##FPKM
library(ggpubr)
library(survival)
library(survminer)
# add.all = TRUE  
down <- readRDS("./Analysis/Tcell/DEG/Down.rds")

gene <-gsub("-", ".",c("event","OS",down))
gene1 <-intersect(gene,rownames(data4))
dat <-as.data.frame(apply(t(data4[gene1,]),2,as.numeric))
covariates <- colnames(dat)[3:ncol(dat)]
univ_formulas <- sapply(covariates,function(x) as.formula(paste('Surv(OS, event)~', x)))
univ_models <- lapply( univ_formulas, function(x){coxph(x, dat)})
# Extract data
univ_results <- lapply(univ_models,function(x){
x <- summary(x)
p.value<-signif(x$wald["pvalue"], digits=2)
wald.test<-signif(x$wald["test"], digits=2)
beta<-signif(x$coef[1], digits=2)
HR <-signif(x$coef[2], digits=2);#exp(beta)
HR.confint.lower <- signif(x$conf.int[,"lower .95"], 2)
HR.confint.upper <- signif(x$conf.int[,"upper .95"],2)
HR <- paste0(HR, " (",HR.confint.lower, "-", HR.confint.upper, ")")
res<-c(beta, HR, wald.test, p.value)
names(res)<-c("beta", "HR (95% CI for HR)", "wald.test",
"p.value")
return(res)
})
res <- as.data.frame(t(as.data.frame(univ_results, check.names = FALSE)))
res1 <- res%>% arrange(p.value)
write.csv(res1,"sur_T_down1.csv",row.names=T,col.names=T,sep="\t")
##
setwd("./Analysis/Tcell/DEG")
Gene <-  rownames(res1%>% filter(p.value < 0.05))
for(i in 1:length(Gene)){
gene <-Gene[i]
index <-which(rownames(data4)== gene)
dat <-as.data.frame(apply(t(data4[c(1,2,index),]),2,as.numeric))
dat$Gene=ifelse(dat[,gene] > median(dat[,gene]),'high','low')
fit <-survfit(Surv(OS, event) ~ Gene, data =dat)
p <-ggsurvplot(fit,pval = TRUE, conf.int = TRUE,
       risk.table = FALSE, # Add risk table
       risk.table.col = "strata", # Change risk table color by groups
       linetype = "strata", # Change line type by groups
       surv.median.line = "hv", ggtheme = theme_bw(), 
       main = "Survival curves",palette = c("#E7B800", "#2E9FDF"),xlim = c(0, max(dat[,2])),break.time.by = 365,axes.offset = T,xlab = "Overall Survival (months)",  #x轴的label
       ylab = "Survival Proportion")    
pdf(paste(gene,"_sur.pdf"),width=6,height=6)
print(p)
dev.off()
}

##
down <- readRDS("./Analysis/Monocell/DEG/Down.rds")

gene <-gsub("-", ".",c("event","OS",down))
gene1 <-intersect(gene,rownames(data4))
dat <-as.data.frame(apply(t(data4[gene1,]),2,as.numeric))
covariates <- colnames(dat)[3:ncol(dat)]
univ_formulas <- sapply(covariates,function(x) as.formula(paste('Surv(OS, event)~', x)))
univ_models <- lapply( univ_formulas, function(x){coxph(x, dat)})
# Extract data
univ_results <- lapply(univ_models,function(x){
x <- summary(x)
p.value<-signif(x$wald["pvalue"], digits=2)
wald.test<-signif(x$wald["test"], digits=2)
beta<-signif(x$coef[1], digits=2)
HR <-signif(x$coef[2], digits=2);#exp(beta)
HR.confint.lower <- signif(x$conf.int[,"lower .95"], 2)
HR.confint.upper <- signif(x$conf.int[,"upper .95"],2)
HR <- paste0(HR, " (",HR.confint.lower, "-", HR.confint.upper, ")")
res<-c(beta, HR, wald.test, p.value)
names(res)<-c("beta", "HR (95% CI for HR)", "wald.test",
"p.value")
return(res)
})
res <- as.data.frame(t(as.data.frame(univ_results, check.names = FALSE)))
res1 <- res%>% arrange(p.value)
write.csv(res1,"sur_mono_down1.csv",row.names=T,col.names=T,sep="\t")

#####
###
library(estimate)
exp <- as.matrix(read.table("./Bulk-seq/OS.htseq_count1.txt",header = T,row.names = 1,check.names = F))
dat=log2(edgeR::cpm(exp)+1)
estimate <- function(dat,pro){
  input.f=paste0(pro,'_estimate_input.txt')
  output.f=paste0(pro,'_estimate_gene.gct')
  output.ds=paste0(pro,'_estimate_score.gct')
  write.table(dat,file = input.f,sep = '\t',quote = F)
  library(estimate)
  filterCommonGenes(input.f=input.f,
                    output.f=output.f ,
                    id="GeneSymbol")
  estimateScore(input.ds = output.f,
                output.ds=output.ds,
                platform="illumina")   
  scores=read.table(output.ds,skip = 2,header = T)
  rownames(scores)=scores[,1]
  scores=t(scores[,3:ncol(scores)])
  return(scores)
}
pro='OS'
scores=estimate(dat,pro)

#
sur_data <-t(cbind(sur4[,c(5,9)],scores))
rownames(sur_data) <-c("OS.time","OS","StromalScore","ImmuneScore","ESTIMATEScore")
gene <-"ImmuneScore"
index <-which(rownames(sur_data)== gene)
dat <-as.data.frame(apply(t(sur_data[c(1,2,index),]),2,as.numeric))
dat$Gene=ifelse(dat[,gene] > median(dat[,gene]),'high','low')
fit <-survfit(Surv(OS.time, OS) ~ Gene, data =dat)
ggsurvplot(fit,pval = TRUE, conf.int = TRUE,
       risk.table = TRUE, # Add risk table
       risk.table.col = "strata", # Change risk table color by groups
       linetype = "strata", # Change line type by groups
       surv.median.line = "hv", ggtheme = theme_bw(), 
       main = "Survival curves",palette = c("#E7B800", "#2E9FDF"),xlim = c(0, max(dat$OS.time)),break.time.by = 365,axes.offset = T,xlab = "Overall Survival (months)",  #x轴的label
       ylab = "Survival Proportion")
gene <-"ESTIMATEScore"
index <-which(rownames(sur_data)== gene)
dat <-as.data.frame(apply(t(sur_data[c(1,2,index),]),2,as.numeric))
dat$Gene=ifelse(dat[,gene] > median(dat[,gene]),'high','low')
fit <-survfit(Surv(OS.time, OS) ~ Gene, data =dat)
ggsurvplot(fit,pval = TRUE, conf.int = TRUE,
       risk.table = TRUE, # Add risk table
       risk.table.col = "strata", # Change risk table color by groups
       linetype = "strata", # Change line type by groups
       surv.median.line = "hv", ggtheme = theme_bw(), 
       main = "Survival curves",palette = c("#E7B800", "#2E9FDF"),xlim = c(0, max(dat$OS.time)),break.time.by = 365,axes.offset = T,xlab = "Overall Survival (months)",  #x轴的label
       ylab = "Survival Proportion")
gene <-"StromalScore"
index <-which(rownames(sur_data)== gene)
dat <-as.data.frame(apply(t(sur_data[c(1,2,index),]),2,as.numeric))
dat$Gene=ifelse(dat[,gene] > median(dat[,gene]),'high','low')
fit <-survfit(Surv(OS.time, OS) ~ Gene, data =dat)
ggsurvplot(fit,pval = TRUE, conf.int = TRUE,
       risk.table = TRUE, # Add risk table
       risk.table.col = "strata", # Change risk table color by groups
       linetype = "strata", # Change line type by groups
       surv.median.line = "hv", ggtheme = theme_bw(), 
       main = "Survival curves",palette = c("#E7B800", "#2E9FDF"),xlim = c(0, max(dat$OS.time)),break.time.by = 365,axes.offset = T,xlab = "Overall Survival (months)",  #x轴的label
       ylab = "Survival Proportion")      

###### 
library(CIBERSORT)
library(ggplot2)
library(dplyr)
library(ggthemes)
library(pheatmap)
library(tibble)
library(tidyr)
library(ggpubr)
library(ggsci)
library(ggthemes)
library(tidyverse)
library(Seurat)
sig_matrix <- system.file("extdata", "LM22.txt", package = "CIBERSORT")
results <- cibersort(sig_matrix, "expr.txt",perm = 500,QN = T)
saveRDS(results,"LM22_500_result.rds")
setwd("./Analysis/Cibersort")
####
OS <- readRDS("./Analysis/data_rPCA_combined.rds") 
Idents(OS) <- "celltype"
X <- AverageExpression(OS)[[1]]  
OS <- FindVariableFeatures(OS, selection.method = "vst", nfeatures = 2000)
gene <- VariableFeatures(OS)
X <-X[gene,]
write.table(X,"sig5.txt",sep = "\t",col.names = T,row.names = T) 
results <- cibersort("sig5.txt", "expr.txt",perm = 500,QN = T)
saveRDS(results,"scRNA_celltype_2000_500_result.rds")
##
results <-readRDS("scRNA_celltype_2000_500_result.rds")
normalize <- function(x) {
  if((max(x) - min(x)) == 0){
    return(mean(x))
  }else{
    return((x - min(x)) / (max(x) - min(x)))
  }
}
heat_map_res <- apply(results[,1:11], 2,normalize)
heat_map_res <- t(as.data.frame(heat_map_res))
mycol<-colorRampPalette(c("navy", "white", "firebrick3"))(100)
library(pheatmap)
pheatmap(heat_map_res,color = mycol,cluster_rows = F,cluster_cols = F,scale = 'none')
##
res4 <- results %>%as.data.frame() %>%rownames_to_column("sample") %>%pivot_longer(cols = 2:11,names_to = "CellType",values_to = "Composition")
res4$CellType <-factor(res4$CellType,levels <-c("Osteoblasts","Chondroblastic","Fibroblast","Pericyte","MSC","Osteoclast","TIL","Myeloid","B cell","Endothelial cell"))
ggbarplot(res4, x = "sample",y = "Composition",size = 0.5,fill = "CellType",color = "CellType") +theme_base() +
  theme(axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1, size = 20),legend.position = "bottom",legend.key.size = unit(10,"pt"))
#
library(survival)
library(survminer)
sur <- read.csv("./Bulk-seq/TARGET-OS.clinical.tsv",sep="\t")
sur[,1]= gsub("-", ".", sur[,1])
sur <-sur[grep("01A",sur[,1]),]
sur1 <-sur[,c(1,6,7,9,10,14,26,30)]
sur2 <-sur1[grep("Alive|Dead",sur1[,4]),]
sur3 <-sur2[which(sur2[,5]!=0),]
event <- c("Alive","Dead")
event1 <- c("0","1")
sur3$event <- plyr::mapvalues(x = sur3[,"Vital.Status"], from = event, to = event1)
rownames(sur3) <-sur3[,1]
sample <- intersect(rownames(results),sur3[,1])
expr <-results[sample,]
sur4 <-sur3[sample,]
sur_data <-as.data.frame(t(cbind(sur4[,c(5,9)],expr)))
rownames(sur_data)[1:2] <-c("OS.time","OS")
data4 <-sur_data
####
gene <-"Osteoblasts"
index <-which(rownames(sur_data)== gene)
dat <-as.data.frame(apply(t(sur_data[c(1,2,index),]),2,as.numeric))
dat$Gene=ifelse(dat[,gene] > median(dat[,gene]),'high','low')
fit <-survfit(Surv(OS.time, OS) ~ Gene, data =dat)
p <- ggsurvplot(fit,pval = TRUE, conf.int = TRUE,
       risk.table = TRUE, # Add risk table
       risk.table.col = "strata", # Change risk table color by groups
       linetype = "strata", # Change line type by groups
       surv.median.line = "hv", ggtheme = theme_bw(), 
       main = "Survival curves",palette = c("#E7B800", "#2E9FDF"),xlim = c(0, max(dat$OS.time)),break.time.by = 365,axes.offset = T,xlab = "Overall Survival (months)",  #x轴的label
       ylab = "Survival Proportion")
pdf(paste(gene,"_sur.pdf"),width=6,height=6)
print(p)
dev.off()
gene <-"Myeloid"
index <-which(rownames(sur_data)== gene)
dat <-as.data.frame(apply(t(sur_data[c(1,2,index),]),2,as.numeric))
dat$Gene=ifelse(dat[,gene] > median(dat[,gene]),'high','low')
fit <-survfit(Surv(OS.time, OS) ~ Gene, data =dat)
p <-  ggsurvplot(fit,pval = TRUE, conf.int = TRUE,
       risk.table = TRUE, # Add risk table
       risk.table.col = "strata", # Change risk table color by groups
       linetype = "strata", # Change line type by groups
       surv.median.line = "hv", ggtheme = theme_bw(), 
       main = "Survival curves",palette = c("#E7B800", "#2E9FDF"),xlim = c(0, max(dat$OS.time)),break.time.by = 365,axes.offset = T,xlab = "Overall Survival (months)",  #x轴的label
       ylab = "Survival Proportion")
pdf(paste(gene,"_sur.pdf"),width=6,height=6)
print(p)
dev.off()
