library(Matrix)
library(patchwork)
library(cowplot)
library(BiocManager)
library(multtest)
library(metap)
library(dplyr)
library(data.table)
library(Seurat)
library(ggplot2)
library(clustree)
library(future)
library(dittoSeq)
library(CellChat)
plan(multicore, workers = 8)

earlya <- read.csv(file = "D:/pfc/early/4341_BA46.csv", sep = ",",header = TRUE, row.names = NULL)
names1 = make.unique(earlya$gene_ID)
rownames(earlya) <- names1
earlya <- earlya[,-1] # get rid of old names
earlya <- CreateSeuratObject(counts = earlya, min.cells = 20, project = "earlya")

earlyb <- read.csv(file = "D:/pfc/early/5387_BA9.csv", sep = ",",header = TRUE, row.names = NULL)

names2 = make.unique(earlyb$gene_ID)
rownames(earlyb) <- names2
earlyb <- earlyb[,-1] # get rid of old names
earlyb <- CreateSeuratObject(counts = earlyb, min.cells = 20, project = "earlyb")

earlyc <- read.csv(file = "D:/pfc/early/5936_PFC_Nova.csv", sep = ",",header = TRUE, row.names = NULL)
names3 = make.unique(earlyc$gene_ID)
rownames(earlyc) <- names3
earlyc <- earlyc[,-1] # get rid of old names
earlyc <- CreateSeuratObject(counts = earlyc, min.cells = 20, project = "earlyc")


latea <- read.csv(file = "D:/pfc/late/5538_PFC_Nova.csv", sep = ",",header = TRUE, row.names = NULL)
names4 = make.unique(latea$gene_ID)
rownames(latea) <- names4
latea <- latea[,-1] # get rid of old names
latea <- CreateSeuratObject(counts = latea, min.cells = 20, project = "latea")

lateb <- read.csv(file = "D:/pfc/late/5577_BA9.csv", sep = ",",header = TRUE, row.names = NULL)
names5 = make.unique(lateb$gene_ID)
rownames(lateb) <- names5
lateb <- lateb[,-1] # get rid of old names
lateb <- CreateSeuratObject(counts = lateb, min.cells = 20, project = "lateb")

latec <- read.csv(file = "D:/pfc/late/5958_BA9.csv", sep = ",",header = TRUE, row.names = NULL)
names6 = make.unique(latec$gene_ID)
rownames(latec) <- names6
latec <- latec[,-1] # get rid of old names
latec <- CreateSeuratObject(counts = latec, min.cells = 20, project = "latec")

earlya[["percent.mt"]] = PercentageFeatureSet(earlya, pattern = "^MT-")
earlyb[["percent.mt"]] = PercentageFeatureSet(earlyb, pattern = "^MT-")
earlyc[["percent.mt"]] = PercentageFeatureSet(earlyc, pattern = "^MT-")
latea[["percent.mt"]] = PercentageFeatureSet(latea, pattern = "^MT-")
lateb[["percent.mt"]] = PercentageFeatureSet(lateb, pattern = "^MT-")
latec[["percent.mt"]] = PercentageFeatureSet(latec, pattern = "^MT-")
earlya <- subset(earlya, subset = nFeature_RNA > 500 & nFeature_RNA < 6000 & nCount_RNA < 7000 & percent.mt < 3)
earlyb <- subset(earlyb, subset = nFeature_RNA > 500 & nFeature_RNA < 6000 & nCount_RNA < 7000 & percent.mt < 3)
earlyc <- subset(earlyc, subset = nFeature_RNA > 500 & nFeature_RNA < 6000 & nCount_RNA < 7000 & percent.mt < 3)

latea <- subset(latea, subset = nFeature_RNA > 500 & nFeature_RNA < 6000 & nCount_RNA < 7000 & percent.mt < 3)
lateb <- subset(lateb, subset = nFeature_RNA > 500 & nFeature_RNA < 6000 & nCount_RNA < 7000 & percent.mt < 3)
latec <- subset(latec, subset = nFeature_RNA > 500 & nFeature_RNA < 6000 & nCount_RNA < 7000 & percent.mt < 3)

VlnPlot(earlya, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"))
VlnPlot(earlyb, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"))
VlnPlot(earlyc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"))
VlnPlot(latea, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"))
VlnPlot(lateb, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"))
VlnPlot(latec, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"))
earlya$sample <- "Early"
earlyb$sample <- "Early"
earlyc$sample <- "Early"
latea$sample <- "Late"
lateb$sample <- "Late"
latec$sample <- "Late"
earlya <- NormalizeData(earlya, verbose = FALSE)
earlya <- FindVariableFeatures(earlya, selection.method = "vst", nfeatures = 3000)
earlyb <- NormalizeData(earlyb, verbose = FALSE)
earlyb <- FindVariableFeatures(earlyb, selection.method = "vst", nfeatures = 3000)
earlyc <- NormalizeData(earlyc, verbose = FALSE)
earlyc <- FindVariableFeatures(earlyc, selection.method = "vst", nfeatures = 3000)
early.anchors <- FindIntegrationAnchors(object.list = list(earlya, earlyb, earlyc), dims = 1:20)
early <- IntegrateData(anchorset = early.anchors, dims = 1:20)
latea <- NormalizeData(latea, verbose = FALSE)
latea <- FindVariableFeatures(latea, selection.method = "vst", nfeatures = 3000)
lateb <- NormalizeData(lateb, verbose = FALSE)
lateb <- FindVariableFeatures(lateb, selection.method = "vst", nfeatures = 3000)
latec <- NormalizeData(latec, verbose = FALSE)
latec <- FindVariableFeatures(latec, selection.method = "vst", nfeatures = 3000)
late.anchors <- FindIntegrationAnchors(object.list = list(latea, lateb, latec), dims = 1:20)
late <- IntegrateData(anchorset = late.anchors, dims = 1:20)


pfc.anchors <- FindIntegrationAnchors(object.list = list(early, late), dims = 1:20)
pfc <- IntegrateData(anchorset = pfc.anchors, dims = 1:20)
DefaultAssay(pfc) <- "integrated"
pfc <- ScaleData(pfc, verbose = FALSE)
pfc <- RunPCA(pfc, npcs = 30, verbose = FALSE)
VizDimLoadings(pfc, dims = 1:2, reduction = "pca")
DimHeatmap(pfc, dims = 1:12, cells = 500, balanced = TRUE)
pfc <- JackStraw(pfc, num.replicate = 100)
pfc <- ScoreJackStraw(pfc, dims = 1:20)
ElbowPlot(pfc)

pfc <- RunUMAP(pfc, reduction = "pca", dims = 1:15)
pfc <- FindNeighbors(pfc, reduction = "pca", dims = 1:15)
res.used <- seq(0.1,2,by=0.2)
res.used
# Loop over and perform clustering of different resolutions
for(i in res.used){
  sce <- FindClusters(object = pfc , verbose = T, resolution = res.used)
}
clus.tree.out <- clustree(sce) +
  theme(legend.position = "bottom") +
  scale_color_brewer(palette = "Set1") +
  scale_edge_color_continuous(low = "grey80", high = "red")
clus.tree.out
library(clustree)
clus.tree.out <- clustree(sce) +
  theme(legend.position = "bottom") +
  scale_color_brewer(palette = "Set1") +
  scale_edge_color_continuous(low = "grey80", high = "red")
clus.tree.out

pfc <- FindClusters(pfc, resolution = 1.4)
DimPlot(pfc, reduction = "umap", label = TRUE)+theme(aspect.ratio = 1)
DimPlot(pfc, reduction = "umap", split.by = "sample")+theme(aspect.ratio = 1)
DefaultAssay(pfc) = "RNA"
FeaturePlot(pfc, features = c("SATB2","RORB", "CAMK2B", "GAD1","GAD2","SST", "VIP"), min.cutoff = "q16")
ggsave("D:/pfc/neuron-markers.png", dpi = 300, width = 90, height =90, unit = "cm")
FeaturePlot(pfc, features = c("AQP4", "GJA1","CD74",  "P2RY12", "MOG","PLP1", "CCND1", "GPR17", "PDGFRA"), min.cutoff = "q16")
ggsave("D:/pfc/glia-marker.png", dpi = 300, width = 60, height =60, unit = "cm")
pfc = readRDS("D:/pfc-for-pub.rds")
pfc <- RenameIdents(pfc, '5'= "Ast", '16' = "Ast", '27'= "Ast", '21' = "Ast", '2' = "Olig", '4' = "Olig", '14' = "Olig", '0' = "ExN", '3' = "ExN", '7' = "ExN", '12' = "ExN", '13' = "ExN",'11' = "ExN", '15' = "ExN",  '26' = "ExN",  '6' = "InN", '9' = "InN",'10' = "InN",'19' = "InN",'20' = "InN",'22' = "InN",'24' = "InN",'8' = "Mic", '29' = "Mic")

DimPlot(pfc, reduction = "umap", label = TRUE, label.size = 8)+xlim(-18, 18)+ylim(-15, 20)+NoLegend()+theme(aspect.ratio = 1)
cells = subset(pfc, idents = c("Ast", "Mic", "Olig", "InN", "ExN"))
DimPlot(cells, reduction = "umap", label = TRUE, label.size = 14)+xlim(-18, 18)+ylim(-18, 20)+NoLegend()+theme(aspect.ratio = 1)

pfc$celltype.sample <- paste(Idents(pfc), pfc$sample, sep = "_")
pfc$celltype <- Idents(pfc)
Idents(pfc) <- "celltype.sample"
pfc = readRDS("D:/pfc-for-pub2.rds")
astmarker<- FindMarkers(pfc, ident.1 = "Ast_late", ident.2 = "Ast_early", min.pct = 0.1, logfc.threshold = 0, verbose = FALSE)
write.csv(astmarker, file = "D:/pfc/astmarker-latevsearly.csv")

oligmarker<- FindMarkers(pfc, ident.1 = "Olig_Late", ident.2 = "Olig_Early", min.pct = 0.1, logfc.threshold = 0, verbose = FALSE)
write.csv(oligmarker, file = "D:/pfc/oligmarker-latevsearly.csv")
exnmarker<- FindMarkers(pfc, ident.1 = "ExN_Late", ident.2 = "ExN_Early", min.pct = 0.1, logfc.threshold = 0, verbose = FALSE)
write.csv(exnmarker, file = "D:/pfc/exnmarker-latevsearly.csv")
Micmarket<- FindMarkers(pfc, ident.1 = "Mic_Late", ident.2 = "Mic_Early", min.pct = 0.1, logfc.threshold = 0, verbose = FALSE)
write.csv(Micmarker, file = "D:/pfc/Micmarker-latevsearly.csv")
InNmarker<- FindMarkers(pfc, ident.1 = "InN_Late", ident.2 = "InN_Early", min.pct = 0.1, logfc.threshold = 0, verbose = FALSE)
write.csv(InNmarker, file = "D:/pfc/InNmarker-latevsearly.csv")

ast = subset(pfc, idents = c("Ast_Early", "Ast_Late"))
dittoPlot(ast, "SAMD4A", group.by = "ident", boxplot.width = 0.6,jitter.width = 1.2, jitter.size = 1,
          plots = c("boxplot", "jitter"))+theme(legend.position = "upper", axis.title.x = element_blank(), axis.text.y = element_text(size = 50), axis.title.y = element_text(size = rel(4), hjust=0.5), axis.text.x = element_text(angle = 50, size = 50), aspect.ratio = 2, plot.title = element_text(size = rel(4), vjust = 1, hjust = 0.4), axis.ticks.x = element_blank())+ ylab("Avg_log2(FC)")
ggsave("D:/ast-samd4a.png", dpi = 300, height = 30, width = 15, units = "cm")

dittoPlot(ast, "CCL2", group.by = "ident", boxplot.width = 0.65,jitter.width = 1.2,jitter.size = 1,
          plots = c("boxplot", "jitter"))+theme(legend.position = "upper", axis.title.x = element_blank(), axis.text.y = element_text(size = 50), axis.title.y = element_text(size = rel(4), hjust=0.5), axis.text.x = element_text(angle = 50, size = 50), aspect.ratio = 2, plot.title = element_text(size = rel(4), vjust = 1, hjust = 0.4), axis.ticks.x = element_blank())+ ylab("Avg_log2(FC)")
ggsave("D:/ast-ccl2.png", dpi = 300, height = 30, width = 15, units = "cm")

dittoPlot(ast, "NRXN1", group.by = "ident", boxplot.width = 0.6,jitter.width = 1.2,jitter.size = 1,
          plots = c("boxplot", "jitter"))+theme(legend.position = "upper", axis.title.x = element_blank(), axis.text.y = element_text(size = 50), axis.title.y = element_text(size = rel(4), hjust=0.5), axis.text.x = element_text(angle = 50, size = 50), aspect.ratio = 2, plot.title = element_text(size = rel(4), vjust = 1, hjust = 0.4), axis.ticks.x = element_blank())+ ylab("Avg_log2(FC)")
ggsave("D:/ast-nrxn1.png", dpi = 300, height = 30, width = 15, units = "cm")

dittoPlot(ast, "GPR37L1", group.by = "ident", boxplot.width = 0.6,jitter.width = 1.2, jitter.size = 1,
          plots = c("boxplot", "jitter"))+theme(legend.position = "upper", axis.title.x = element_blank(), axis.text.y = element_text(size = 50), axis.title.y = element_text(size = rel(4), hjust=0.5), axis.text.x = element_text(angle = 50, size = 50), aspect.ratio = 2, plot.title = element_text(size = rel(4), vjust = 1, hjust = 0.4), axis.ticks.x = element_blank())+ ylab("Avg_log2(FC)")
ggsave("D:/ast-GPR37L1.png", dpi = 300, height = 30, width = 15, units = "cm")

dittoPlot(ast, "PSAP", group.by = "ident", boxplot.width = 0.65,jitter.width = 1.2,jitter.size = 1.2,
          plots = c("boxplot", "jitter"))+theme(legend.position = "upper", axis.title.x = element_blank(), axis.text.y = element_text(size = 50), axis.title.y = element_text(size = rel(4), hjust=0.5), axis.text.x = element_text(angle = 50, size = 50), aspect.ratio = 2, plot.title = element_text(size = rel(4), vjust = 1, hjust = 0.4), axis.ticks.x = element_blank())+ ylab("Avg_log2(FC)")
ggsave("D:/ast-PASP.png", dpi = 300, height = 30, width = 15, units = "cm")

dittoPlot(ast, "ERBB4", group.by = "ident", boxplot.width = 0.6,jitter.width = 1.2,jitter.size = 1,
          plots = c("boxplot", "jitter"))+theme(legend.position = "upper", axis.title.x = element_blank(), axis.text.y = element_text(size = 50), axis.title.y = element_text(size = rel(4), hjust=0.5), axis.text.x = element_text(angle = 50, size = 50), aspect.ratio = 2, plot.title = element_text(size = rel(4), vjust = 1, hjust = 0.4), axis.ticks.x = element_blank())+ ylab("Avg_log2(FC)")
ggsave("D:/ast-ERBB4.png", dpi = 300, height = 30, width = 15, units = "cm")

dittoPlot(ast, "NRG3", group.by = "ident", boxplot.width = 0.6,jitter.width = 1.2, jitter.size = 1,
          plots = c("boxplot", "jitter"))+theme(legend.position = "upper", axis.title.x = element_blank(), axis.text.y = element_text(size = 50), axis.title.y = element_text(size = rel(4), hjust=0.5), axis.text.x = element_text(angle = 50, size = 50), aspect.ratio = 2, plot.title = element_text(size = rel(4), vjust = 1, hjust = 0.4), axis.ticks.x = element_blank())+ ylab("Avg_log2(FC)")
ggsave("D:/ast-NRG3.png", dpi = 300, height = 30, width = 15, units = "cm")

dittoPlot(ast, "SDC3", group.by = "ident", boxplot.width = 0.65,jitter.width = 1.2,jitter.size = 1.2,
          plots = c("boxplot", "jitter"))+theme(legend.position = "upper", axis.title.x = element_blank(), axis.text.y = element_text(size = 50), axis.title.y = element_text(size = rel(4), hjust=0.5), axis.text.x = element_text(angle = 50, size = 50), aspect.ratio = 2, plot.title = element_text(size = rel(4), vjust = 1, hjust = 0.4), axis.ticks.x = element_blank())+ ylab("Avg_log2(FC)")
ggsave("D:/ast-SDC3.png", dpi = 300, height = 30, width = 15, units = "cm")

dittoPlot(ast, "TGFB2", group.by = "ident", boxplot.width = 0.6,jitter.width = 1.2,jitter.size = 1,
          plots = c("boxplot", "jitter"))+theme(legend.position = "upper", axis.title.x = element_blank(), axis.text.y = element_text(size = 50), axis.title.y = element_text(size = rel(4), hjust=0.5), axis.text.x = element_text(angle = 50, size = 50), aspect.ratio = 2, plot.title = element_text(size = rel(4), vjust = 1, hjust = 0.4), axis.ticks.x = element_blank())+ ylab("Avg_log2(FC)")
ggsave("D:/ast-TGFB2.png", dpi = 300, height = 30, width = 15, units = "cm")

dittoPlot(ast, "FGF2", group.by = "ident", boxplot.width = 0.6,jitter.width = 1.2, jitter.size = 1.2,
          plots = c("boxplot", "jitter"))+theme(legend.position = "upper", axis.title.x = element_blank(), axis.text.y = element_text(size = 50), axis.title.y = element_text(size = rel(4), hjust=0.5), axis.text.x = element_text(angle = 50, size = 50), aspect.ratio = 2, plot.title = element_text(size = rel(4), vjust = 1, hjust = 0.4), axis.ticks.x = element_blank())+ ylab("Avg_log2(FC)")
ggsave("D:/ast-FGF2.png", dpi = 300, height = 30, width = 15, units = "cm")

dittoPlot(ast, "IGF1R", group.by = "ident", boxplot.width = 0.65,jitter.width = 1.2,jitter.size = 1,
          plots = c("boxplot", "jitter"))+theme(legend.position = "upper", axis.title.x = element_blank(), axis.text.y = element_text(size = 50), axis.title.y = element_text(size = rel(4), hjust=0.5), axis.text.x = element_text(angle = 50, size = 50), aspect.ratio = 2, plot.title = element_text(size = rel(4), vjust = 1, hjust = 0.4), axis.ticks.x = element_blank())+ ylab("Avg_log2(FC)")
ggsave("D:/ast-IGF1R.png", dpi = 300, height = 30, width = 15, units = "cm")

dittoPlot(ast, "ANGPTL4", group.by = "ident", boxplot.width = 0.6,jitter.width = 1.2,jitter.size = 1,
          plots = c("boxplot", "jitter"))+theme(legend.position = "upper", axis.title.x = element_blank(), axis.text.y = element_text(size = 50), axis.title.y = element_text(size = rel(4), hjust=0.5), axis.text.x = element_text(angle = 50, size = 50), aspect.ratio = 2, plot.title = element_text(size = rel(4), vjust = 1, hjust = 0.4), axis.ticks.x = element_blank())+ ylab("Avg_log2(FC)")
ggsave("D:/ast-ANGPTL4.png", dpi = 300, height = 30, width = 15, units = "cm")

dittoPlot(ast, "FGFR3", group.by = "ident", boxplot.width = 0.6,jitter.width = 1.2,jitter.size = 1.2,
          plots = c("boxplot", "jitter"))+theme(legend.position = "upper", axis.title.x = element_blank(), axis.text.y = element_text(size = 50), axis.title.y = element_text(size = rel(4), hjust=0.5), axis.text.x = element_text(angle = 50, size = 50), aspect.ratio = 2, plot.title = element_text(size = rel(4), vjust = 1, hjust = 0.4), axis.ticks.x = element_blank())+ ylab("Avg_log2(FC)")
ggsave("D:/ast-FGFR3.png", dpi = 300, height = 30, width = 15, units = "cm")

exn = subset(pfc, idents = c("ExN_Early", "ExN_Late"))
dittoPlot(exn, "EPHA6", group.by = "ident", boxplot.width = 0.6,jitter.width = 1.2,jitter.size = 1,
          plots = c("boxplot", "jitter"))+theme(legend.position = "upper", axis.title.x = element_blank(), axis.text.y = element_text(size = 50), axis.title.y = element_text(size = rel(5), hjust=0.5), axis.text.x = element_text(angle = 50, size = 50), aspect.ratio = 2, plot.title = element_text(size = rel(4), vjust = 1, hjust = 0.4), axis.ticks.x = element_blank())+ ylab("Avg_log2(FC)")
ggsave("D:/exm-epha6b.png", dpi = 300, height = 30, width = 15, units = "cm")

dittoPlot(exn, "NPAS4", group.by = "ident", boxplot.width = 0.6,jitter.width = 1.2,jitter.size = 1,
          plots = c("boxplot", "jitter"))+theme(legend.position = "upper", axis.title.x = element_blank(), axis.text.y = element_text(size = 50), axis.title.y = element_text(size = rel(5), hjust=0.5), axis.text.x = element_text(angle = 50, size = 50), aspect.ratio = 2, plot.title = element_text(size = rel(4), vjust = 1, hjust = 0.4), axis.ticks.x = element_blank())+ ylab("Avg_log2(FC)")
ggsave("D:/exm-npas4.png", dpi = 300, height = 30, width = 15, units = "cm")

dittoPlot(exn, "KIRREL3", group.by = "ident", boxplot.width = 0.6,jitter.width = 1.2,jitter.size = 1,
          plots = c("boxplot", "jitter"))+theme(legend.position = "upper", axis.title.x = element_blank(), axis.text.y = element_text(size = 50), axis.title.y = element_text(size = rel(5), hjust=0.5), axis.text.x = element_text(angle = 50, size = 50), aspect.ratio = 2, plot.title = element_text(size = rel(4), vjust = 1, hjust = 0.4), axis.ticks.x = element_blank())+ ylab("Avg_log2(FC)")
ggsave("D:/exm-kirrel3.png", dpi = 300, height = 30, width = 15, units = "cm")
dittoPlot(exn, "NWD2", group.by = "ident", boxplot.width = 0.6,jitter.width = 1.2,jitter.size = 1,
          plots = c("boxplot", "jitter"))+theme(legend.position = "upper", axis.title.x = element_blank(), axis.text.y = element_text(size = 50), axis.title.y = element_text(size = rel(5), hjust=0.5), axis.text.x = element_text(angle = 50, size = 50), aspect.ratio = 2, plot.title = element_text(size = rel(4), vjust = 1, hjust = 0.4), axis.ticks.x = element_blank())+ ylab("Avg_log2(FC)")
ggsave("D:/exm-nwd2.png", dpi = 300, height = 30, width = 15, units = "cm")



dittoPlot(exn, "CX3CL1", group.by = "ident", boxplot.width = 0.6,jitter.width = 1.2,jitter.size = 1,
          plots = c("boxplot", "jitter"))+theme(legend.position = "upper", axis.title.x = element_blank(), axis.text.y = element_text(size = 50), axis.title.y = element_text(size = rel(5), hjust=0.5), axis.text.x = element_text(angle = 50, size = 50), aspect.ratio = 2, plot.title = element_text(size = rel(4), vjust = 1, hjust = 0.4), axis.ticks.x = element_blank())+ ylab("Avg_log2(FC)")
ggsave("D:/exm-CX3CL1.png", dpi = 300, height = 30, width = 15, units = "cm")

dittoPlot(exn, "SDC3", group.by = "ident", boxplot.width = 0.6,jitter.width = 1.2,jitter.size = 1,
          plots = c("boxplot", "jitter"))+theme(legend.position = "upper", axis.title.x = element_blank(), axis.text.y = element_text(size = 50), axis.title.y = element_text(size = rel(5), hjust=0.5), axis.text.x = element_text(angle = 50, size = 50), aspect.ratio = 2, plot.title = element_text(size = rel(4), vjust = 1, hjust = 0.4), axis.ticks.x = element_blank())+ ylab("Avg_log2(FC)")
ggsave("D:/exm-SDC3.png", dpi = 300, height = 30, width = 15, units = "cm")

dittoPlot(exn, "PTPRZ1", group.by = "ident", boxplot.width = 0.6,jitter.width = 1.2,jitter.size = 1,
          plots = c("boxplot", "jitter"))+theme(legend.position = "upper", axis.title.x = element_blank(), axis.text.y = element_text(size = 50), axis.title.y = element_text(size = rel(5), hjust=0.5), axis.text.x = element_text(angle = 50, size = 50), aspect.ratio = 2, plot.title = element_text(size = rel(4), vjust = 1, hjust = 0.4), axis.ticks.x = element_blank())+ ylab("Avg_log2(FC)")
ggsave("D:/exm-PTPRZ14.png", dpi = 300, height = 30, width = 15, units = "cm")

dittoPlot(exn, "CDH11", group.by = "ident", boxplot.width = 0.6,jitter.width = 1.2,jitter.size = 1,
          plots = c("boxplot", "jitter"))+theme(legend.position = "upper", axis.title.x = element_blank(), axis.text.y = element_text(size = 50), axis.title.y = element_text(size = rel(5), hjust=0.5), axis.text.x = element_text(angle = 50, size = 50), aspect.ratio = 2, plot.title = element_text(size = rel(4), vjust = 1, hjust = 0.4), axis.ticks.x = element_blank())+ ylab("Avg_log2(FC)")
ggsave("D:/exm-CDH11.png", dpi = 300, height = 30, width = 15, units = "cm")

dittoPlot(exn, "NRG1", group.by = "ident", boxplot.width = 0.6,jitter.width = 1.2,jitter.size = 1,
          plots = c("boxplot", "jitter"))+theme(legend.position = "upper", axis.title.x = element_blank(), axis.text.y = element_text(size = 50), axis.title.y = element_text(size = rel(5), hjust=0.5), axis.text.x = element_text(angle = 50, size = 50), aspect.ratio = 2, plot.title = element_text(size = rel(4), vjust = 1, hjust = 0.4), axis.ticks.x = element_blank())+ ylab("Avg_log2(FC)")
ggsave("D:/exm-NRG1.png", dpi = 300, height = 30, width = 15, units = "cm")
dittoPlot(exn, "NRG3", group.by = "ident", boxplot.width = 0.6,jitter.width = 1.2,jitter.size = 1,
          plots = c("boxplot", "jitter"))+theme(legend.position = "upper", axis.title.x = element_blank(), axis.text.y = element_text(size = 50), axis.title.y = element_text(size = rel(5), hjust=0.5), axis.text.x = element_text(angle = 50, size = 50), aspect.ratio = 2, plot.title = element_text(size = rel(4), vjust = 1, hjust = 0.4), axis.ticks.x = element_blank())+ ylab("Avg_log2(FC)")
ggsave("D:/exm-NRG3.png", dpi = 300, height = 30, width = 15, units = "cm")

mic = subset(pfc, idents = c("Mic_Early", "Mic_Late"))
dittoPlot(mic, "RGS1", group.by = "ident", boxplot.width = 0.6,jitter.width = 1.2,jitter.size = 1.2,
          plots = c("boxplot", "jitter"))+theme(legend.position = "upper", axis.title.x = element_blank(), axis.text.y = element_text(size = 50), axis.title.y = element_text(size = rel(5), hjust=0.5), axis.text.x = element_text(angle = 50, size = 50), aspect.ratio = 2, plot.title = element_text(size = rel(4), vjust = 1, hjust = 0.4), axis.ticks.x = element_blank())+ ylab("Avg_log2(FC)")
ggsave("D:/mic-rgs1.png", dpi = 300, height = 30, width = 15, units = "cm")

dittoPlot(mic, "RASGEF1B", group.by = "ident", boxplot.width = 0.6,jitter.width = 1.2,jitter.size = 1.2,
          plots = c("boxplot", "jitter"))+theme(legend.position = "upper", axis.title.x = element_blank(), axis.text.y = element_text(size = 50), axis.title.y = element_text(size = rel(5), hjust=0.5), axis.text.x = element_text(angle = 50, size = 50), aspect.ratio = 2, plot.title = element_text(size = rel(4), vjust = 1, hjust = 0.4), axis.ticks.x = element_blank())+ ylab("Avg_log2(FC)")
ggsave("D:/mic-RASGEF1B.png", dpi = 300, height = 30, width = 15, units = "cm")

dittoPlot(mic, "SPP1", group.by = "ident", boxplot.width = 0.6,jitter.width = 1.2,jitter.size = 1.2,
          plots = c("boxplot", "jitter"))+theme(legend.position = "upper", axis.title.x = element_blank(), axis.text.y = element_text(size = 50), axis.title.y = element_text(size = rel(5), hjust=0.5), axis.text.x = element_text(angle = 50, size = 50), aspect.ratio = 2, plot.title = element_text(size = rel(4), vjust = 1, hjust = 0.4), axis.ticks.x = element_blank())+ ylab("Avg_log2(FC)")
ggsave("D:/mic-SPP1.png", dpi = 300, height = 30, width = 15, units = "cm")


olig = subset(pfc, idents = c("Olig_Early", "Olig_Late"))
dittoPlot(olig, "CNDP1", group.by = "ident", boxplot.width = 0.6,jitter.width = 1.2,jitter.size = 1.2,
          plots = c("boxplot", "jitter"))+theme(legend.position = "upper", axis.title.x = element_blank(), axis.text.y = element_text(size = 50), axis.title.y = element_text(size = rel(5), hjust=0.5), axis.text.x = element_text(angle = 50, size = 50), aspect.ratio = 2, plot.title = element_text(size = rel(4), vjust = 1, hjust = 0.4), axis.ticks.x = element_blank())+ ylab("Avg_log2(FC)")
ggsave("D:/olig-cndp1.png", dpi = 300, height = 30, width = 15, units = "cm")

dittoPlot(olig, "TTLL7", group.by = "ident", boxplot.width = 0.6,jitter.width = 1.2,jitter.size = 1,
          plots = c("boxplot", "jitter"))+theme(legend.position = "upper", axis.title.x = element_blank(), axis.text.y = element_text(size = 50), axis.title.y = element_text(size = rel(5), hjust=0.5), axis.text.x = element_text(angle = 50, size = 50), aspect.ratio = 2, plot.title = element_text(size = rel(4), vjust = 1, hjust = 0.4), axis.ticks.x = element_blank())+ ylab("Avg_log2(FC)")
ggsave("D:/olig-ttll7.png", dpi = 300, height = 30, width = 15, units = "cm")

dittoPlot(olig, "RTN4", group.by = "ident", boxplot.width = 0.6,jitter.width = 1.2,jitter.size = 1.2,
          plots = c("boxplot", "jitter"))+theme(legend.position = "upper", axis.title.x = element_blank(), axis.text.y = element_text(size = 50), axis.title.y = element_text(size = rel(5), hjust=0.5), axis.text.x = element_text(angle = 50, size = 50), aspect.ratio = 2, plot.title = element_text(size = rel(4), vjust = 1, hjust = 0.4), axis.ticks.x = element_blank())+ ylab("Avg_log2(FC)")
ggsave("D:/olig-rtn4.png", dpi = 300, height = 30, width = 15, units = "cm")

dittoPlot(olig, "ERBB4", group.by = "ident", boxplot.width = 0.6,jitter.width = 1.2,jitter.size = 1.2,
          plots = c("boxplot", "jitter"))+theme(legend.position = "upper", axis.title.x = element_blank(), axis.text.y = element_text(size = 50), axis.title.y = element_text(size = rel(5), hjust=0.5), axis.text.x = element_text(angle = 50, size = 50), aspect.ratio = 2, plot.title = element_text(size = rel(4), vjust = 1, hjust = 0.4), axis.ticks.x = element_blank())+ ylab("Avg_log2(FC)")
ggsave("D:/olig-ERBB4.png", dpi = 300, height = 30, width = 15, units = "cm")

inn = subset(pfc, idents = c("InN_Early", "InN_Late"))
dittoPlot(inn, "HSPH1", group.by = "ident", boxplot.width = 0.6,jitter.width = 1.2,jitter.size = 1.2,
          plots = c("boxplot", "jitter"))+theme(legend.position = "upper", axis.title.x = element_blank(), axis.text.y = element_text(size = 50), axis.title.y = element_text(size = rel(5), hjust=0.5), axis.text.x = element_text(angle = 50, size = 50), aspect.ratio = 2, plot.title = element_text(size = rel(4), vjust = 1, hjust = 0.4), axis.ticks.x = element_blank())+ ylab("Avg_log2(FC)")
ggsave("D:/inn-hsph1.png", dpi = 300, height = 30, width = 15, units = "cm")

dittoPlot(inn, "FTH1", group.by = "ident", boxplot.width = 0.6,jitter.width = 1.2,jitter.size = 1.2,
          plots = c("boxplot", "jitter"))+theme(legend.position = "upper", axis.title.x = element_blank(), axis.text.y = element_text(size = 50), axis.title.y = element_text(size = rel(5), hjust=0.5), axis.text.x = element_text(angle = 50, size = 50), aspect.ratio = 2, plot.title = element_text(size = rel(4), vjust = 1, hjust = 0.4), axis.ticks.x = element_blank())+ ylab("Avg_log2(FC)")
ggsave("D:/inn-fth1.png", dpi = 300, height = 30, width = 15, units = "cm")

dittoPlot(inn, "CCDC85B", group.by = "ident", boxplot.width = 0.6,jitter.width = 1.2,jitter.size = 1.2,
          plots = c("boxplot", "jitter"))+theme(legend.position = "upper", axis.title.x = element_blank(), axis.text.y = element_text(size = 50), axis.title.y = element_text(size = rel(5), hjust=0.5), axis.text.x = element_text(angle = 50, size = 50), aspect.ratio = 2, plot.title = element_text(size = rel(4), vjust = 1, hjust = 0.4), axis.ticks.x = element_blank())+ ylab("Avg_log2(FC)")
ggsave("D:/inn-ccdc85b.png", dpi = 300, height = 30, width = 15, units = "cm")

dittoPlot(inn, "ERBB4", group.by = "ident", boxplot.width = 0.6,jitter.width = 1.2,jitter.size = 1.2,
          plots = c("boxplot", "jitter"))+theme(legend.position = "upper", axis.title.x = element_blank(), axis.text.y = element_text(size = 50), axis.title.y = element_text(size = rel(5), hjust=0.5), axis.text.x = element_text(angle = 50, size = 50), aspect.ratio = 2, plot.title = element_text(size = rel(4), vjust = 1, hjust = 0.4), axis.ticks.x = element_blank())+ ylab("Avg_log2(FC)")
ggsave("D:/inn-ERBB4.png", dpi = 300, height = 30, width = 15, units = "cm")

dittoPlot(inn, "NRG3", group.by = "ident", boxplot.width = 0.6,jitter.width = 1.2,jitter.size = 1.2,
          plots = c("boxplot", "jitter"))+theme(legend.position = "upper", axis.title.x = element_blank(), axis.text.y = element_text(size = 50), axis.title.y = element_text(size = rel(5), hjust=0.5), axis.text.x = element_text(angle = 50, size = 50), aspect.ratio = 2, plot.title = element_text(size = rel(4), vjust = 1, hjust = 0.4), axis.ticks.x = element_blank())+ ylab("Avg_log2(FC)")
ggsave("D:/inn-NRG3.png", dpi = 300, height = 30, width = 15, units = "cm")

dittoPlot(inn, "ALK", group.by = "ident", boxplot.width = 0.6,jitter.width = 1.2,jitter.size = 1.2,
          plots = c("boxplot", "jitter"))+theme(legend.position = "upper", axis.title.x = element_blank(), axis.text.y = element_text(size = 50), axis.title.y = element_text(size = rel(5), hjust=0.5), axis.text.x = element_text(angle = 50, size = 50), aspect.ratio = 2, plot.title = element_text(size = rel(4), vjust = 1, hjust = 0.4), axis.ticks.x = element_blank())+ ylab("Avg_log2(FC)")
ggsave("D:/inn-ALK.png", dpi = 300, height = 30, width = 15, units = "cm")



earlyado = subset(pfc, idents = c("Olig_Early", "Ast_Early", "InN_Early", "ExN_Early", "Mic_Early"))
lateado = subset(pfc, idents = c("Olig_Late", "Ast_Late", "InN_Late", "ExN_Late", "Mic_Late"))

data.input <- GetAssayData(earlyado, assay = "RNA", slot = "data") # normalized data matrix
labels <- Idents(earlyado)
identity <- data.frame(group = labels, row.names = names(labels)) # create a dataframe of the cell labels
CCearly <- createCellChat(data.input)
CCearly <- addMeta(CCearly, meta = identity, meta.name = "labels")
CCearly <- setIdent(CCearly, ident.use = "labels") # set "labels" as default cell identity
CellChatDB <- CellChatDB.human # use CellChatDB.human if running on human data
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling")
CCearly@DB <- CellChatDB.use # set the used database in the object
CCearly <- subsetData(CCearly) # subset the expression data of signaling genes for saving computation cost
future::plan("multiprocess", workers = 8) # do parallel
CCearly <- identifyOverExpressedGenes(CCearly)
CCearly <- identifyOverExpressedInteractions(CCearly)
CCearly <- projectData(CCearly, PPI.human)
groupSize <- as.numeric(table(CCearly@idents))
CCearly <- computeCommunProb(CCearly)
CCearly <- computeCommunProbPathway(CCearly)
CClate <- aggregateNet(CCearly)


data.input <- GetAssayData(lateado, assay = "RNA", slot = "data") # normalized data matrix
labels <- Idents(lateado)
identity <- data.frame(group = labels, row.names = names(labels)) # create a dataframe of the cell labels
CClate <- createCellChat(data.input)
CClate <- addMeta(CClate, meta = identity, meta.name = "labels")
CClate <- setIdent(CClate, ident.use = "labels") # set "labels" as default cell identity
CellChatDB <- CellChatDB.human # use CellChatDB.human if running on human data
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling")
CClate@DB <- CellChatDB.use # set the used database in the object
CClate <- subsetData(CClate) # subset the expression data of signaling genes for saving computation cost
future::plan("multiprocess", workers = 8) # do parallel
CClate <- identifyOverExpressedGenes(CClate)
CClate <- identifyOverExpressedInteractions(CClate)
CClate <- projectData(CClate, PPI.human)
groupSize <- as.numeric(table(CClate@idents))
CClate <- computeCommunProb(CClate)
CClate <- computeCommunProbPathway(CClate)
CClate <- aggregateNet(CClate)

pathways.show <- c("FGF") 
par(mfrow=c(1,1))
netVisual_aggregate(CClate, signaling = pathways.show, layout = "chord")

pathways.show <- c("SPP1") 
par(mfrow=c(1,1))
netVisual_aggregate(CClate, signaling = pathways.show, layout = "chord")

pathways.show <- c("ANGPTL") 
par(mfrow=c(1,1))
netVisual_aggregate(CClate, signaling = pathways.show, layout = "chord")

pathways.show <- c("TGFb") 
par(mfrow=c(1,1))
netVisual_aggregate(CClate, signaling = pathways.show, layout = "chord")


CCearly.ad<-netAnalysis_computeCentrality(CCearly)
CClate.ad<-netAnalysis_computeCentrality(CClate)
object.list <- list(early = CCearly.ad, late = CClate.ad)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))

gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2))
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight")
gg1 + gg2



rankNet(cellchat, mode = "comparison", stacked = T, do.stat = TRUE)

pathways.show <- c("IGF") 
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "chord", signaling.name = paste(pathways.show, names(object.list)[i]))
}

pathways.show <- c("PSAP") 
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "chord", signaling.name = paste(pathways.show, names(object.list)[i]))
}

pathways.show <- c("NRG") 
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "chord", signaling.name = paste(pathways.show, names(object.list)[i]))
}

pathways.show <- c("PTN") 
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "chord", signaling.name = paste(pathways.show, names(object.list)[i]))
}

pathways.show <- c("CX3C") 
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "chord", signaling.name = paste(pathways.show, names(object.list)[i]))
}

par(mfrow = c(2,3), xpd=TRUE)
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "count.merged", label.edge = T)
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight.merged", label.edge = T)

