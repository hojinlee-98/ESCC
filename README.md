# ESCC


## Environment

### Seurat and harmony
```
Seurat : 4.4.0
SeuratObject : 4.1.4
harmony : 0.1.1 
# install.packages("https://cran.r-project.org/src/contrib/Archive/harmony/harmony_0.1.1.tar.gz", repos = NULL, type="source")
```

### font - Arial
```
library(extrafont)
font_import(pattern = "Arial", prompt = FALSE)
loadfonts(device = "pdf")

date = "20250722"
project_name = "escc"

# example
fileidentity <- "no_cellcycle_regressout_umap_res0.9"
cairo_pdf(paste0(date, "_", project_name, "_", fileidentity, "_hj.pdf"), width = 5, height = 5, family="Arial")
set.seed(1234); DimPlot(escc, label = T, group.by = "RNA_snn_res.0.9") & theme_void(base_family = "Arial") & theme(aspect.ratio = 1)
dev.off()

```

## Batch correction for subclustering step
```
library(Seurat)
library(dplyr)
library(ggplot2)
library(harmony)

seurat_obj <- readRDS("20250722_escc_myeloid_hj.rds")

harmony_process <- function(seurat_object, dim, nldr){
  seurat_object <- RunPCA(seurat_object, npcs = dim, verbose = T)
  seurat_object <- RunHarmony(seurat_object, c("donor", "platform"))
  if(nldr == "umap"){
    set.seed(1234); seurat_object <- RunUMAP(seurat_object, reduction = "harmony", dims = 1:dim)
  }
  else if(nldr == "tsne"){
    set.seed(1234); seurat_object <- RunTSNE(seurat_object, reduction = "harmony", dims = 1:dim)
  }
  else{
    stop("Choose a method for nonlinear dimensionality reduction!")
  }
  set.seed(1234); seurat_object <- FindNeighbors(seurat_object, reduction = "harmony", dims = 1:dim)
  return(seurat_object)
}

set.seed(1234); seurat_obj <- NormalizeData(seurat_obj, normalization.method = "LogNormalize", scale.factor = 10000) # log norm
set.seed(1234); seurat_obj@meta.data <- seurat_obj@meta.data %>% dplyr::mutate(donor_platform = paste(donor, platform, sep = "_")) # batch

set.seed(1234); obj_list <- SplitObject(seurat_obj, split.by = "donor_platform") # split
set.seed(1234); obj_list <- obj_list[sapply(obj_list, ncol) >= 10] # exclude small size samples
# find variable gene for each object
set.seed(1234); obj_list <- map(obj_list, function(obj) { 
  set.seed(1234); obj <- NormalizeData(obj)
  set.seed(1234); obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2000)
  return(obj)
})

hvg_all <- obj_list %>% map(VariableFeatures) %>% unlist() %>% unique() # pull hvg
VariableFeatures(seurat_obj) <- hvg_all # set hvg
set.seed(1234); seurat_obj <- ScaleData(seurat_obj, features = hvg_all, split.by = "donor_platform") # scaling
set.seed(1234); seurat_obj <- RunPCA(seurat_obj, features = hvg_all) # pca

ElbowPlot(seurat_obj, ndims = 50) # elbowplot

set.seed(1234); seurat_obj <- harmony_process(seurat_obj, 20, "umap") # PCA, harmony, neighbor, umap
set.seed(1234); seurat_obj <- FindClusters(seurat_obj, resolution = 0.5) # clustering

seurat_obj$RNA_snn_res.0.5 <- factor(seurat_obj$RNA_snn_res.0.5, levels = as.numeric(levels(seurat_obj$RNA_snn_res.0.5)) %>% sort() %>% as.character())

# after subclustering step, the following figures should be checked

mygenes <- c("PTPRC", "TP63", "EGFR", "DSG3", "MUC5B", "KRT7", "AGR2", "CD2", "CD3D", "CD3E", "MS4A1", "CD79A", "BANK1", "CD68", "CD14", "FCER1G", "DCN", "LUM", "COL1A2", "MYH11", "ACTA2", "LMOD1", "PECAM1", "CLDN5", "VWF")

fileidentity <- "lympoid_globalmarkergenes_fp"
cairo_pdf(paste0(date, "_", project_name, "_", fileidentity, "_hj.pdf"), width = 4*5, height = 4*8, family="Arial")
FeaturePlot(seurat_obj, features = mygenes, cols = c("grey", "red"), ncol = 5) & theme(aspect.ratio = 1)
dev.off()

fileidentity <- "lympoid_globalmarkergenes_vp"
cairo_pdf(paste0(date, "_", project_name, "_", fileidentity, "_hj.pdf"), width = 7*5, height = 4*8, family="Arial")
VlnPlot(seurat_obj, features = mygenes, pt.size = 1, ncol = 5) & theme(aspect.ratio = 0.5)
dev.off()

fileidentity <- "quality_check"
cairo_pdf(paste0(date, "_", project_name, "_", fileidentity, "_hj.pdf"), width = 25, height = 4, family = "Arial")
VlnPlot(escc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), pt.size = 0, log = T)
dev.off()

```
