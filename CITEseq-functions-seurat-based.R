library(Seurat)
library(tidyverse)
library(SingleR)
library(SingleCellExperiment)
#==============================================
#ST normalize data: Seurat workflow - ST transform for RNA and CLR protein and do PCA for both
#this is what is recommend in the vignette . Regress ccell cycle scores too
#===================================================

Seurat.STnorm.pca.CITEseq <- function(SeuratObj){
  
  s.genes <- cc.genes$s.genes
  g2m.genes <- cc.genes$g2m.genes

  
  DefaultAssay(SeuratObj) <-"RNA"
  SeuratObj <- CellCycleScoring(SeuratObj, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
  SeuratObj.ST <- SCTransform(SeuratObj, assay = "RNA",vars.to.regress =c( "percent.mt", "S.Score", "G2M.Score" ), return.only.var.genes = FALSE)
  SeuratObj.ST <- FindVariableFeatures(SeuratObj.ST) 
  SeuratObj.ST <- ScaleData(SeuratObj.ST)
  SeuratObj.ST <- RunPCA(SeuratObj.ST)
  
  DefaultAssay(SeuratObj.ST) <- 'ADT'
  # we will use all ADT features for dimensional reduction
  # we set a dimensional reduction name to avoid overwriting the 
  VariableFeatures(SeuratObj.ST) <- rownames(SeuratObj.ST[["ADT"]])
  SeuratObj.ST <- NormalizeData(SeuratObj.ST, assay = "ADT",normalization.method = 'CLR', margin = 2) %>% 
    ScaleData() %>% RunPCA(reduction.name = 'apca')
  return(SeuratObj.ST)
}

#==============================================
#WNN multidimensional clustering and workflow - requires prior computed pca
#===================================================

Seurat.wnn.CITEseq <- function(SeuratObj) {
  
 SeuratObj <- FindMultiModalNeighbors(
    SeuratObj, reduction.list = list("pca", "apca"), 
    dims.list = list(1:30, 1:18), modality.weight.name = "RNA.weight"
  )
  
 SeuratObj <- RunUMAP(SeuratObj, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
 SeuratObj <- FindClusters(SeuratObj, graph.name = "wsnn", algorithm = 3, resolution = 2, verbose = FALSE)
 return(SeuratObj)
 
}


#==============================================
#call SingleR on clusters to assign cell types
#===================================================

Seurat.Singler <- function(SeuratObj, wmm) {
  
  
  df.encode <- readRDS("C:/Users/aflorescu/Molecular Partners AG/DEV_TP_ExVivo - Ana/ToTransfer/encode-personalized-celltypes.rds")
  tt <-SeuratObj@assays[["SCT"]]@data
  singler1<-SingleR(tt, ref=df.encode, labels = df.encode$cell.type,  method = c("cluster"),
                    clusters=SeuratObj@meta.data$seurat_cluster,genes = "de", quantile = 0.8, fine.tune = TRUE,
                    tune.thresh = 0.05, sd.thresh = 1, prune = TRUE,
                    check.missing = TRUE)
  
 
  ClusterCellTypes <- data.frame(singler1@listData[["pruned.labels"]])
  Clusters <- data.frame(seurat_cluster=rownames(ClusterCellTypes), cell.type=ClusterCellTypes[1])
  colnames(Clusters)[2]<-"cell.type"
  rm(tt)
  
  
  
  MetaDataM <- data.frame(CellID=rownames(SeuratObj@meta.data),SeuratObj@meta.data)
  #colnames(MetaDataM)[29]<-"cluster_subtype"
  MetaDataM$cluster.type.singler <- "unclassified"
  for (i in 1: nrow(MetaDataM)){
    for (k in 1: nrow(Clusters)){
      if (MetaDataM$seurat_cluster[i]==Clusters$seurat_cluster[k]) MetaDataM$cluster.type.singler[i]<-Clusters$cell.type[k]
      
      
    } 
  }
  SeuratObj <- AddMetaData(SeuratObj, MetaDataM)

  
  return(SeuratObj)
  
}


#==============================================
#conevrt a Seurat multimodal object to sce. Multimodal includes ADTs and HTOs
#===================================================

Seurat.multimodal.sce <- function(SeuratObj) {

  DefaultAssay(SeuratObj)<-"RNA" 
  sceObj <- as.SingleCellExperiment(SeuratObj)
  DefaultAssay(SeuratObj)<-"ADT" 
  ADT.sce <- as.SingleCellExperiment(SeuratObj)
  altExp(sceObj, "ADT")<-ADT.sce
  #altExpNames(sceObj)<-"ADT"
  DefaultAssay(SeuratObj)<-"HTO" 
  HTO.sce <- as.SingleCellExperiment(SeuratObj)
  altExp(sceObj,"HTO")<-HTO.sce
  #altExpNames(sceObj, e=2)<-"HTO"
  
  return(sceObj)
}