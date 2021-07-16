## Load package
library(Seurat)
library(SummarizedExperiment) 
library(AnnotationDbi)
library(org.Mm.eg.db)
library('org.Hs.eg.db')
library(Hmisc)
library(monocle3)

###### Convert Monocle3 Object to Seurat Object ######

Monocle3_To_Seurat <- function(cds,FigureName) {
CCdata <- cds@assays@data@listData$counts
colnames(CCdata) = cds@assays@data@listData[["counts"]]@Dimnames[[2]]
rownames(CCdata) = cds@assays@data@listData[["counts"]]@Dimnames[[1]]

DataCellcycle <- CCdata
DataCellcycle <- as(as.matrix(DataCellcycle), "dgCMatrix")

marrow <- CreateSeuratObject(counts = DataCellcycle)
marrow <- NormalizeData(marrow)
marrow <- FindVariableFeatures(marrow, selection.method = "vst")
marrow <- ScaleData(marrow, features = rownames(marrow))

marrow@assays[["RNA"]]@counts@Dimnames[[1]] <- cds@assays@data@listData[["counts"]]@Dimnames[[1]]
marrow@assays[["RNA"]]@data@Dimnames[[1]] <- cds@assays@data@listData[["counts"]]@Dimnames[[1]]
marrow@commands[["ScaleData.RNA"]]@params[["features"]] <- cds@assays@data@listData[["counts"]]@Dimnames[[1]]

## Run PCA
marrow <- RunPCA(marrow, features = VariableFeatures(marrow), ndims.print = 1:10, nfeatures.print = 10)
DimHeatmap(marrow, dims = c(1, 2))

png(paste0(PathName,"/",RVersion,"/",RVersion,"_",FigureName,"_CC_DimHP1to10.png")) # 設定輸出圖檔
DimHeatmap(marrow, dims = c(1: 10))
dev.off() # 關閉輸出圖檔

png(paste0(PathName,"/",RVersion,"/",RVersion,"_",FigureName,"_CC_DimHP1&2.png")) # 設定輸出圖檔
DimHeatmap(marrow, dims = c(1, 2))
dev.off() # 關閉輸出圖檔

png(paste0(PathName,"/",RVersion,"/",RVersion,"_",FigureName,"_CC_DimHP3&4.png")) # 設定輸出圖檔
DimHeatmap(marrow, dims = c(3, 4))
dev.off() # 關閉輸出圖檔

png(paste0(PathName,"/",RVersion,"/",RVersion,"_",FigureName,"_CC_DimHP5&6.png")) # 設定輸出圖檔
DimHeatmap(marrow, dims = c(5, 6))
dev.off() # 關閉輸出圖檔

png(paste0(PathName,"/",RVersion,"/",RVersion,"_",FigureName,"_CC_DimHP7&8.png")) # 設定輸出圖檔
DimHeatmap(marrow, dims = c(7, 8))
dev.off() # 關閉輸出圖檔

png(paste0(PathName,"/",RVersion,"/",RVersion,"_",FigureName,"_CC_DimHP9&10.png")) # 設定輸出圖檔
DimHeatmap(marrow, dims = c(9, 10))
dev.off() # 關閉輸出圖檔

return(marrow)
}

