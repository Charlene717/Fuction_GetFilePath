## Load package
library(Seurat)
library(SummarizedExperiment) 
library(AnnotationDbi)
library(org.Mm.eg.db)
library('org.Hs.eg.db')
library(Hmisc)
library(tidyverse)
library(monocle3)

###### Cell-Cycle Scoring and Regression (for Monocle3) ######
# GeneNAFMT <- c("HuGSymbol") # Gene names format  of data: HuGSymbol,MouGSymbol,HuENSEMBL,MouENSEMBL

CCScorReg <- function(GeneNAFMT,marrow) {
## (Version 2)

## Cell cycle genes file
cc.genes_list <- read.csv(paste0(PathName,"/Cell cycle/regev_lab_cell_cycle_genesCh.csv")) # A list of cell cycle markers, from Tirosh et al, 2015, is loaded with Seurat. 
# cc: Cell-Cycle


###### Assign Cell-Cycle Scores ######
cc.genes <- cc.genes_list # cc: Cell-Cycle

## Segregate "cc.genes_list" into markers of G2/M phase and markers of S phase
s.genes <- cc.genes[1:43,]
g2m.genes <- cc.genes[44:97,]

s.genes <- as.character(s.genes)
g2m.genes <- as.character(g2m.genes)

## Assign Cell-Cycle Scores
if (GeneNAFMT=="HuGSymbol") {
  # For Human gene symbol # HuGSymbol
  s.genes <- as.character(s.genes)
  g2m.genes <- as.character(g2m.genes)
  
} else if (GeneNAFMT=="MouGSymbol") {
  # For Mouse gene symbol # MouGSymbol
  s.genes <- capitalize(tolower(s.genes))
  g2m.genes <- capitalize(tolower(g2m.genes))
  
} else if (GeneNAFMT=="HuENSEMBL") {
  # For Human ENSEMBL # HuENSEMBL
  s.genes3_0 <- AnnotationDbi::select(org.Hs.eg.db, keys=s.genes, columns='ENSEMBL', keytype='SYMBOL')
  s.genes3_1 <- s.genes3_0[!duplicated(s.genes3_0[2]),]
  s.genes3_2 <- na.omit(s.genes3_1[2])
  s.genes <- s.genes3_2[,1]
  
  g2m.genes3_0 <- AnnotationDbi::select(org.Hs.eg.db, keys=g2m.genes, columns='ENSEMBL', keytype='SYMBOL')
  g2m.genes3_1 <- g2m.genes3_0[!duplicated(g2m.genes3_0[2]),]
  g2m.genes3_2 <- na.omit(g2m.genes3_1[2])
  g2m.genes <- g2m.genes3_2[,1]
  
  
  
} else if (GeneNAFMT=="MouENSEMBL") {
  # For Mouse ENSEMBL # MouENSEMBL
  s.genes4_0 <- AnnotationDbi::select(org.Mm.eg.db, keys=s.genes2, columns='ENSEMBL', keytype='SYMBOL')
  s.genes4_1 <- s.genes4_0[!duplicated(s.genes4_0[2]),]
  s.genes4_2 <- na.omit(s.genes4_1[2])
  s.genes <- s.genes4_2[,1]
  
  g2m.genes4_0 <- AnnotationDbi::select(org.Mm.eg.db, keys=g2m.genes2, columns='ENSEMBL', keytype='SYMBOL')
  g2m.genes4_1 = g2m.genes4_0[!duplicated(g2m.genes4_0[2]),]
  g2m.genes4_2 <- na.omit(g2m.genes4_1[2])
  g2m.genes <- g2m.genes4_2[,1]
  
  
  
} else {
  c("Wrong input from GeneNAFMT!")
}


marrow <- CellCycleScoring(marrow, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

 ## view cell cycle scores and phase assignments
 head(marrow[[]])

return(marrow)
}





