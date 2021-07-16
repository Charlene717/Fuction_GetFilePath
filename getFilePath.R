#********************************************************************************************************************#

##################  Function setting ################## 

## Call function
filePath <- ""
#匯入 同一個資料夾中的R檔案
getFilePath <- function(fileName) {
  # path <- setwd("~")  #專案資料夾絕對路徑
  path <- setwd(getwd()) 
  #字串合併無間隔
  # 「<<-」為全域變數給值的指派
  filePath <<- paste0(path ,"/" , fileName)  
  # 載入檔案
  sourceObj <- source(filePath)
  return(sourceObj)
}

#********************************************************************************************************************#

###### Test function ######

## Current path and new folder setting
PathName = setwd(getwd())
RVersion = "20210716V1"
dir.create(paste0(PathName,"/",RVersion))

###### Convert Monocle3 Object to Seurat Object ######
getFilePath("Monocle3_To_Seurat.R")
marrow <- Monocle3_To_Seurat(cds,"cds") #這個function存在於Monocle3_To_Seurat.R裡面

###### Assign Cell-Cycle Scores ######
getFilePath("Cell-Cycle Scoring and Regression.R")
marrow <- CCScorReg(GeneNAFMT,marrow) #這個function存在於Cell-Cycle Scoring and Regression.R裡面

## Color setting
colors_cc <-c("#FF9912B3", "#2e6087", "#417034")  ## Color for Cell-Cycle
RidgePlot(marrow,cols = colors_cc, features = c(Main), ncol = 1)

###### Insert the cell cycle results from Seurat into the  Monocle3 cds object ######
cds@colData@listData$cell_cycle <- marrow@active.ident
# cds@colData@listData$cell_cycle <- marrow@meta.data[["Phase"]]

plot_cells(cds, color_cells_by="cell_cycle", label_cell_groups=FALSE ,show_trajectory_graph = F) + scale_color_manual(values = colors_cc)

## Plot the violin diagram
cds_marrow_cc <- cds[rowData(cds)$gene_short_name %in% Main,]

plot_genes_violin(cds_marrow_cc, group_cells_by="cell_cycle", ncol=2, log_scale = FALSE)+ scale_fill_manual(values = colors_cc)
plot_genes_violin(cds_marrow_cc, group_cells_by="cell_cycle", ncol=2, log_scale = T)+ scale_fill_manual(values = colors_cc)
plot_genes_violin(cds_marrow_cc, group_cells_by="cell_cycle", ncol=2, log_scale = T)+ scale_fill_manual(values = colors_cc)+
  geom_boxplot(width=0.1, fill="white",alpha = 0.7) + theme(axis.text.x=element_text(angle=45, hjust=1))


##############  Annotate your cells according to type (Custom Marker)  ##############

####################    Cell discrimination by AddModuleScore    ####################
getFilePath("Monocle3_AddModuleScore.R")
set.seed(1) # Fix the seed

Marker_PDAC_file_Name <- c("GRUETZMANN_PANCREATIC_CANCER_UP")
Marker_PDAC_Name <- c("PDAC")
cds <- Monocle3_AddModuleScore(Marker_PDAC_file_Name,Marker_PDAC_Name,marrow,cds)
plot_cells(cds, color_cells_by= Marker_PDAC_Name, label_cell_groups=FALSE, show_trajectory_graph = FALSE)
plot_cells(cds, color_cells_by= Marker_PDAC_Name, label_cell_groups=FALSE, show_trajectory_graph = FALSE,cell_size = 1.2) +
  scale_colour_gradient2(low = "#440075", mid = "#ffd261", high = "#4aff8c", 
                         guide = "colourbar",midpoint = 0.2, labs(fill = Marker_PDAC_Name))


