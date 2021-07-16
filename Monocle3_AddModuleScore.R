# https://rdrr.io/github/satijalab/seurat/man/AddModuleScore.html
# https://academic.oup.com/nar/article/47/21/e133/5531181#185740462

# Marker_file_Name <- c("GRUETZMANN_PANCREATIC_CANCER_UP")
# Marker_Name <- c("PDAC_Marker")

Monocle3_AddModuleScore <- function(Marker_file_Name,Marker_Name,marrow,cds) {
  PathName = setwd(getwd())
  #PDAC_Marker
    PDAC_Marker_file <- paste0(PathName,"/",Marker_file_Name,".txt")
  
  PDAC_Marker <- read.delim(PDAC_Marker_file, header = TRUE, stringsAsFactors = FALSE)
  PDAC_Marker <- as.data.frame(PDAC_Marker[-1,])
  
  marrow_PDAC_Marker <- AddModuleScore(
    object = marrow,
    features = PDAC_Marker,
    ctrl = 5,
    name = Marker_Name)
  
  cds@colData@listData[[paste0(Marker_Name)]] <- marrow_PDAC_Marker@meta.data[[paste0(Marker_Name,"1")]]
  
  return(cds)
}


# # Test
#  Marker_file_Name <- c("GRUETZMANN_PANCREATIC_CANCER_UP")
#  Marker_Name <- c("PDAC_Marker")
#  cds <- Monocle3_AddModuleScore(Marker_file_Name,Marker_Name,marrow,cds)
#  plot_cells(cds, color_cells_by= Marker_Name, label_cell_groups=FALSE, show_trajectory_graph = FALSE)
#  
#  cds_sub_DucT2 <- Monocle3_AddModuleScore(Marker_file_Name,Marker_Name,marrow_sub_DucT2,cds_sub_DucT2)
#  plot_cells(cds_sub_DucT2, color_cells_by= Marker_Name, label_cell_groups=FALSE, show_trajectory_graph = FALSE)
