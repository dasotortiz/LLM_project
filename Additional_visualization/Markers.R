library(Seurat)
library(viridis)
library(scCustomize)
library(ggplot2)
library(patchwork)

# Define the markers for each cell population
markers <- c(
  #MEP
 'GATA2','FCER1A', 'PBX1',
  # MDP
  'IRF8','IRF7','SPIB',
  # GMP
  'ELANE', 'CTSG', 'PRTN3',
  # Eo/B/Mast
 'CSF2RB', 'MS4A2', 'MS4A3', 'CLC',
  # ERP
  'GATA1', 'CA1', 'APOC1', 'HBB',
  # PreB
  'CD24', 'SMC4', 'VPREB1', 'TOP2A', 'IGHM',
  # ProB
  'RAG1', 'AKAP12', 'ARPP21', 'TOP2B',
  #CLP
  'ADA', 'JCHAIN', 'DNTT',
  #MPP
  'SPINK2', 'SMIM24',
  # HSC
  'AVP', 'CD164', 'CRHBP'
)

# Load and preprocess the datasets
xpand <- readRDS("/ibex/user/sotoorda/masterh1/public_data_analysis/3_datasets_analysis/objects/individual_datasets/final_version/xpand_final_annot.rds")
xpand <- NormalizeData(xpand)
celltypes <- unique(xpand$annotation_2)

hca_ref <- readRDS('/ibex/user/sotoorda/masterh1/public_data_analysis/hca/objects/hca_OriData_integrated_Marcal_analysis.rds')
hca_subset <- subset(hca_ref, annotation_2 %in% celltypes)

# Ensure consistent factor order across both datasets
celltype_order <- c("HSC", "MPP", "CLP", "ProB", "PreB", "ERP", "Eo/B/Mast", "GMP", "MDP", 'MEP')
xpand$annotation_2 <- factor(xpand$annotation_2, levels = celltype_order)
hca_subset$annotation_2 <- factor(hca_subset$annotation_2, levels = celltype_order)

# Generate the dot plots for both datasets
pal <- viridis(n = 10, option = "D")
Idents(xpand) <- "annotation_2"
Idents(hca_subset) <- "annotation_2"
pdf("/ibex/user/sotoorda/LLM_project/Additional_visualization/markers_expression.pdf",width = 8, height = 8)
DotPlot_scCustom(seurat_object = xpand, features = markers, x_lab_rotate = TRUE, colors = pal) +
  coord_flip() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    plot.title = element_text(hjust = 0.5)
  ) +
  ggtitle("Xpand")
DotPlot_scCustom(seurat_object = hca_subset, features = markers, x_lab_rotate = TRUE, colors = pal) +
  coord_flip() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    plot.title = element_text(hjust = 0.5)
  ) +
  ggtitle("HCA")
dev.off()
