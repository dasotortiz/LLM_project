library(Seurat)
library(viridis)
library(scCustomize)
library(ggplot2)
library(patchwork)
library(dplyr)

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

# populations barplots

hca_ref$ref.cell.type <- hca_ref$CellType
hca_ref$annotation_2 <- as.character(hca_ref$annotation_2)
hca_ref$ref.cell.type <- ifelse(!is.na(hca_ref$annotation_2), hca_ref$annotation_2, hca_ref$ref.cell.type)

xpand_counts <- FetchData(xpand, vars = "annotation_2") %>%
  group_by(annotation_2) %>%
  summarise(Freq = n()) %>%
  mutate(dataset = "Xpand", cell_type = annotation_2) %>%
  select(-annotation_2
  
  )
hca_counts <- FetchData(hca_ref, vars = "ref.cell.type") %>%
  group_by(ref.cell.type) %>%
  summarise(Freq = n()) %>%
  mutate(dataset = "HCA", cell_type = ref.cell.type) %>%
  select(-ref.cell.type)

cell_counts <- rbind(xpand_counts, hca_counts)
cell_type_levels <- c(celltype_order, setdiff(unique(cell_counts$cell_type), celltype_order))
cell_counts$cell_type <- factor(cell_counts$cell_type, levels = cell_type_levels)

pdf("/ibex/user/sotoorda/LLM_project/Additional_visualization/celltype_counts.pdf", width = 9, height = 5)
ggplot(cell_counts, aes(x = cell_type, y = Freq, fill = dataset)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  labs(x = "Cell type", y = "Number of cells", fill = "Dataset") +
  scale_fill_manual(values = c("HCA" = "steelblue", "Xpand" = "red")) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(hjust = 0.5)
  ) +
  ggtitle("Cell counts by dataset")
dev.off()

