library(Seurat)

# Cell population markers
HSC_markers <- c("CD34", "CD164", "BEX1", "BEX2", "AVP", "CRHBP", "HLF")
MPP_markers <- c("CD34", "CD33", "MPO", "FLT3", "MZB1")
GMP_markers <- c("MPO", "AZU1", "PRTN3", "ELANE", "LYZ")
MDP_markers <- c("MPO", "LYZ", "IRF8", "LY86", "RUNX2", "LILRB4")
MEP_markers <- c("GATA2", "FCER1A", "ITGA2B", "CSF2RB")
EP_markers <- c("GATA1", "EPOR", "CA1", "CA2", "EPCAM", "KLF1", "BLVRB", "APOC1", "APOE")
MKP_markers <- c("GATA2", "FCER1A", "ITGA2B", "PLEK", "PPBP", "PF4", "GP9")
BMP_markers <- c("MS4A2", "MS4A3", "TPSAB1", "TPSB2", "HDC", "CLC", "PRG2")
CLP_markers <- c("FLT3", "TRBC2", "MZB1", "LTB", "JCHAIN", "ADA", "BCL2")
Pro_B_markers <- c("DNTT", "MME", "PAX5", "RAG1", "RAG2")
Pre_B_markers <- c("VPREB1", "IGLL1", "JCHAIN", "TCL1A")
# Nuria's markers
CLP_nuria <- c('FLT3', 'UMODL1', 'LTB', 'CYGB', 'RGL4', 'ADA', 'IL7R', 'TESC')
pro_b <- c('MME', 'ARPP21', 'TOP2B', 'VPREB3', 'VPREB1', 'BTG2', 'CD79B')
cycling_pre_b <- c('HMGB1', 'HMGB2', 'EBF1', 'XRCC5', 'PRKDC')
Pre_b <- c('IGHM', 'CD74', 'RPS27', 'IGHD')
all_markers <- c(HSC_markers, MPP_markers, MEP_markers, EP_markers, MKP_markers, CLP_markers,
Pro_B_markers, Pre_B_markers, MDP_markers, GMP_markers, BMP_markers, CLP_nuria, pro_b, cycling_pre_b, Pre_b)
unique_markers <- unique(all_markers)

xpand <- readRDS("/ibex/user/sotoorda/masterh1/public_data_analysis/3_datasets_analysis/objects/individual_datasets/final_version/xpand_final_annot.rds")
hca_ref <- readRDS('/ibex/user/sotoorda/masterh1/public_data_analysis/hca/objects/hca_OriData_integrated_Marcal_analysis.rds')

# Getting the variable features for both datasets
xpand <- NormalizeData(xpand)
xpand <- FindVariableFeatures(xpand, nfeatures = 2000)
vf_xpand <- VariableFeatures(xpand)
hca_ref <- FindVariableFeatures(hca_ref, nfeatures = 2000)
vf_hca_ref <- VariableFeatures(hca_ref)

# Taking the union of variable features from both datasets
features_2_use <- intersect(vf_xpand, vf_hca_ref)
features_2_use <- unique(c(features_2_use, unique_markers))

hca_ref$ref.cell.type <- hca_ref$CellType
hca_ref$annotation_2 <- as.character(hca_ref$annotation_2)
hca_ref$ref.cell.type <- ifelse(!is.na(hca_ref$annotation_2), hca_ref$annotation_2, hca_ref$ref.cell.type)

# Extracting the expression data for the features Xpand
xpand_subset <- subset(xpand, features = features_2_use)
data <- as.data.frame(t(GetAssayData(xpand_subset, layer = "data")))
data$CellType <- xpand$annotation_2
data$CellName <- rownames(data)
colnames(data)

# Extracting the expression data for the features HCA reference
hca_ref_subset <- subset(hca_ref, features = features_2_use)
data_ref <- as.data.frame(t(GetAssayData(hca_ref_subset, layer = "data")))
data_ref$CellType <- hca_ref$ref.cell.type
data_ref$CellName <- rownames(data_ref)
colnames(data_ref)

# Ensuring that both datasets have the same features (genes)
features_2_keep <- intersect(colnames(data), colnames(data_ref))
data <- data[, features_2_keep]
data_ref <- data_ref[, features_2_keep]

# Exporting expression data to Python
write.csv(data, file = "/ibex/user/sotoorda/masterh1/LLM project/Random_forest/test_cell_data.csv", row.names = FALSE)
write.csv(data_ref, file = "/ibex/user/sotoorda/masterh1/LLM project/Random_forest/ref_cell_data.csv", row.names = FALSE)




