library(Seurat)

# Load and normalize the data
xpand <- readRDS("/ibex/user/sotoorda/masterh1/public_data_analysis/3_datasets_analysis/objects/individual_datasets/final_version/xpand_final_annot.rds")
xpand <- NormalizeData(xpand)
hca_ref <- readRDS('/ibex/user/sotoorda/masterh1/public_data_analysis/hca/objects/hca_OriData_integrated_Marcal_analysis.rds')

# Homogenize labels in the reference dataset
hca_ref$ref.cell.type <- hca_ref$CellType
hca_ref$annotation_2 <- as.character(hca_ref$annotation_2)
hca_ref$ref.cell.type <- ifelse(!is.na(hca_ref$annotation_2), hca_ref$annotation_2, hca_ref$ref.cell.type)

# Classify query cells using the reference dataset
anchors <- FindTransferAnchors(reference = hca_ref, query = xpand, dims = 1:30, reference.reduction = "pca")
predictions <- TransferData(anchorset = anchors, refdata = hca_ref$ref.cell.type, dims = 1:30)

# Get the probability matrix
predictions$Cell_id <- rownames(predictions)
# Add the ground truth labels to the predictions data frame
predictions$Ground_truth <- xpand$annotation_2[match(predictions$Cell_id, rownames(xpand@meta.data))]
colnames(predictions)[1] <- 'Predicted_label'

# Save the matrix
write.csv(predictions, file = "/ibex/user/sotoorda/masterh1/LLM project/Label_transfer/predictions_seurat.csv", row.names = FALSE)
