library(Seurat)
library(ggrepel)

seurat <- readRDS("/ibex/user/sotoorda/masterh1/public_data_analysis/3_datasets_analysis/objects/individual_datasets/final_version/xpand_final_annot.rds")

# Remove redundant raw counts
seurat[["RNA"]]$countsraw_counts <- NULL

seurat[["RNA"]] <- split(seurat[["RNA"]], f = seurat$sample)

# Preprocessing functions
npcs<-30
seurat <- NormalizeData(seurat)
seurat <- FindVariableFeatures(seurat)
seurat <- ScaleData(seurat)
seurat <- RunPCA(seurat, npcs = npcs)
seurat <-RunUMAP(object = seurat, dims = 1:npcs, reduction = "pca", reduction.name = "umap_unintegrated", seed.use = 123)

pdf("/ibex/user/sotoorda/masterh1/public_data_analysis/hca/plots/OriginalData/unintegrated.pdf")
DimPlot(seurat, reduction = "umap_unintegrated", group.by = "orig.ident")
dev.off()


# Integration
options(future.globals.maxSize = 16 * 1024^3) # Increase memory allocation for the integration
# seurat_int <- IntegrateLayers(object = seurat, method = CCAIntegration, orig.reduction = "pca", new.reduction = "integrated.cca", verbose = FALSE, assay = "RNA")
seurat_int <- IntegrateLayers(object = seurat, method = HarmonyIntegration, orig.reduction = "pca", new.reduction = "harmony", verbose = FALSE, assay = "RNA")
# seurat_int <- IntegrateLayers(object = seurat_int, method = RPCAIntegration, orig.reduction = "pca", new.reduction = "integrated.rpca", verbose = FALSE, assay = "RNA")

# re-join layers after integration
seurat_int[["RNA"]] <- JoinLayers(seurat_int[["RNA"]])


# UMAP
npcs<-30
seurat_int <- RunUMAP(object = seurat_int, reduction = "integrated.cca", dims = 1:npcs, reduction.name = "umap_integrated_cca", seed.use = 123)
seurat_int <- RunUMAP(object = seurat_int, reduction = "harmony", dims = 1:npcs, reduction.name = "umap_integrated_harmony", seed.use = 123)
seurat_int <- RunUMAP(object = seurat_int, reduction = "integrated.rpca", dims = 1:npcs, reduction.name = "umap_integrated_rpca", seed.use = 123)


# Adding predictions as metadata
cm_celltypist <- read.table('/ibex/user/sotoorda/LLM_project/CellTypist/predictions_celltypist.csv', header = TRUE, sep = ',')
cm_rf <- read.table('/ibex/user/sotoorda/LLM_project/Random_forest/prediction_results.csv', header = TRUE, sep = ',')
cm_seurat <- read.table('/ibex/user/sotoorda/LLM_project/Label_transfer/predictions_seurat.csv', header = TRUE, sep = ',')
cm_scgpt <- read.table('/ibex/user/sotoorda/LLM_project/scGPT/scgpt_prob_with_labels.csv', header = TRUE, sep = ',')
colnames(cm_scgpt)[colnames(cm_scgpt) == "scGPT_prediction"] <- "Predicted_label"

rownames(cm_celltypist) <- cm_celltypist$Cell_id
rownames(cm_rf) <- cm_rf$Cell_id
rownames(cm_seurat) <- cm_seurat$Cell_id
rownames(cm_scgpt) <- cm_scgpt$Cell_id

intersect(rownames(cm_celltypist), rownames(seurat_int@meta.data))
rownames(seurat@meta.data)
seurat_int$pred_celltypist <- cm_celltypist$Predicted_label[match(rownames(seurat_int@meta.data), cm_celltypist$Cell_id)]
seurat_int$pred_rf <- cm_rf$Predicted_label[match(rownames(seurat_int@meta.data), cm_rf$Cell_id)]
seurat_int$pred_seurat <- cm_seurat$Predicted_label[match(rownames(seurat_int@meta.data), cm_seurat$Cell_id)]
seurat_int$pred_scgpt <- cm_scgpt$Predicted_label[match(rownames(seurat_int@meta.data), cm_scgpt$Cell_id)]

celltype_colors <- c(
  "HSC" = "steelblue",
  "Stromal" = "black",

  # Stem/progenitors
  "MPP" = "#9AC0CD",
  "CLP" = "#98FB98",

  # T/NK lineage (greens)
  "NK cells" = "#88EB88",
  "pre-T" = "#79DC79",
  "Naive T-cell" = "#6ACD6A",
  "CD8 T-cell" = "#5BBE5B",

  # B lineage (greens, darker)
  "ProB" = "#4CAF4C",
  "PreB" = "#3CA03C",
  "Follicular B cell" = "#1E821E",
  "pre-PC" = "#0F730F",
  "Plasma Cell" = "#006400",

  # Myeloid progenitors
  "GMP" = "#F5CD94",

  # Erythroid / MEP / MKP (purples)
  "MEP" = "#FFE1FF",
  "MKP" = "#E2BFEB",
  "ERP" = "#AA7DC4",
  "Early-Erythroblast" = "#8D5CB1",
  "Erythroblast" = "#713B9E",
  "Platelet" = "#551A8B",

  # Eo/B/Mast & granulocytic (oranges)
  "Eo/B/Mast" = "#FFE7BA",
  "Granulocytic-UNK" = "#EBB36F",
  "Immature-Neutrophil" = "#E1994A",
  "Neutrophil" = "#D77F25",
  "Eosinophil" = "#CD6600",

  # MDP / dendritic / monocyte (blues)
  "MDP" = "#5C73AE",
  "Pre-Dendritic" = "#3D4C9E",
  "Dendritic Cell" = "#1E268F",
  "Monocyte" = "#000080"
)
unique(seurat_int$annotation_2)
# UMAPS for integrated data
pdf("/ibex/user/sotoorda/LLM_project/Additional_visualization/xpand_umap.pdf", width = 15, height = 10)
DimPlot(seurat_int, group.by = "annotation_2", reduction = "umap_integrated_harmony", label = TRUE, cols = celltype_colors, pt.size = 0.4, label.size = 5)
DimPlot(seurat_int, group.by = "pred_celltypist", reduction = "umap_integrated_harmony", label = TRUE, cols = celltype_colors, pt.size = 0.4, label.size = 5)
DimPlot(seurat_int, group.by = "pred_rf", reduction = "umap_integrated_harmony", label = TRUE, cols = celltype_colors, pt.size = 0.4, label.size = 5)
DimPlot(seurat_int, group.by = "pred_seurat", reduction = "umap_integrated_harmony", label = TRUE, cols = celltype_colors, pt.size = 0.4, label.size = 5)
DimPlot(seurat_int, group.by = "pred_scgpt", reduction = "umap_integrated_harmony", label = TRUE, cols = celltype_colors, pt.size = 0.4, label.size = 5)
dev.off()
