# load libraries
library(Seurat)
library(Signac)
library(Matrix)

output <- "/ibex/user/sotoorda/masterh1/LLM project/CellTypist/Objects"

xpand <- readRDS("/ibex/user/sotoorda/masterh1/public_data_analysis/3_datasets_analysis/objects/individual_datasets/final_version/xpand_final_annot.rds")
hca_ref <- readRDS('/ibex/user/sotoorda/masterh1/public_data_analysis/hca/objects/hca_OriData_integrated_Marcal_analysis.rds')

# add ref.cell.type column to HCA
hca_ref$ref.cell.type <- hca_ref$CellType
hca_ref$annotation_2 <- as.character(hca_ref$annotation_2)
hca_ref$ref.cell.type <- ifelse(!is.na(hca_ref$annotation_2), hca_ref$annotation_2, hca_ref$ref.cell.type)
hca_ref$annotation_2 <- hca_ref$ref.cell.type

objects_list <- list(xpand = xpand, hca_ref = hca_ref)

for(name in names(objects_list)) {
    obj <- objects_list[[name]]
    # Get expression data
    counts_matrix <- GetAssayData(obj, assay = "RNA", slot = "counts")
    writeMM(counts_matrix, file = paste0(output, '/', name, '_matrix.mtx'))
    # write gene names
    write.table(data.frame('gene'=rownames(counts_matrix)),
            file=paste0(output, '/', name, '_gene_names.csv', sep=''),
            quote=F,row.names=F,col.names=F)

    # save metadata table:
    dataframe <- data.frame('barcode' = colnames(obj),
                        'CellType' = obj$annotation_2)

    write.csv(dataframe, file=paste0(output, '/', name, '_metadata.csv'), 
          quote=F, row.names=F)
}

