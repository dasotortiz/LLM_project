# LLM Project

This repository collects a small set of notebook-driven experiments for cell type annotation and benchmarking on single-cell data. It compares multiple annotation approaches, including scGPT, CellTypist, Random Forest, and Seurat label transfer, then combines their outputs into consensus and benchmark summaries.

## What's in the repo

- [LLM_model/LLM.ipynb](LLM_model/LLM.ipynb) builds a consensus annotation workflow from CellTypist, Random Forest, and Seurat predictions.
- [scGPT/ZeroShot_Reference_Mapping.ipynb](scGPT/ZeroShot_Reference_Mapping.ipynb) runs a zero-shot reference mapping workflow with scGPT.
- [benchmarking/benchmarking.ipynb](benchmarking/benchmarking.ipynb) compares model predictions and calculates benchmark metrics.
- [Label_transfer/Seurat_Label_Transfer.R](Label_transfer/Seurat_Label_Transfer.R) contains the Seurat label transfer workflow.
- The CSV and JSON files in the subfolders are generated outputs from those workflows.

## Repository layout

- [Additional_visualization/](Additional_visualization/) - plots and visual checks
- [CellTypist/](CellTypist/) - CellTypist prediction outputs and related artifacts
- [LLM_model/](LLM_model/) - consensus notebook and results
- [Label_transfer/](Label_transfer/) - Seurat label transfer code and outputs
- [Random_forest/](Random_forest/) - random forest prediction outputs
- [benchmarking/](benchmarking/) - benchmark notebook and summary files
- [scGPT/](scGPT/) - scGPT notebook, model assets, and prediction outputs

## Outputs

Common generated files include:

- consensus result tables and metrics from the LLM notebook
- scGPT prediction probabilities and labels
- benchmark comparison CSV and JSON summaries
- Seurat label transfer predictions
