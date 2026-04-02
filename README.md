# README

## Citation
Savage N., Grewal S., Shaikh V., Zemp F., McKenna D., Mikolajewicz N., Najem H., Pyczek J., Wei J., Taleb M., Asselstine L., Anand A., Chafe S., Zhai K., Maich W., Chokshi C., Patel H., Korman T., Subapanditha M., Shapovalova Z., Tatari N., Miltec P., Chen D., Pacheco S., Han H., Chan J., Brown K., Venugopal C., Kislinger T., Heimberger A., Moffat J., Mahoney D., Singh S. (2026). 	Dual tumor-myeloid targeting of glioblastoma with GPNMB CAR-T cells. *Nature*. 

## Data Availability
- Myeloid Meta-Atlas: [https://doi.org/10.6084/m9.figshare.31926300](https://doi.org/10.6084/m9.figshare.31926300)
- GL261 WT vs KO tumor dataset: [https://doi.org/10.6084/m9.figshare.31926462](https://doi.org/10.6084/m9.figshare.31926462)

## Files
- `01_export_gpnmb_wtko_publishable.Rmd`: exports the GL261 WT vs KO tumor dataset, and reproduces tumor-side UMAP, tumor-state, DEG, and Neftel GSEA outputs.
- `02_export_myeloid_meta_atlas_publishable.Rmd`: exports the CNS myeloid meta-atlas, and reproduces atlas UMAP and GPNMB/state summary outputs.
- `03_gpnmb_scRNA_methods_reproducibility.Rmd`: methods-oriented notebook that consolidates figure and statistics reproduction for the tumor and myeloid analyses.
- `04_gpnmb_myeloid_meta_atlas_nmf_reproducibility.R`: myeloid atlas NMF workflow. 

## Software requirements

- R with `rmarkdown` and `knitr` for rendering the notebooks
- Core packages used across the release scripts: `Seurat`, `Matrix`, `dplyr`, `tidyr`, `tibble`, `readr`, `stringr`, `ggplot2`, `patchwork`
- Additional packages used in specific workflows: `fgsea`, `foreach`, `doParallel`, `readxl`
- Optional acceleration for `04`: `glmGamPoi` if you enable atlas rebuild mode and want faster `SCTransform`
