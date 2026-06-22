# NLRC5 Proteomics Heatmaps

This repository contains R scripts used to generate pathway-level proteomics heatmaps from proteomic profiling of B16 melanoma cells expressing different NLRC5 constructs.

## Project overview

The repository includes scripts for generating publication-quality heatmaps related to:

1) Autophagy
2) Proteasome and immunoproteasome
3) ER protein processing
4) Lysosome
5) Melanosome

Heatmaps were generated using the `ComplexHeatmap` R package following log2 transformation and row-wise Z-score normalization.

## Input data

The input dataset consists of protein abundance values obtained from quantitative proteomics experiments.

Expected sample groups:

1) V (Vector control)
2) FL (Full-length NLRC5)
3) SA (NLRC5 variant)

## Repository structure

```text
scripts/      R analysis scripts
data/         Processed proteomics matrices
outputs/      Generated heatmaps and figures
```

## Main R packages

1) ComplexHeatmap
2) circlize
3) dplyr
4) tidyr
5) readr
6) KEGGREST
7) svglite

## Citation

If you use these scripts, please cite:

Shukla A, Moradzad M, et al.

NLRC5 expression within tumor cells is critical to activate adaptive and innate antitumor immune responses.

Immunotherapy Advances.

## Author

Mohammad Moradzad
