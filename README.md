# NLRC5 Proteomics Heatmaps

This repository contains R scripts used to generate pathway-level proteomics heatmaps from quantitative proteomic profiling of B16 melanoma cells expressing different NLRC5 constructs.


## Project overview

Publication-quality heatmaps were generated for the following biological pathways:

1. Autophagy
2. Proteasome / Immunoproteasome
3. ER protein processing
4. Lysosome
5. Melanosome

Heatmaps shown in Figure 5 and Supplementary Figure S12 of the manuscript were generated using these scripts.


## Experimental groups

The study includes three experimental groups:

- **V** : Vector control
- **N5FL** : Full-length NLRC5
- **N5CA** : NLRC5 SA variant

Heatmaps were generated using group means (average abundance across biological replicates).


## Data processing workflow

1. Import quantitative proteomics matrix.
2. Select pathway-associated proteins.
3. Calculate mean abundance for each experimental group.
4. Perform log2 transformation.
5. Calculate row-wise Z-scores.
6. Generate publication-quality heatmaps using ComplexHeatmap.

---

## Repository structure

```text
scripts/        R scripts
data/           Processed proteomics matrix
outputs/        Generated heatmaps and result files
```

---

## Main R packages

- ComplexHeatmap
- circlize
- dplyr
- tidyr
- readr
- svglite

---

## Citation

If you use these scripts, please cite:

Shukla A, Moradzad M, et al.

NLRC5 expression within tumor cells is critical to activate adaptive and innate antitumor immune responses.


## Author

**Mohammad Moradzad**
