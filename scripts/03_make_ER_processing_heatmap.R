# ER protein processing heatmap
# Input: protein abundance matrix with Genes, V1-4, N5FL1-4, N5CA1-4
# Output: group-mean heatmap as PDF and SVG

required_packages <- c(
  "readr", "dplyr", "tidyr", "tibble", "stringr",
  "ComplexHeatmap", "circlize", "grid", "svglite"
)

missing_packages <- required_packages[
  !required_packages %in% rownames(installed.packages())
]

if (length(missing_packages) > 0) {
  install.packages(missing_packages)
}

library(readr)
library(dplyr)
library(tidyr)
library(tibble)
library(stringr)
library(ComplexHeatmap)
library(circlize)
library(grid)
library(svglite)

input_file <- "data/processed/data/processed/protein_abundance_matrix_B16_V_N5FL_N5CA.tsv"
output_dir <- "outputs/ER_processing"

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

data <- read_tsv(
  input_file,
  na = c("", "NA", "NaN"),
  show_col_types = FALSE,
  progress = FALSE
)

samples <- c(
  "V1", "V2", "V3", "V4",
  "N5FL1", "N5FL2", "N5FL3", "N5FL4",
  "N5CA1", "N5CA2", "N5CA3", "N5CA4"
)

required_cols <- c("Genes", samples)
missing_cols <- setdiff(required_cols, colnames(data))

if (length(missing_cols) > 0) {
  stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
}

data <- data %>%
  select(Genes, all_of(samples)) %>%
  mutate(across(-Genes, as.numeric))

target_genes <- c(
  "Fbxo2", "Eif2s1", "Hspbp1", "Cryab", "Ckap4", "Ganab", "Sec23a", "Mapk9",
  "Capn1", "Hsph1", "Hspa4l", "Syvn1", "Derl1", "Ube2g1", "Edem3", "Atxn3",
  "Ube4b", "Sec24a", "Sec24c", "Man1a2", "Sar1b", "Capn2", "Sil1", "Dnaja1",
  "Ero1a", "Sec61a1", "Sec23b", "Stt3b", "Sel1l", "Uggt1", "Canx", "Stt3a",
  "Fbxo6", "Rpn1", "Dad1", "Krtcap2", "Bak1", "Lman2", "Bag1", "Uggt2",
  "Sec13", "Skp1", "Rad23b", "Wfs1", "Sec61a2", "Erp29", "Hyou1", "Erlec1",
  "Bcl2", "Sec31a", "Sec24b", "Rrbp1", "Ostc", "Sec63", "Ssr4", "Amfr",
  "Mogs", "Ddost"
)

annotation <- tibble(
  Genes = target_genes,
  Component = case_when(
    Genes %in% c(
      "Sec23a", "Sec23b", "Sec24a", "Sec24b", "Sec24c",
      "Sec31a", "Sec13", "Sar1b"
    ) ~ "COPII trafficking",

    Genes %in% c(
      "Sec61a1", "Sec61a2", "Sec63", "Ssr4", "Rrbp1"
    ) ~ "ER translocation",

    Genes %in% c(
      "Stt3a", "Stt3b", "Dad1", "Ddost", "Ostc", "Rpn1"
    ) ~ "N-glycosylation / OST",

    Genes %in% c(
      "Mogs", "Ganab", "Man1a2", "Uggt1", "Uggt2",
      "Canx", "Erlec1", "Lman2"
    ) ~ "Folding / QC",

    Genes %in% c(
      "Derl1", "Sel1l", "Syvn1", "Ube2g1", "Ube4b",
      "Atxn3", "Amfr", "Rad23b", "Fbxo2", "Fbxo6", "Edem3"
    ) ~ "ERAD",

    Genes %in% c(
      "Sil1", "Dnaja1", "Hspbp1", "Hsph1", "Hspa4l",
      "Cryab", "Erp29", "Hyou1", "Wfs1", "Ero1a", "Ckap4"
    ) ~ "Chaperones",

    Genes %in% c(
      "Bak1", "Bcl2", "Capn1", "Capn2", "Mapk9",
      "Bag1", "Skp1", "Krtcap2", "Eif2s1"
    ) ~ "Stress / apoptosis",

    TRUE ~ NA_character_
  )
)

df <- data %>%
  separate_rows(Genes, sep = ";") %>%
  mutate(Genes = str_trim(Genes)) %>%
  filter(Genes %in% target_genes) %>%
  group_by(Genes) %>%
  summarise(across(all_of(samples), ~ mean(.x, na.rm = TRUE)), .groups = "drop") %>%
  left_join(annotation, by = "Genes")

component_levels <- c(
  "COPII trafficking",
  "ER translocation",
  "N-glycosylation / OST",
  "Folding / QC",
  "ERAD",
  "Chaperones",
  "Stress / apoptosis"
)

df <- df %>%
  mutate(Component = factor(Component, levels = component_levels)) %>%
  arrange(Component, Genes)

df_mean <- df %>%
  mutate(
    V = rowMeans(across(c(V1, V2, V3, V4)), na.rm = TRUE),
    N5FL = rowMeans(across(c(N5FL1, N5FL2, N5FL3, N5FL4)), na.rm = TRUE),
    N5CA = rowMeans(across(c(N5CA1, N5CA2, N5CA3, N5CA4)), na.rm = TRUE)
  ) %>%
  select(Genes, Component, V, N5FL, N5CA)

mat <- df_mean %>%
  select(Genes, V, N5FL, N5CA) %>%
  column_to_rownames("Genes") %>%
  as.matrix()

mat <- log2(mat + 1)

row_zscore <- function(x) {
  z <- t(scale(t(x)))
  z[is.na(z)] <- 0
  z
}

mat_z <- row_zscore(mat)

row_components <- df_mean$Component
names(row_components) <- df_mean$Genes
row_components <- row_components[rownames(mat_z)]

component_colors <- c(
  "COPII trafficking" = "#1f78b4",
  "ER translocation" = "#33a02c",
  "N-glycosylation / OST" = "#6a3d9a",
  "Folding / QC" = "#ff7f00",
  "ERAD" = "#e31a1c",
  "Chaperones" = "#a6cee3",
  "Stress / apoptosis" = "#b15928"
)

row_anno <- rowAnnotation(
  Component = row_components,
  col = list(Component = component_colors),
  show_annotation_name = FALSE
)

col_fun <- colorRamp2(
  c(-2, 0, 2),
  c("#2166AC", "white", "#B2182B")
)

ht <- Heatmap(
  mat_z,
  name = "Row Z-score",
  col = col_fun,
  left_annotation = row_anno,
  row_split = row_components,
  cluster_rows = TRUE,
  cluster_row_slices = FALSE,
  cluster_columns = FALSE,
  row_names_side = "left",
  column_labels = c("V", "N5FL", "N5CA"),
  row_names_gp = gpar(fontsize = 8),
  column_names_gp = gpar(fontsize = 11, fontface = "bold"),
  column_title = NULL,
  row_title = NULL
)

pdf(
  file.path(output_dir, "ER_processing_heatmap_group_means.pdf"),
  width = 7,
  height = 10
)
draw(ht)
dev.off()

svglite(
  file.path(output_dir, "ER_processing_heatmap_group_means.svg"),
  width = 7,
  height = 10
)
draw(ht)
dev.off()

write_tsv(
  df_mean,
  file.path(output_dir, "ER_processing_genes_used_for_heatmap.tsv")
)

sessionInfo()
