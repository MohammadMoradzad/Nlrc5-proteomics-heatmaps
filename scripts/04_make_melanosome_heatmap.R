# Melanosome / pigment granule heatmap
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
output_dir <- "outputs/melanosome"

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

melanin_synthesis <- c(
  "Tyr", "Tyrp1", "Dct"
)

melanosome_structure_biogenesis <- c(
  "Pmel", "Mlph", "Mlana", "Oca2", "Gpnmb"
)

melanosome_trafficking <- c(
  "Rab27a", "Rab38", "Rab32", "Myo5a",
  "Rilp", "Ap3b1", "Ap3d1", "Ap3m1", "Ap3s1"
)

ion_transport_pH <- c(
  "Slc45a2", "Slc24a5",
  "Atp6v0d1", "Atp6v1a", "Atp6v1b2",
  "Atp6v1c1", "Atp6v1d"
)

regulatory_melanocyte_markers <- c(
  "Mitf", "Sox10", "Tfap2a", "Kit", "Mc1r"
)

target_genes <- unique(c(
  melanin_synthesis,
  melanosome_structure_biogenesis,
  melanosome_trafficking,
  ion_transport_pH,
  regulatory_melanocyte_markers
))

annotation <- tibble(
  Genes = target_genes,
  Component = case_when(
    Genes %in% melanin_synthesis ~ "Melanin synthesis",
    Genes %in% melanosome_structure_biogenesis ~ "Melanosome structure / biogenesis",
    Genes %in% melanosome_trafficking ~ "Melanosome trafficking",
    Genes %in% ion_transport_pH ~ "Ion transport / pH regulation",
    Genes %in% regulatory_melanocyte_markers ~ "Regulatory / melanocyte markers",
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
  "Melanin synthesis",
  "Melanosome structure / biogenesis",
  "Melanosome trafficking",
  "Ion transport / pH regulation",
  "Regulatory / melanocyte markers"
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
  "Melanin synthesis" = "#8c510a",
  "Melanosome structure / biogenesis" = "#bf812d",
  "Melanosome trafficking" = "#35978f",
  "Ion transport / pH regulation" = "#01665e",
  "Regulatory / melanocyte markers" = "#762a83"
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
  file.path(output_dir, "melanosome_heatmap_group_means.pdf"),
  width = 7,
  height = 9
)
draw(ht)
dev.off()

svglite(
  file.path(output_dir, "melanosome_heatmap_group_means.svg"),
  width = 7,
  height = 9
)
draw(ht)
dev.off()

write_tsv(
  df_mean,
  file.path(output_dir, "melanosome_genes_used_for_heatmap.tsv")
)

sessionInfo()
