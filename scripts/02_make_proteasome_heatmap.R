# Proteasome / immunoproteasome heatmap
# Input: protein abundance matrix with Genes, SA1-4, V1-4, FL1-4
# Output: PDF and SVG heatmaps

required_packages <- c(
  "readr", "dplyr", "tidyr", "tibble",
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
library(ComplexHeatmap)
library(circlize)
library(grid)
library(svglite)


# Input and output paths
input_file <- "data/processed/data/processed/protein_abundance_matrix_B16_V_FL_SA.tsv"
output_dir <- "outputs/proteasome"

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Read input data

data <- read_tsv(
  input_file,
  na = c("", "NA", "NaN"),
  show_col_types = FALSE,
  progress = FALSE
)

# Sample columns

samples <- c(
  "SA1", "SA2", "SA3", "SA4",
  "V1", "V2", "V3", "V4",
  "FL1", "FL2", "FL3", "FL4"
)

required_cols <- c("Genes", samples)
missing_cols <- setdiff(required_cols, colnames(data))

if (length(missing_cols) > 0) {
  stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
}

data <- data %>%
  select(Genes, all_of(samples)) %>%
  mutate(across(-Genes, as.numeric))

# Proteasome gene sets

alpha20S <- c(
  "Psma1", "Psma2", "Psma3", "Psma4",
  "Psma5", "Psma6", "Psma7"
)

beta20S <- c(
  "Psmb1", "Psmb2", "Psmb3", "Psmb4",
  "Psmb5", "Psmb6", "Psmb7"
)

ATPase19S <- c(
  "Psmc1", "Psmc2", "Psmc3",
  "Psmc4", "Psmc5", "Psmc6"
)

nonATPase19S <- c(
  "Psmd1", "Psmd2", "Psmd3", "Psmd4", "Psmd5", "Psmd6", "Psmd7",
  "Psmd8", "Psmd9", "Psmd10", "Psmd11", "Psmd12", "Psmd13", "Psmd14"
)

immunoproteasome <- c(
  "Psmb8", "Psmb9", "Psmb10"
)

PA28 <- c(
  "Psme1", "Psme2"
)

genes <- unique(c(
  alpha20S,
  beta20S,
  ATPase19S,
  nonATPase19S,
  immunoproteasome,
  PA28
))

# Component annotation

annotation <- tibble(
  Genes = genes,
  Component = case_when(
    Genes %in% immunoproteasome ~ "Immunoproteasome core",
    Genes %in% PA28 ~ "PA28 activator",
    Genes %in% alpha20S ~ "20S core alpha",
    Genes %in% beta20S ~ "20S core beta",
    Genes %in% ATPase19S ~ "19S ATPase",
    Genes %in% nonATPase19S ~ "19S non-ATPase",
    TRUE ~ NA_character_
  )
)

# Extract proteasome proteins from dataset

df <- data %>%
  separate_rows(Genes, sep = ";") %>%
  mutate(Genes = trimws(Genes)) %>%
  filter(Genes %in% genes) %>%
  group_by(Genes) %>%
  summarise(across(all_of(samples), ~ mean(.x, na.rm = TRUE)), .groups = "drop") %>%
  left_join(annotation, by = "Genes")

# Keep biologically defined component order

component_levels <- c(
  "Immunoproteasome core",
  "PA28 activator",
  "20S core alpha",
  "20S core beta",
  "19S ATPase",
  "19S non-ATPase"
)

df <- df %>%
  mutate(Component = factor(Component, levels = component_levels)) %>%
  arrange(Component, Genes)

# Calculate group means

df_mean <- df %>%
  mutate(
    V  = rowMeans(across(c(V1, V2, V3, V4)), na.rm = TRUE),
    FL = rowMeans(across(c(FL1, FL2, FL3, FL4)), na.rm = TRUE),
    SA = rowMeans(across(c(SA1, SA2, SA3, SA4)), na.rm = TRUE)
  ) %>%
  select(Genes, Component, V, FL, SA)

# Build matrix

mat <- df_mean %>%
  select(Genes, V, FL, SA) %>%
  column_to_rownames("Genes") %>%
  as.matrix()

mat <- log2(mat + 1)

row_zscore <- function(x) {
  z <- t(scale(t(x)))
  z[is.na(z)] <- 0
  z
}

mat_z <- row_zscore(mat)

# Row annotation

row_components <- df_mean$Component
names(row_components) <- df_mean$Genes
row_components <- row_components[rownames(mat_z)]

component_colors <- c(
  "Immunoproteasome core" = "#b2182b",
  "PA28 activator" = "#ef8a62",
  "20S core alpha" = "#1f78b4",
  "20S core beta" = "#33a02c",
  "19S ATPase" = "#6a3d9a",
  "19S non-ATPase" = "#a6cee3"
)

row_anno <- rowAnnotation(
  Component = row_components,
  col = list(Component = component_colors),
  show_annotation_name = FALSE
)

# Z-score color scale

col_fun <- colorRamp2(
  c(-2, 0, 2),
  c("#2166AC", "white", "#B2182B")
)

# Build heatmap

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
  column_labels = c("V", "FL", "SA"),
  row_names_gp = gpar(fontsize = 8),
  column_names_gp = gpar(fontsize = 11, fontface = "bold"),
  column_title = NULL,
  row_title = NULL
)

# Save outputs

pdf(
  file.path(output_dir, "proteasome_heatmap_group_means.pdf"),
  width = 7,
  height = 10
)
draw(ht)
dev.off()

svglite(
  file.path(output_dir, "proteasome_heatmap_group_means.svg"),
  width = 7,
  height = 10
)
draw(ht)
dev.off()

write_tsv(
  df_mean,
  file.path(output_dir, "proteasome_genes_used_for_heatmap.tsv")
)

sessionInfo()
