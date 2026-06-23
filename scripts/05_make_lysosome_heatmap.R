# Lysosome / lytic vacuole heatmap
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
output_dir <- "outputs/lysosome"

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

cathepsins <- c(
  "Ctsb", "Ctsd", "Ctsl", "Ctss", "Ctsz", "Ctsc", "Ctsh", "Ctsk", "Ctso"
)

lysosomal_membrane <- c(
  "Lamp1", "Lamp2", "Laptm4a", "Laptm4b", "Laptm5", "Cd63", "Scarb2", "Mcoln1"
)

v_atpase <- c(
  "Atp6v0a1", "Atp6v0a2", "Atp6v0a4", "Atp6v0b", "Atp6v0c", "Atp6v0d1",
  "Atp6v0e", "Atp6v1a", "Atp6v1b2", "Atp6v1c1", "Atp6v1d", "Atp6v1e1",
  "Atp6v1f", "Atp6v1g1", "Atp6v1h"
)

glycosidases_hydrolases <- c(
  "Gba", "Gaa", "Hexa", "Hexb", "Man2b1", "Neu1", "Gusb", "Glb1", "Arsa", "Arsb",
  "Acp2", "Acp5", "Dnase2a", "Psap", "Smpd1", "Tpp1"
)

trafficking_sorting <- c(
  "Clta", "Cltb", "Cltc", "Ap1b1", "Ap1g1", "Ap1m1", "Ap1s1",
  "Ap3b1", "Ap3d1", "Ap3m1", "Ap3s1", "Vps11", "Vps16", "Vps18",
  "Vps33a", "Vps39", "Vps41", "Rab7", "Rab7b", "Rab9a", "Rab5a"
)

target_genes <- unique(c(
  cathepsins,
  lysosomal_membrane,
  v_atpase,
  glycosidases_hydrolases,
  trafficking_sorting
))

annotation <- tibble(
  Genes = target_genes,
  Component = case_when(
    Genes %in% cathepsins ~ "Cathepsins / proteases",
    Genes %in% lysosomal_membrane ~ "Lysosomal membrane proteins",
    Genes %in% v_atpase ~ "V-ATPase / acidification",
    Genes %in% glycosidases_hydrolases ~ "Glycosidases / hydrolases",
    Genes %in% trafficking_sorting ~ "Trafficking / sorting",
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
  "Cathepsins / proteases",
  "Lysosomal membrane proteins",
  "V-ATPase / acidification",
  "Glycosidases / hydrolases",
  "Trafficking / sorting"
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
  "Cathepsins / proteases" = "#b2182b",
  "Lysosomal membrane proteins" = "#ef8a62",
  "V-ATPase / acidification" = "#2166ac",
  "Glycosidases / hydrolases" = "#1b7837",
  "Trafficking / sorting" = "#762a83"
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
  file.path(output_dir, "lysosome_heatmap_group_means.pdf"),
  width = 7,
  height = 10
)
draw(ht)
dev.off()

svglite(
  file.path(output_dir, "lysosome_heatmap_group_means.svg"),
  width = 7,
  height = 10
)
draw(ht)
dev.off()

write_tsv(
  df_mean,
  file.path(output_dir, "lysosome_genes_used_for_heatmap.tsv")
)

sessionInfo()
