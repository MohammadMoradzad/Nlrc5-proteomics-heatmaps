# Autophagy heatmap for B16 proteomics data
# Input: protein abundance matrix with Genes, SA1-4, V1-4, FL1-4
# Output: PDF and SVG heatmaps

required_packages <- c(
  "readr", "dplyr", "tidyr", "stringr", "tibble",
  "ComplexHeatmap", "circlize", "grid", "svglite", "KEGGREST"
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
library(stringr)
library(tibble)
library(ComplexHeatmap)
library(circlize)
library(grid)
library(svglite)
library(KEGGREST)

input_file <- "data/processed/protein_abundance_matrix_B16_V_FL_SA.tsv"
output_dir <- "outputs/autophagy"

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

sample_cols <- c(
  "V1", "V2", "V3", "V4",
  "FL1", "FL2", "FL3", "FL4",
  "SA1", "SA2", "SA3", "SA4"
)

raw_df <- read_tsv(
  input_file,
  na = c("", "NA", "NaN"),
  show_col_types = FALSE,
  progress = FALSE
)

required_cols <- c("Genes", sample_cols)
missing_cols <- setdiff(required_cols, colnames(raw_df))

if (length(missing_cols) > 0) {
  stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
}

df <- raw_df %>%
  select(Genes, all_of(sample_cols)) %>%
  mutate(across(all_of(sample_cols), as.numeric))

get_kegg_symbols <- function(pathway_id) {
  pathway <- KEGGREST::keggGet(pathway_id)[[1]]
  gene_vec <- pathway$GENE

  if (is.null(gene_vec)) {
    stop("No GENE field returned for ", pathway_id)
  }

  symbols <- unname(gene_vec)
  symbols <- symbols[grepl(";", symbols)]
  symbols <- symbols %>%
    str_replace(";.*$", "") %>%
    str_trim() %>%
    unique() %>%
    sort()

  symbols
}

autophagy_genes <- get_kegg_symbols("path:mmu04140")

annotate_autophagy <- function(gene) {
  case_when(
    gene %in% c("Ulk1", "Rb1cc1", "Atg13") ~
      "ULK initiation complex",

    gene %in% c("Pik3c3", "Atg14", "Ambra1", "Rubcn", "Zfyve1", "Becn1") ~
      "PI3K-BECN1 nucleation",

    gene %in% c(
      "Atg3", "Atg4b", "Atg4d", "Atg5", "Atg7", "Atg9a", "Atg9b",
      "Atg10", "Atg12", "Atg16l1", "Map1lc3a", "Map1lc3b",
      "Gabarap", "Gabarapl1", "Gabarapl2", "Wipi1", "Wipi2",
      "Uba52", "Ubb", "Ubc"
    ) ~
      "ATG conjugation / elongation",

    gene %in% c("Sqstm1", "Nbr1", "Optn", "Tax1bp1", "Calcoco2") ~
      "Cargo receptors / selective autophagy",

    gene %in% c(
      "Rab1a", "Rab7", "Rab7b", "Rab8a", "Rab33b", "Vamp8",
      "Epg5", "Plekhm1", "Vps18", "Vps41", "Wdr41", "Smcr8",
      "Stx17", "Snap29"
    ) ~
      "Maturation / fusion",

    gene %in% c("Lamp1", "Lamp2", "Ctsb", "Ctsd", "Ctsl", "Ctss", "Ctsz", "Lgmn") ~
      "Lysosomal degradation arm",

    gene %in% c(
      "Akt1", "Akt2", "Akt3", "Prkaa1", "Prkaa2", "Pik3ca", "Pik3cd",
      "Pik3r1", "Pik3r2", "Pik3r3", "Pten", "Pdpk1", "Rheb", "Tsc2",
      "Stk11", "Camkk2", "Raf1", "Hras", "Kras", "Nras", "Igf1r",
      "Irs1", "Irs3", "Irs4", "Prkaca", "Prkacb", "Ppp2ca", "Ppp2cb",
      "Eif2ak3", "Eif2s1", "Hif1a", "Traf6", "Tank", "Mtor", "Rragc"
    ) ~
      "Upstream signaling",

    TRUE ~ "Other autophagy KEGG genes"
  )
}

df_long <- df %>%
  separate_rows(Genes, sep = ";") %>%
  mutate(Genes = str_trim(Genes)) %>%
  filter(!is.na(Genes), Genes != "")

autophagy_df <- df_long %>%
  filter(Genes %in% autophagy_genes) %>%
  group_by(Genes) %>%
  summarise(across(all_of(sample_cols), ~ mean(.x, na.rm = TRUE)), .groups = "drop") %>%
  mutate(Component = annotate_autophagy(Genes)) %>%
  filter(!Component %in% c("Upstream signaling", "Other autophagy KEGG genes"))

component_levels <- c(
  "ULK initiation complex",
  "PI3K-BECN1 nucleation",
  "ATG conjugation / elongation",
  "Cargo receptors / selective autophagy",
  "Maturation / fusion",
  "Lysosomal degradation arm"
)

autophagy_df <- autophagy_df %>%
  mutate(Component = factor(Component, levels = component_levels)) %>%
  arrange(Component, Genes)

mat <- autophagy_df %>%
  select(Genes, all_of(sample_cols)) %>%
  column_to_rownames("Genes") %>%
  as.matrix()

mat <- log2(mat + 1)

row_zscore <- function(x) {
  z <- t(scale(t(x)))
  z[is.na(z)] <- 0
  z
}

mat_z <- row_zscore(mat)

row_components <- autophagy_df$Component
names(row_components) <- autophagy_df$Genes
row_components <- row_components[rownames(mat_z)]

component_colors <- c(
  "ULK initiation complex" = "#d95f02",
  "PI3K-BECN1 nucleation" = "#7570b3",
  "ATG conjugation / elongation" = "#e7298a",
  "Cargo receptors / selective autophagy" = "#66a61e",
  "Maturation / fusion" = "#e6ab02",
  "Lysosomal degradation arm" = "#a6761d"
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
  column_labels = sample_cols,
  row_names_gp = gpar(fontsize = 8),
  column_names_gp = gpar(fontsize = 10, fontface = "bold"),
  column_title = NULL,
  row_title = NULL
)

pdf(file.path(output_dir, "autophagy_heatmap_replicates.pdf"), width = 8, height = 10)
draw(ht)
dev.off()

svglite(file.path(output_dir, "autophagy_heatmap_replicates.svg"), width = 8, height = 10)
draw(ht)
dev.off()

write_tsv(
  autophagy_df,
  file.path(output_dir, "autophagy_genes_used_for_heatmap.tsv")
)

sessionInfo()
