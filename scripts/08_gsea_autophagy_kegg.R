# KEGG Autophagy GSEA plot
# Input: protein abundance matrix with Genes, V1-4, N5FL1-4, N5CA1-4
# Output: GSEA plots for N5FL vs V and N5CA vs V

required_packages <- c(
  "readr", "dplyr", "tidyr", "tibble", "stringr",
  "ggplot2", "ggrepel", "limma", "fgsea",
  "KEGGREST", "org.Mm.eg.db", "AnnotationDbi"
)

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

cran_packages <- c(
  "readr", "dplyr", "tidyr", "tibble", "stringr",
  "ggplot2", "ggrepel"
)

bioc_packages <- c(
  "limma", "fgsea", "KEGGREST", "org.Mm.eg.db", "AnnotationDbi"
)

missing_cran <- cran_packages[
  !cran_packages %in% rownames(installed.packages())
]

if (length(missing_cran) > 0) {
  install.packages(missing_cran)
}

missing_bioc <- bioc_packages[
  !bioc_packages %in% rownames(installed.packages())
]

if (length(missing_bioc) > 0) {
  BiocManager::install(missing_bioc)
}

library(readr)
library(dplyr)
library(tidyr)
library(tibble)
library(stringr)
library(ggplot2)
library(ggrepel)
library(limma)
library(fgsea)
library(KEGGREST)
library(org.Mm.eg.db)
library(AnnotationDbi)

input_file <- "data/processed/data/processed/protein_abundance_matrix_B16_V_N5FL_N5CA.tsv"
output_dir <- "outputs/gsea/autophagy"

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

gene_df <- data %>%
  separate_rows(Genes, sep = ";") %>%
  mutate(Genes = str_trim(Genes)) %>%
  filter(!is.na(Genes), Genes != "") %>%
  mutate(Genes = toupper(Genes)) %>%
  group_by(Genes) %>%
  summarise(across(all_of(samples), ~ mean(.x, na.rm = TRUE)), .groups = "drop")

expr <- gene_df %>%
  column_to_rownames("Genes") %>%
  as.matrix()

expr <- log2(expr + 1)

for (i in seq_len(nrow(expr))) {
  nas <- is.na(expr[i, ])
  if (any(nas)) {
    expr[i, nas] <- mean(expr[i, ], na.rm = TRUE)
  }
}

row_sd <- apply(expr, 1, sd, na.rm = TRUE)
expr <- expr[row_sd > 0, , drop = FALSE]

group <- factor(
  c(
    rep("V", 4),
    rep("N5FL", 4),
    rep("N5CA", 4)
  ),
  levels = c("V", "N5FL", "N5CA")
)

design <- model.matrix(~ 0 + group)
colnames(design) <- levels(group)

fit <- lmFit(expr, design)

highlight_genes <- c(
  "ULK1", "RB1CC1", "ATG13",
  "BECN1", "PIK3C3", "ATG14",
  "ATG3", "ATG5", "ATG7", "ATG12", "ATG16L1",
  "MAP1LC3A", "MAP1LC3B",
  "GABARAP", "GABARAPL1", "GABARAPL2",
  "SQSTM1", "NBR1", "OPTN",
  "STX17", "SNAP29", "VAMP8"
)

highlight_genes <- toupper(highlight_genes)

kegg_links <- KEGGREST::keggLink("mmu", "path:mmu04140")
kegg_entrez <- unique(sub("^mmu:", "", unname(kegg_links)))

map_tbl <- AnnotationDbi::select(
  org.Mm.eg.db,
  keys = kegg_entrez,
  keytype = "ENTREZID",
  columns = c("ENTREZID", "SYMBOL")
)

map_tbl <- as.data.frame(map_tbl, stringsAsFactors = FALSE)
map_tbl <- map_tbl[!is.na(map_tbl$SYMBOL), , drop = FALSE]
map_tbl$SYMBOL <- toupper(map_tbl$SYMBOL)

autophagy_genes <- unique(map_tbl$SYMBOL)

cat("Autophagy genes in KEGG:", length(autophagy_genes), "\n")
cat("Highlighted genes present in expression matrix:\n")
print(intersect(highlight_genes, rownames(expr)))

make_ranks_from_contrast <- function(fit, design, contrast_formula) {
  contrast_matrix <- makeContrasts(contrasts = contrast_formula, levels = design)
  fit2 <- contrasts.fit(fit, contrast_matrix)
  fit2 <- eBayes(fit2)

  tt <- topTable(fit2, number = Inf, sort.by = "none")

  ranks <- tt$t
  names(ranks) <- rownames(tt)

  keep <- !is.na(ranks) & !is.na(names(ranks)) & names(ranks) != ""
  ranks <- ranks[keep]

  ranks <- ranks + seq_along(ranks) * 1e-12
  ranks <- sort(ranks, decreasing = TRUE)

  list(
    ranks = ranks,
    fit2 = fit2,
    table = tt
  )
}

compute_running_es <- function(stats, pathway_genes, gsea_weight = 1) {
  stats <- stats[!is.na(stats)]
  stats <- sort(stats, decreasing = TRUE)

  genes <- names(stats)
  hits <- genes %in% pathway_genes

  N <- length(stats)
  Nh <- sum(hits)
  Nm <- N - Nh

  if (Nh == 0) stop("No pathway genes found in ranked list.")
  if (Nm == 0) stop("All genes are pathway genes.")

  hit_weights <- abs(stats[hits])^gsea_weight
  hit_weights <- hit_weights / sum(hit_weights)

  increments <- rep(-1 / Nm, N)
  increments[hits] <- hit_weights

  running_es <- cumsum(increments)

  data.frame(
    Rank = seq_len(N),
    Gene = genes,
    RunningES = running_es,
    Hit = hits,
    stringsAsFactors = FALSE
  )
}

make_autophagy_gsea_plot <- function(
  ranks,
  pathway_genes,
  highlight_genes,
  comparison_name
) {
  pathway_genes <- intersect(pathway_genes, names(ranks))

  if (length(pathway_genes) < 10) {
    stop(paste("Too few KEGG autophagy genes overlap for", comparison_name))
  }

  es_df <- compute_running_es(ranks, pathway_genes, gsea_weight = 1)

  fg <- fgseaMultilevel(
    pathways = list(Autophagy_KEGG = pathway_genes),
    stats = ranks,
    minSize = 10,
    maxSize = 500
  )

  fg <- as.data.frame(fg)

  NES_val <- round(fg$NES[1], 3)
  FDR_val <- signif(fg$padj[1], 3)
  SIZE_val <- fg$size[1]

  barcode_df <- es_df[es_df$Hit, c("Rank", "Gene"), drop = FALSE]

  highlight_present <- intersect(highlight_genes, names(ranks))
  highlight_in_pathway <- intersect(highlight_present, pathway_genes)

  cat("\n", comparison_name, "\n", sep = "")
  cat("Highlighted genes present in ranked list:\n")
  print(highlight_present)
  cat("Highlighted genes present in KEGG autophagy pathway:\n")
  print(highlight_in_pathway)

  highlight_df <- barcode_df[
    barcode_df$Gene %in% highlight_in_pathway,
    ,
    drop = FALSE
  ]

  if (nrow(highlight_df) > 0) {
    highlight_df$order_id <- match(highlight_df$Gene, highlight_genes)
    highlight_df <- highlight_df[order(highlight_df$order_id, highlight_df$Rank), ]
  }

  y_min <- min(es_df$RunningES)
  y_max <- max(es_df$RunningES)
  y_range <- y_max - y_min

  black_y0 <- y_min - 0.11 * y_range
  black_y1 <- y_min - 0.02 * y_range

  red_y0 <- y_min - 0.17 * y_range
  red_y1 <- y_min - 0.02 * y_range

  if (nrow(highlight_df) > 0) {
    highlight_df$label_y <- red_y0 - 0.05 * y_range
  }

  stats_line <- paste0(
    "NES = ", NES_val,
    "   |   FDR = ", FDR_val,
    "   |   pathway size = ", SIZE_val
  )

  tick_line <- "Black ticks = KEGG autophagy genes   |   Red ticks = highlighted autophagy genes"

  p <- ggplot(es_df, aes(x = Rank, y = RunningES)) +
    geom_line(
      color = "#00C853",
      linewidth = 1.15,
      lineend = "round"
    ) +
    geom_hline(
      yintercept = 0,
      linetype = "dashed",
      linewidth = 0.7,
      color = "black"
    ) +
    geom_segment(
      data = barcode_df,
      aes(x = Rank, xend = Rank, y = black_y0, yend = black_y1),
      inherit.aes = FALSE,
      color = "black",
      linewidth = 0.28,
      alpha = 0.9
    ) +
    geom_segment(
      data = highlight_df,
      aes(x = Rank, xend = Rank, y = red_y0, yend = red_y1),
      inherit.aes = FALSE,
      color = "#D62728",
      linewidth = 0.8
    ) +
    ggrepel::geom_text_repel(
      data = highlight_df,
      aes(x = Rank, y = label_y, label = Gene),
      inherit.aes = FALSE,
      color = "#D62728",
      size = 3.2,
      angle = 90,
      direction = "y",
      force = 2.5,
      box.padding = 0.18,
      point.padding = 0.05,
      min.segment.length = 0,
      segment.color = "#D62728",
      segment.size = 0.25,
      max.overlaps = Inf,
      seed = 123
    ) +
    annotate(
      "text",
      x = min(es_df$Rank),
      y = y_max + 0.18 * y_range,
      label = stats_line,
      hjust = 0,
      vjust = 0,
      size = 4.7
    ) +
    annotate(
      "text",
      x = min(es_df$Rank),
      y = y_max + 0.10 * y_range,
      label = tick_line,
      hjust = 0,
      vjust = 0,
      size = 3.7
    ) +
    labs(
      title = paste0("KEGG GSEA: Autophagy (", comparison_name, ")"),
      x = "Rank in ordered proteins",
      y = "Running enrichment score"
    ) +
    coord_cartesian(
      ylim = c(red_y0 - 0.20 * y_range, y_max + 0.26 * y_range),
      clip = "off"
    ) +
    theme_classic(base_size = 16) +
    theme(
      plot.title = element_text(size = 22, face = "plain", hjust = 0),
      axis.title = element_text(size = 18),
      axis.text = element_text(size = 13),
      plot.margin = margin(t = 28, r = 26, b = 95, l = 26)
    )

  list(
    plot = p,
    fgsea = fg,
    highlighted = highlight_df,
    running_es = es_df
  )
}

res_N5FL_vs_V <- make_ranks_from_contrast(
  fit = fit,
  design = design,
  contrast_formula = "N5FL - V"
)

res_N5CA_vs_V <- make_ranks_from_contrast(
  fit = fit,
  design = design,
  contrast_formula = "N5CA - V"
)

plot_N5FL <- make_autophagy_gsea_plot(
  ranks = res_N5FL_vs_V$ranks,
  pathway_genes = autophagy_genes,
  highlight_genes = highlight_genes,
  comparison_name = "N5FL_vs_V"
)

plot_N5CA <- make_autophagy_gsea_plot(
  ranks = res_N5CA_vs_V$ranks,
  pathway_genes = autophagy_genes,
  highlight_genes = highlight_genes,
  comparison_name = "N5CA_vs_V"
)

ggsave(
  file.path(output_dir, "Autophagy_KEGG_GSEA_N5FL_vs_V.pdf"),
  plot_N5FL$plot,
  width = 11,
  height = 7.2
)

ggsave(
  file.path(output_dir, "Autophagy_KEGG_GSEA_N5FL_vs_V.png"),
  plot_N5FL$plot,
  width = 11,
  height = 7.2,
  dpi = 600
)

ggsave(
  file.path(output_dir, "Autophagy_KEGG_GSEA_N5CA_vs_V.pdf"),
  plot_N5CA$plot,
  width = 11,
  height = 7.2
)

ggsave(
  file.path(output_dir, "Autophagy_KEGG_GSEA_N5CA_vs_V.png"),
  plot_N5CA$plot,
  width = 11,
  height = 7.2,
  dpi = 600
)

write_tsv(
  plot_N5FL$fgsea,
  file.path(output_dir, "Autophagy_KEGG_GSEA_N5FL_vs_V_fgsea_results.tsv")
)

write_tsv(
  plot_N5CA$fgsea,
  file.path(output_dir, "Autophagy_KEGG_GSEA_N5CA_vs_V_fgsea_results.tsv")
)

write_tsv(
  plot_N5FL$highlighted,
  file.path(output_dir, "Autophagy_KEGG_GSEA_N5FL_vs_V_highlighted_genes.tsv")
)

write_tsv(
  plot_N5CA$highlighted,
  file.path(output_dir, "Autophagy_KEGG_GSEA_N5CA_vs_V_highlighted_genes.tsv")
)

print(plot_N5FL$plot)
print(plot_N5CA$plot)

sessionInfo()
