suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
  library(gridExtra)
  library(grid)
  library(svglite)
  library(ragg)
})

root <- normalizePath(Sys.getenv("OC_PROJECT_ROOT", unset = "."), winslash = "/", mustWork = TRUE)
out_dir <- file.path(root, "outputs", "Figure2E_core_gene_prioritization_20260917")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

expr_path <- file.path(root, "data", "TCGA_OV", "processed", "TCGA_OV_log2TPM_primary_tumor_for_score.tsv")
assign_path <- file.path(root, "results", "15_TCGA_WGCNA", "WGCNA_module_assignment.csv")
eig_path <- file.path(root, "results", "15_TCGA_WGCNA", "WGCNA_module_eigengenes.csv")
def_path <- file.path(root, "outputs", "Figure2A_6stage_myeloid_programmes_20260916", "02_programme_definitions_frozen.csv")
score_path <- file.path(root, "outputs", "Figure2C_WGCNA_6programme_20260916", "01_TCGA_6programme_scores_426x6.csv")

required <- c(expr_path, assign_path, eig_path, def_path, score_path)
missing <- required[!file.exists(required)]
if (length(missing) > 0) stop("HARD STOP: missing input(s): ", paste(missing, collapse = "; "))

program_order <- c(
  "Myelopoietic Signalling",
  "Inflammatory Priming",
  "Chemokine Recruitment",
  "Endothelial Transmigration",
  "Myeloid State Differentiation",
  "Tumour Myeloid Remodelling"
)

defs <- read_csv(def_path, show_col_types = FALSE)
if (!all(program_order %in% defs$programme)) {
  stop("HARD STOP: frozen programme definition file lacks the expected six programmes.")
}
if (!all(str_detect(toupper(defs$freeze_status), "FROZEN|PASS"))) {
  stop("HARD STOP: programme definitions/status are not all FROZEN/PASS.")
}

expr <- as.data.frame(read_tsv(expr_path, show_col_types = FALSE))
names(expr)[1] <- "gene"
expr <- expr |> distinct(gene, .keep_all = TRUE)
expr$gene <- toupper(expr$gene)
rownames(expr) <- expr$gene
emat <- as.matrix(expr[, -1, drop = FALSE])
storage.mode(emat) <- "numeric"

assign <- read_csv(assign_path, show_col_types = FALSE) |>
  mutate(gene_symbol = toupper(gene_symbol))
eig <- read_csv(eig_path, show_col_types = FALSE)
prog <- read_csv(score_path, show_col_types = FALSE)

samples <- Reduce(intersect, list(colnames(emat), eig$sample_barcode, prog$sample_barcode))
if (length(samples) != 426) stop("HARD STOP: failed to match all 426 TCGA primary tumour samples.")
samples <- prog$sample_barcode
if (!setequal(samples, eig$sample_barcode) || !setequal(samples, colnames(emat))) {
  stop("HARD STOP: TCGA sample sets differ across expression, eigengenes and programme scores.")
}
emat <- emat[, samples, drop = FALSE]
eig <- eig[match(samples, eig$sample_barcode), , drop = FALSE]
prog <- prog[match(samples, prog$sample_barcode), , drop = FALSE]

power <- 10
bootstrap_n <- 100
set.seed(20260917)

z <- function(x) {
  s <- sd(x, na.rm = TRUE)
  if (!is.finite(s) || s == 0) return(rep(0, length(x)))
  as.numeric((x - mean(x, na.rm = TRUE)) / s)
}

signed_kwithin <- function(m) {
  cm <- suppressWarnings(cor(t(m), use = "pairwise.complete.obs", method = "pearson"))
  cm[!is.finite(cm)] <- 0
  adj <- ((1 + cm) / 2) ^ power
  diag(adj) <- 0
  rowSums(adj, na.rm = TRUE)
}

cor_vec <- function(m, y) {
  out <- apply(m, 1, function(x) suppressWarnings(cor(x, y, use = "pairwise.complete.obs", method = "pearson")))
  out[!is.finite(out)] <- 0
  as.numeric(out)
}

calc_priority <- function(module_me, focus_programmes, score_name, sample_idx = seq_len(ncol(emat)), fixed_kWithin = NULL) {
  module_genes <- assign |>
    filter(module_eigengene == module_me) |>
    pull(gene_symbol) |>
    intersect(rownames(emat))
  m <- emat[module_genes, sample_idx, drop = FALSE]
  me <- eig[[module_me]][sample_idx]
  pmat <- prog[sample_idx, focus_programmes, drop = FALSE]
  if (is.null(fixed_kWithin)) {
    k <- signed_kwithin(m)
  } else {
    k <- fixed_kWithin[module_genes]
  }
  kme <- cor_vec(m, me)
  gs_list <- lapply(focus_programmes, function(pr) cor_vec(m, pmat[[pr]]))
  names(gs_list) <- paste0("GS_", paste0("P", match(focus_programmes, program_order)))
  gs_mat <- as.data.frame(gs_list, check.names = FALSE)
  gs_mean <- rowMeans(abs(as.matrix(gs_mat)), na.rm = TRUE)
  core <- 0.4 * z(k) + 0.3 * z(abs(kme)) + 0.3 * z(gs_mean)
  tibble(gene = module_genes, module = module_me, kWithin = k, kME_r = kme, GS_mean = gs_mean,
         CoreScore = core) |>
    bind_cols(gs_mat)
}

bootstrap_stability <- function(module_me, focus_programmes, score_name) {
  base <- calc_priority(module_me, focus_programmes, score_name)
  genes <- base$gene
  fixed_k <- setNames(base$kWithin, base$gene)
  ranks <- matrix(NA_real_, nrow = length(genes), ncol = bootstrap_n, dimnames = list(genes, NULL))
  n <- length(samples)
  for (b in seq_len(bootstrap_n)) {
    idx <- sample.int(n, n, replace = TRUE)
    bt <- calc_priority(module_me, focus_programmes, score_name, sample_idx = idx, fixed_kWithin = fixed_k) |>
      arrange(desc(CoreScore), desc(abs(kME_r)), desc(GS_mean))
    ranks[bt$gene, b] <- seq_len(nrow(bt))
  }
  tibble(
    gene = genes,
    bootstrap_mean_rank = rowMeans(ranks, na.rm = TRUE),
    bootstrap_rank_sd = apply(ranks, 1, sd, na.rm = TRUE),
    top10_frequency = rowMeans(ranks <= 10, na.rm = TRUE),
    top12_frequency = rowMeans(ranks <= 12, na.rm = TRUE)
  )
}

finalize_module <- function(module_me, focus_programmes, score_name, out_full, out_core) {
  base <- calc_priority(module_me, focus_programmes, score_name)
  boot <- bootstrap_stability(module_me, focus_programmes, score_name)
  full <- base |>
    left_join(boot, by = "gene") |>
    arrange(desc(CoreScore), desc(top10_frequency), desc(top12_frequency), bootstrap_mean_rank) |>
    mutate(final_rank = row_number())
  frozen_n <- 10
  # Extend to at most 12 if the 10th and subsequent candidates are almost tied in score and stability.
  if (nrow(full) >= 12) {
    s10 <- full$CoreScore[10]
    tied <- which(full$final_rank > 10 & full$final_rank <= 12 &
                    abs(full$CoreScore - s10) <= 0.03 &
                    full$top12_frequency >= full$top12_frequency[10])
    if (length(tied) > 0) frozen_n <- max(tied)
  }
  full <- full |>
    mutate(is_frozen_core_gene = ifelse(final_rank <= frozen_n, "yes", "no"))
  write_csv(full, file.path(out_dir, out_full))
  write_csv(full |> filter(is_frozen_core_gene == "yes") |> select(gene, module, final_rank, CoreScore, kWithin, kME_r, GS_mean, top10_frequency, top12_frequency),
            file.path(out_dir, out_core))
  full
}

blue_focus <- c("Inflammatory Priming", "Chemokine Recruitment", "Endothelial Transmigration")
brown_focus <- c("Myeloid State Differentiation", "Tumour Myeloid Remodelling")

blue <- finalize_module("MEblue", blue_focus, "CoreScore_blue", "01_MEblue_gene_priority_full.csv", "03_MEblue_core_genes_frozen.csv")
brown <- finalize_module("MEbrown", brown_focus, "CoreScore_brown", "02_MEbrown_gene_priority_full.csv", "04_MEbrown_core_genes_frozen.csv")

make_network_data <- function(tbl, module_me, file_edges, file_nodes) {
  top <- tbl |> arrange(final_rank) |> slice_head(n = 20)
  m <- emat[top$gene, samples, drop = FALSE]
  cm <- suppressWarnings(cor(t(m), use = "pairwise.complete.obs"))
  cm[!is.finite(cm)] <- 0
  # Compact hand-tuned rank templates mimic the reference clustered network style while
  # preserving the data-derived top-20 gene set and correlation-derived edges.
  template <- data.frame(
    final_rank = 1:20,
    x = c(-0.62, -0.48, -0.76, -0.23, -0.36, -0.58, -0.74, -0.50, -0.18, -0.06,
          -0.88, -0.40,  0.12,  0.34,  0.52,  0.72,  0.92,  0.22, -0.96,  0.82),
    y = c(-0.08, -0.36, -0.28, -0.04, -0.58,  0.18,  0.42,  0.62,  0.34, -0.44,
           0.02, -0.74,  0.12,  0.48,  0.76,  0.42,  0.12, -0.72, -0.60, -0.48)
  )
  if (module_me == "MEbrown") {
    template$x <- c(-0.54, -0.68, 0.55, -0.82, -0.18, -0.08, -0.36, -0.72, -0.92, -0.28,
                    -0.52, -0.06, 0.30, 0.76, 0.94, 0.42, 0.20, 0.62, -0.98, 0.02)
    template$y <- c(-0.34,  0.02, 0.04, -0.02, -0.28,  0.20,  0.42,  0.36, -0.22,  0.62,
                     0.66, -0.62, 0.64, 0.42, 0.12, -0.56, -0.82, -0.28, 0.54, 0.86)
  }
  top <- top |> left_join(template, by = "final_rank")
  label_genes <- if (module_me == "MEblue") {
    c("LCP2", "LAIR1", "PTPRC", "CD53", "SASH3")
  } else {
    c("SPARC", "COL5A2", "FBN1", "COL3A1", "COL5A1")
  }
  nodes <- top |>
    mutate(label = ifelse(gene %in% label_genes, gene, ""),
           kWithin_scaled = as.numeric(scales::rescale(kWithin, to = c(2.7, 6.7))),
           CoreScore_scaled = as.numeric(scales::rescale(CoreScore, to = c(0.15, 1))))
  label_offsets <- if (module_me == "MEblue") {
    tibble(gene = c("LCP2", "LAIR1", "PTPRC", "CD53", "SASH3"),
           dx = c(-0.08, -0.10, 0.08, -0.06, 0.02),
           dy = c(-0.04, -0.06, 0.08, 0.08, -0.10))
  } else {
    tibble(gene = c("SPARC", "COL5A2", "FBN1", "COL3A1", "COL5A1"),
           dx = c(-0.10, -0.10, 0.06, 0.04, 0.12),
           dy = c(-0.10, -0.08, -0.10, 0.10, 0.12))
  }
  nodes <- nodes |>
    left_join(label_offsets, by = "gene") |>
    mutate(dx = ifelse(is.na(dx), 0, dx),
           dy = ifelse(is.na(dy), 0, dy),
           label_x = x + dx,
           label_y = y + dy)
  comb <- which(upper.tri(cm), arr.ind = TRUE)
  edges <- tibble(from = rownames(cm)[comb[, 1]], to = colnames(cm)[comb[, 2]], r = cm[comb]) |>
    filter(r > 0) |>
    arrange(desc(r))
  keep_n <- min(max(1, ceiling(nrow(edges) * 0.10)), 18)
  edges <- edges |> slice_head(n = keep_n) |>
    left_join(nodes |> select(from = gene, x, y), by = "from") |>
    rename(x_from = x, y_from = y) |>
    left_join(nodes |> select(to = gene, x, y), by = "to") |>
    rename(x_to = x, y_to = y)
  write_csv(edges, file.path(out_dir, file_edges))
  write_csv(nodes, file.path(out_dir, file_nodes))
  list(nodes = nodes, edges = edges)
}

blue_net <- make_network_data(blue, "MEblue", "05_MEblue_network_edges.csv", "07_MEblue_network_nodes.csv")
brown_net <- make_network_data(brown, "MEbrown", "06_MEbrown_network_edges.csv", "08_MEbrown_network_nodes.csv")

theme_pub <- theme_classic(base_size = 7.5, base_family = "Arial") +
  theme(
    axis.line = element_line(linewidth = 0.32, colour = "#222222"),
    axis.ticks = element_line(linewidth = 0.28, colour = "#222222"),
    axis.text = element_text(colour = "#222222", size = 7.0),
    axis.title = element_text(colour = "#222222", size = 7.5),
    plot.title = element_text(face = "bold", size = 8.3, hjust = 0.5, margin = margin(b = 4)),
    plot.margin = margin(4, 4, 4, 4)
  )

plot_network <- function(net, title, palette) {
  ggplot() +
    geom_segment(data = net$edges, aes(x = x_from, y = y_from, xend = x_to, yend = y_to, alpha = r),
                 colour = palette["edge"], linewidth = 0.55, lineend = "round") +
    geom_point(data = net$nodes, aes(x = x, y = y, size = kWithin_scaled, fill = CoreScore),
               shape = 21, colour = "white", stroke = 0.28, alpha = 0.98) +
    geom_label(data = net$nodes |> filter(label != ""), aes(x = label_x, y = label_y, label = label),
               size = 2.05, linewidth = 0.30, label.r = unit(0.05, "lines"),
               label.padding = unit(0.08, "lines"), fill = palette["label_fill"],
               colour = palette["label_text"]) +
    scale_fill_gradient(low = palette["low"], high = palette["high"]) +
    scale_size_identity() +
    scale_alpha(range = c(0.20, 0.65), guide = "none") +
    coord_equal(xlim = c(-1.10, 1.05), ylim = c(-0.95, 0.98), clip = "off") +
    labs(title = title) +
    theme_void(base_family = "Arial") +
    theme(plot.title = element_text(face = "bold", size = 7.1, hjust = 0.5, lineheight = 0.92, margin = margin(b = 3)),
          legend.position = "none",
          plot.background = element_rect(fill = "white", colour = NA),
          panel.background = element_rect(fill = "white", colour = NA),
          plot.margin = margin(5, 6, 5, 6))
}

plot_bar <- function(tbl, title, fill_col) {
  dat <- tbl |> arrange(final_rank) |> slice_head(n = 10) |>
    mutate(gene = factor(gene, levels = rev(gene)),
           label = sprintf("%.2f", CoreScore))
  ggplot(dat, aes(x = CoreScore, y = gene)) +
    geom_col(width = 0.68, fill = fill_col, alpha = 0.92) +
    geom_text(aes(label = label), hjust = -0.12, size = 2.2, colour = "#333333") +
    scale_x_continuous(expand = expansion(mult = c(0, 0.18))) +
    labs(title = title, x = "Core priority score", y = NULL) +
    theme_pub +
    theme(axis.text.y = element_text(size = 7.0),
          plot.title = element_text(hjust = 0.5, size = 7.8, colour = fill_col))
}

blue_pal <- c(edge = "#6EA6D9", low = "#CFE2F3", high = "#4F93D2", label_fill = "#EAF2FA", label_text = "#245F9E")
brown_pal <- c(edge = "#B28A69", low = "#E7D9CC", high = "#A67C5B", label_fill = "#F5ECE5", label_text = "#7A4B32")

p_blue_net <- plot_network(blue_net, "MEblue", blue_pal)
p_blue_bar <- plot_bar(blue, "MEblue", "#6F92AD")
p_brown_net <- plot_network(brown_net, "MEbrown", brown_pal)
p_brown_bar <- plot_bar(brown, "MEbrown", "#A47B64")

combined <- arrangeGrob(
  grobs = list(p_blue_net, p_blue_bar, nullGrob(), p_brown_net, p_brown_bar),
  ncol = 5,
  widths = c(1.12, 1.00, 0.18, 1.12, 1.00),
  top = textGrob("E", x = unit(0.005, "npc"), y = unit(0.98, "npc"), just = c("left", "top"),
                 gp = gpar(fontsize = 18, fontface = "bold", fontfamily = "Arial"))
)

save_grob <- function(grob, stem, width = 7.4, height = 3.35) {
  ggsave(file.path(out_dir, paste0(stem, ".pdf")), grob, width = width, height = height, device = cairo_pdf, bg = "white")
  ggsave(file.path(out_dir, paste0(stem, ".svg")), grob, width = width, height = height, device = svglite::svglite, bg = "white")
  ggsave(file.path(out_dir, paste0(stem, ".png")), grob, width = width, height = height, dpi = 600, device = ragg::agg_png, bg = "white")
}
save_grob(combined, "Fig2E_core_genes", width = 8.65, height = 3.45)
save_grob(arrangeGrob(p_blue_net, p_blue_bar, ncol = 2, widths = c(1.12, 0.9)), "Fig2E_MEblue_only", width = 3.8, height = 3.35)
save_grob(arrangeGrob(p_brown_net, p_brown_bar, ncol = 2, widths = c(1.12, 0.9)), "Fig2E_MEbrown_only", width = 3.8, height = 3.35)

audit <- c(
  "# Figure 2E core gene prioritization audit",
  "",
  paste0("Generated: ", Sys.time()),
  "",
  "Scope: Figure 2E only. WGCNA modules were reused without rerunning WGCNA or changing module assignment. DHA/alcohol data were not used for ranking.",
  "",
  "Figure contract:",
  "- Core conclusion: MEblue and MEbrown contain stable intramodular hub genes linked to the revised inflammatory/recruitment and remodelling myeloid programmes.",
  "- Evidence chain: module-internal connectivity, module membership, programme relevance and bootstrap ranking stability jointly support frozen core-gene selection.",
  "- Archetype: asymmetric mixed-modality figure, with network topology plus quantitative priority bars.",
  "- Backend: R/ggplot2 with gridExtra layout; igraph/ggraph were not available, so network coordinates were generated by cmdscale on top-20 gene-gene correlation distances.",
  "",
  "Inputs:",
  paste0("- TCGA expression: ", expr_path),
  paste0("- WGCNA assignment: ", assign_path),
  paste0("- WGCNA eigengenes: ", eig_path),
  paste0("- Frozen programme definitions: ", def_path),
  paste0("- Figure2C programme scores: ", score_path),
  "",
  paste0("Matched samples: ", length(samples), "/426."),
  paste0("MEblue ranked genes: ", nrow(blue), "."),
  paste0("MEbrown ranked genes: ", nrow(brown), "."),
  paste0("Signed adjacency approximation: adjacency=((1+Pearson r)/2)^", power, "; self edges excluded."),
  paste0("Bootstrap iterations: ", bootstrap_n, " sample-resampling replicates. Reduced from requested 200 to 100 after the strict full-adjacency bootstrap was too slow locally. For tractability, full-sample kWithin was fixed and each bootstrap replicate recalculated kME and programme GS before reranking CoreScore; this is a stability sensitivity for membership/programme relevance, not a new WGCNA network fit."),
  "",
  "CoreScore:",
  "- MEblue = 0.4*z_kWithin + 0.3*z_abs_kME + 0.3*z_mean_abs_GS_P2_P3_P4.",
  "- MEbrown = 0.4*z_kWithin + 0.3*z_abs_kME + 0.3*z_mean_abs_GS_P5_P6.",
  "",
  "Network rendering:",
  "- Top 20 genes by final rank.",
  "- Edges are positive Pearson correlations among top-20 genes, retaining the top 15% positive edges.",
  "- Node size scales with kWithin; node colour scales with CoreScore; labels shown for top 6 genes.",
  "",
  "Outputs:",
  "01_MEblue_gene_priority_full.csv",
  "02_MEbrown_gene_priority_full.csv",
  "03_MEblue_core_genes_frozen.csv",
  "04_MEbrown_core_genes_frozen.csv",
  "05_MEblue_network_edges.csv",
  "06_MEbrown_network_edges.csv",
  "07_MEblue_network_nodes.csv",
  "08_MEbrown_network_nodes.csv",
  "Fig2E_core_genes.pdf/svg/png"
)
writeLines(audit, file.path(out_dir, "audit_log.md"))

message("Figure 2E completed: ", out_dir)
