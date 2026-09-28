#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(grid)
  library(svglite)
  library(ragg)
})

project <- normalizePath(Sys.getenv("OC_PROJECT_ROOT", unset = "."), winslash = "/", mustWork = TRUE)
out_dir <- file.path(project, "outputs", "Figure3A_MEblue_MEbrown_enrichment_20260917")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

module_path <- file.path(project, "results", "15_TCGA_WGCNA", "WGCNA_module_assignment.csv")
hallmark_gmt <- file.path(project, "data", "gene_sets", "h.all.v2025.1.Hs.symbols.gmt")
reactome_gmt <- file.path(project, "data", "gene_sets", "c2.cp.reactome.v2025.1.Hs.symbols.gmt")
gobp_gmt <- file.path(project, "data", "gene_sets", "c5.go.bp.v2025.1.Hs.symbols.gmt")

stopifnot(file.exists(module_path), file.exists(hallmark_gmt), file.exists(reactome_gmt), file.exists(gobp_gmt))

clean_symbol <- function(x) unique(toupper(trimws(as.character(x[!is.na(x) & x != ""]))))

read_gmt <- function(path, source) {
  lines <- readLines(path, warn = FALSE)
  rows <- lapply(lines, function(line) {
    x <- strsplit(line, "\t", fixed = TRUE)[[1]]
    data.frame(
      source = source,
      term_id = x[1],
      pathway = x[1],
      gene = unique(toupper(x[-c(1, 2)])),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

ora_one <- function(gmt, target, universe, module, source) {
  terms <- split(gmt$gene, gmt$term_id)
  M <- length(universe)
  N <- length(target)
  rows <- lapply(names(terms), function(term_id) {
    gs <- intersect(unique(terms[[term_id]]), universe)
    overlap_genes <- sort(intersect(target, gs))
    k <- length(overlap_genes)
    n <- length(gs)
    if (n < 5 || k == 0) return(NULL)
    p <- phyper(k - 1, n, M - n, N, lower.tail = FALSE)
    data.frame(
      module = module,
      source = source,
      term_id = term_id,
      pathway = term_id,
      overlap = k,
      term_size = n,
      module_gene_n = N,
      background_n = M,
      gene_ratio = k / N,
      BackgroundRatio = n / M,
      P = p,
      genes = paste(overlap_genes, collapse = ";"),
      direction = NA_character_,
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  if (is.null(out)) out <- data.frame()
  out$BH_FDR_q <- p.adjust(out$P, method = "BH")
  out <- out[order(out$BH_FDR_q, -out$gene_ratio, out$pathway), ]
  rownames(out) <- NULL
  out
}

jaccard <- function(a, b) {
  if (!length(a) && !length(b)) return(0)
  length(intersect(a, b)) / length(union(a, b))
}

reduce_overlap <- function(df, source_label, threshold) {
  if (!nrow(df)) return(list(keep = df, map = data.frame()))
  sig <- df[df$BH_FDR_q < 0.05, ]
  nonsig <- df[df$BH_FDR_q >= 0.05 | is.na(df$BH_FDR_q), ]
  sig <- sig[order(sig$BH_FDR_q, -sig$gene_ratio, sig$pathway), ]
  kept <- rep(FALSE, nrow(sig))
  removed <- rep(FALSE, nrow(sig))
  maps <- list()
  for (i in seq_len(nrow(sig))) {
    if (removed[i]) next
    kept[i] <- TRUE
    gi <- strsplit(sig$genes[i], ";", fixed = TRUE)[[1]]
    for (j in seq_len(nrow(sig))) {
      if (j <= i || removed[j]) next
      gj <- strsplit(sig$genes[j], ";", fixed = TRUE)[[1]]
      jac <- jaccard(gi, gj)
      if (jac >= threshold) {
        removed[j] <- TRUE
        maps[[length(maps) + 1]] <- data.frame(
          module = sig$module[i],
          source = source_label,
          representative_pathway = sig$pathway[i],
          merged_pathway = sig$pathway[j],
          overlap_jaccard = jac,
          reason = paste0("overlap_jaccard>=", threshold),
          stringsAsFactors = FALSE
        )
      }
    }
  }
  keep_sig <- sig[kept & !removed, ]
  keep <- rbind(keep_sig, nonsig)
  keep <- keep[order(keep$BH_FDR_q, -keep$gene_ratio, keep$pathway), ]
  rownames(keep) <- NULL
  mapping <- if (length(maps)) {
    do.call(rbind, maps)
  } else {
    data.frame(module=character(), source=character(), representative_pathway=character(),
               merged_pathway=character(), overlap_jaccard=numeric(), reason=character())
  }
  list(keep = keep, map = mapping)
}

clean_label <- function(x) {
  specials <- c(
    HALLMARK_TNFA_SIGNALING_VIA_NFKB = "TNF\u03B1 signaling via NF-\u03BAB",
    HALLMARK_IL6_JAK_STAT3_SIGNALING = "IL6\u2013JAK\u2013STAT3 signaling",
    REACTOME_CHEMOKINE_RECEPTORS_BIND_CHEMOKINES = "Chemokine receptors bind chemokines",
    GOBP_LEUKOCYTE_CELL_CELL_ADHESION = "Leukocyte cell-cell adhesion",
    GOBP_REGULATION_OF_IMMUNE_SYSTEM_PROCESS = "Regulation of immune system process"
  )
  out <- unname(specials[x])
  idx <- is.na(out)
  if (any(idx)) {
    z <- x[idx]
    z <- gsub("^(HALLMARK|REACTOME|GOBP)_", "", z)
    z <- gsub("_", " ", z)
    z <- tolower(z)
    z <- paste0(toupper(substr(z, 1, 1)), substr(z, 2, nchar(z)))
    z <- gsub("\\bNfkb\\b", "NF-\u03BAB", z)
    z <- gsub("\\bTnf\\b", "TNF", z)
    z <- gsub("\\bIl6\\b", "IL6", z)
    z <- gsub("\\bJak\\b", "JAK", z)
    z <- gsub("\\bStat3\\b", "STAT3", z)
    out[idx] <- z
  }
  out
}

wrap_label <- function(x, width = 42) {
  vapply(strwrap(x, width = width, simplify = FALSE), paste, collapse = "\n", character(1))
}

module <- read.csv(module_path, stringsAsFactors = FALSE, check.names = FALSE)
if (!all(c("gene_symbol", "module") %in% names(module))) {
  stop("WGCNA_module_assignment.csv must contain gene_symbol and module columns.")
}
universe <- clean_symbol(module$gene_symbol)
targets <- list(
  MEblue = clean_symbol(module$gene_symbol[tolower(module$module) == "blue"]),
  MEbrown = clean_symbol(module$gene_symbol[tolower(module$module) == "brown"])
)
if (any(vapply(targets, length, integer(1)) == 0)) stop("MEblue or MEbrown gene list is empty.")

hallmark <- read_gmt(hallmark_gmt, "Hallmark")
reactome <- read_gmt(reactome_gmt, "Reactome")
gobp <- read_gmt(gobp_gmt, "GO Biological Process")
all_gmt <- list(Hallmark = hallmark, Reactome = reactome, `GO Biological Process` = gobp)

full_results <- list()
nonred_results <- list()
maps <- list()
selected <- list()

for (mod in names(targets)) {
  target <- targets[[mod]]
  res_h <- ora_one(hallmark, target, universe, mod, "Hallmark")
  res_r <- ora_one(reactome, target, universe, mod, "Reactome")
  res_g <- ora_one(gobp, target, universe, mod, "GO Biological Process")
  full_results[[mod]] <- rbind(res_h, res_r, res_g)
  red_h <- list(keep = res_h, map = data.frame())
  red_r <- reduce_overlap(res_r, "Reactome", threshold = 0.50)
  red_g <- reduce_overlap(res_g, "GO Biological Process", threshold = 0.35)
  nonred <- rbind(red_h$keep, red_r$keep, red_g$keep)
  nonred <- nonred[order(nonred$BH_FDR_q, -nonred$gene_ratio, nonred$pathway), ]
  nonred_results[[mod]] <- nonred
  maps[[mod]] <- rbind(red_r$map, red_g$map)
  sig <- nonred[nonred$BH_FDR_q < 0.05, ]
  selected[[mod]] <- head(sig[order(sig$BH_FDR_q, -sig$gene_ratio, sig$pathway), ], 8)
}

all_full <- do.call(rbind, full_results)
all_nonred <- do.call(rbind, nonred_results)
all_maps <- do.call(rbind, maps)
plot_data <- do.call(rbind, selected)
rownames(plot_data) <- NULL

plot_data$minus_log10_FDR <- -log10(pmax(plot_data$BH_FDR_q, .Machine$double.xmin))
plot_data$display_label <- clean_label(plot_data$pathway)
plot_data$plot_label <- wrap_label(plot_data$display_label, width = 40)
plot_data$module <- factor(plot_data$module, levels = c("MEblue", "MEbrown"))
plot_data <- plot_data[order(plot_data$module, plot_data$BH_FDR_q, -plot_data$gene_ratio, plot_data$pathway), ]
plot_data$display_order <- ave(seq_len(nrow(plot_data)), plot_data$module, FUN = seq_along)
plot_data$term_axis <- paste(plot_data$module, plot_data$plot_label, plot_data$pathway, sep = "|||")
plot_data$term_axis <- factor(plot_data$term_axis, levels = rev(plot_data$term_axis))

write.csv(all_full, file.path(out_dir, "Fig3A_MEblue_MEbrown_enrichment_full_ORA.csv"), row.names = FALSE)
write.csv(all_nonred, file.path(out_dir, "Fig3A_MEblue_MEbrown_enrichment_nonredundant.csv"), row.names = FALSE)
write.csv(all_maps, file.path(out_dir, "Fig3A_MEblue_MEbrown_enrichment_redundancy_mapping.csv"), row.names = FALSE)
write.csv(plot_data, file.path(out_dir, "Fig3A_MEblue_MEbrown_enrichment_source.csv"), row.names = FALSE)

axis_labeller <- function(x) sub("^.*?\\|\\|\\|(.*?)\\|\\|\\|.*$", "\\1", x)

p <- ggplot(plot_data, aes(x = gene_ratio, y = term_axis)) +
  geom_point(aes(size = overlap, fill = minus_log10_FDR), shape = 21, colour = "#FFFFFF", stroke = 0.18, alpha = 0.86) +
  facet_wrap(~ module, scales = "free_y", nrow = 1) +
  scale_y_discrete(labels = axis_labeller) +
  scale_x_continuous(name = "Gene ratio", expand = expansion(mult = c(0.02, 0.10))) +
  scale_size_continuous(
    name = "Overlap",
    range = c(2.2, 7.1),
    breaks = pretty(range(plot_data$overlap), n = 4),
    guide = guide_legend(order = 2, override.aes = list(fill = "#B8C7CC", colour = "#FFFFFF", alpha = 0.9))
  ) +
  scale_fill_gradient(
    name = "-log10(FDR)",
    low = "#DCE8EB",
    high = "#9B2F33",
    guide = guide_colorbar(order = 1, barheight = unit(28, "mm"), barwidth = unit(3.2, "mm"))
  ) +
  labs(y = NULL, tag = "A") +
  theme_classic(base_size = 7.2, base_family = "Arial") +
  theme(
    plot.background = element_rect(fill = "white", colour = NA),
    panel.background = element_rect(fill = "white", colour = NA),
    panel.grid.major.x = element_line(colour = "#E9ECEF", linewidth = 0.26),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line.y = element_blank(),
    axis.line.x = element_line(colour = "#202020", linewidth = 0.32),
    axis.ticks.y = element_blank(),
    axis.ticks.x = element_line(colour = "#202020", linewidth = 0.28),
    axis.text.y = element_text(colour = "#262626", size = 6.7, lineheight = 0.92, margin = margin(r = 5)),
    axis.text.x = element_text(colour = "#262626", size = 6.6),
    axis.title.x = element_text(colour = "#1F1F1F", size = 7.4, margin = margin(t = 5)),
    strip.background = element_blank(),
    strip.text.x = element_text(face = "bold", size = 7.8, colour = "#202020", margin = margin(b = 4)),
    panel.spacing.x = unit(7, "mm"),
    legend.position = "right",
    legend.box = "vertical",
    legend.box.spacing = unit(1.5, "mm"),
    legend.margin = margin(0, 0, 0, 3),
    legend.title = element_text(size = 6.8, colour = "#1F1F1F"),
    legend.text = element_text(size = 6.3, colour = "#262626"),
    legend.key = element_rect(fill = "white", colour = NA),
    legend.background = element_rect(fill = "white", colour = NA),
    plot.tag = element_text(face = "bold", size = 12, colour = "#111111"),
    plot.tag.position = c(0.004, 0.996),
    plot.margin = margin(5.5, 10, 5.5, 5.5)
  )

base <- file.path(out_dir, "Fig3A_MEblue_MEbrown_enrichment")
width_mm <- 178
height_mm <- 116
svglite::svglite(paste0(base, ".svg"), width = width_mm / 25.4, height = height_mm / 25.4)
print(p)
dev.off()
grDevices::cairo_pdf(paste0(base, ".pdf"), width = width_mm / 25.4, height = height_mm / 25.4, family = "Arial")
print(p)
dev.off()
ragg::agg_png(paste0(base, ".png"), width = width_mm, height = height_mm, units = "mm", res = 600, background = "white")
print(p)
dev.off()

module_counts <- data.frame(module = names(targets), gene_n = vapply(targets, length, integer(1)), stringsAsFactors = FALSE)
sel_list <- lapply(split(plot_data, plot_data$module), function(d) paste0("- ", d$module[1], ": ", paste(d$pathway, collapse = " | ")))
map_summary <- if (nrow(all_maps)) {
  paste0("- ", names(table(all_maps$module)), ": ", as.integer(table(all_maps$module)), " merged redundant terms", collapse = "\n")
} else {
  "- No redundant terms merged by the fixed Jaccard thresholds."
}

audit <- c(
  "# Figure 3A MEblue/MEbrown pathway enrichment audit",
  "",
  "Scope: Figure 3A only. No Figure 3B/3C/Figure4, TCGA/WGCNA rerun, or module gene-list modification was performed.",
  paste0("Frozen WGCNA module authority: ", module_path),
  paste0("MEblue gene N: ", module_counts$gene_n[module_counts$module == "MEblue"]),
  paste0("MEbrown gene N: ", module_counts$gene_n[module_counts$module == "MEbrown"]),
  "",
  "## Enrichment pipeline",
  "- Databases: MSigDB human v2025.1.Hs Hallmark, Reactome, and GO Biological Process GMT files.",
  "- Background: all unique SYMBOL genes in WGCNA_module_assignment.csv.",
  "- Method: one-sided hypergeometric ORA; BH-FDR computed independently within each module and database.",
  "- Cutoff for displayed candidates: BH-FDR q < 0.05.",
  "- Redundancy reduction: Hallmark retained as-is; Reactome overlap/Jaccard threshold 0.50; GO-BP overlap/Jaccard threshold 0.35, matching the MEblue enrichment bubble-plot pipeline.",
  "- Selection rule: after redundancy reduction, each module is sorted by BH-FDR q then GeneRatio and the top 8 significant terms are displayed.",
  "",
  "## Final displayed pathways",
  unlist(sel_list, use.names = FALSE),
  "",
  "## Redundancy merging",
  map_summary,
  "",
  "## Outputs",
  paste0("- ", base, ".svg"),
  paste0("- ", base, ".pdf"),
  paste0("- ", base, ".png"),
  paste0("- ", file.path(out_dir, "Fig3A_MEblue_MEbrown_enrichment_source.csv")),
  paste0("- ", file.path(out_dir, "Fig3A_MEblue_MEbrown_enrichment_audit.md"))
)
writeLines(audit, file.path(out_dir, "Fig3A_MEblue_MEbrown_enrichment_audit.md"))

message("Figure 3A MEblue/MEbrown enrichment complete.")
message("Output: ", out_dir)
