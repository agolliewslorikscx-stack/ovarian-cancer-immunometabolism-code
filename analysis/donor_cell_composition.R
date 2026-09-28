#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(grid)
  library(svglite)
  library(ragg)
})

project <- normalizePath(Sys.getenv("OC_PROJECT_ROOT", unset = "."), winslash = "/", mustWork = TRUE)
out_dir <- file.path(project, "outputs", "Figure4_cell_composition_violin_20260917")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

meta_path <- file.path(project, "outputs", "FIG4A_SINGLE_CELL_ECOSYSTEM_ATLAS_PUBLICATION_REDRAW_V1",
                       "02_PLOT_DATA", "FIG4A_umap_plot_ready.csv")
stopifnot(file.exists(meta_path))

cell_order <- c(
  "Malignant",
  "Normal epithelial",
  "Myeloid",
  "T/NK",
  "B/plasma",
  "Stromal",
  "Endothelial",
  "unresolved"
)

map_group <- function(x) {
  out <- ifelse(x == "HGSOC", "HGSOC",
                ifelse(x == "NON_HGSOC_NORMAL", "Normal", as.character(x)))
  if (any(!out %in% c("Normal", "HGSOC"))) {
    stop("Unexpected disease_group values: ", paste(unique(out[!out %in% c("Normal", "HGSOC")]), collapse = ", "))
  }
  out
}

map_identity <- function(primary, display) {
  x <- as.character(display)
  p <- as.character(primary)
  out <- ifelse(x == "Fibroblast/stromal" | p == "Stromal", "Stromal",
         ifelse(x == "Unresolved" | p == "Other_unresolved", "unresolved",
         ifelse(x %in% cell_order, x, x)))
  if (any(!out %in% cell_order)) {
    stop("Unexpected major cell identity values after mapping: ", paste(unique(out[!out %in% cell_order]), collapse = ", "))
  }
  out
}

fmt_p <- function(x) {
  ifelse(is.na(x), "NA", ifelse(x < 0.001, "<0.001", sprintf("%.3f", x)))
}

meta <- read.csv(meta_path, stringsAsFactors = FALSE, check.names = FALSE)
required <- c("cell_id", "patient_id", "disease_group", "major_compartment_primary", "display_label")
missing <- setdiff(required, names(meta))
if (length(missing) > 0) stop("Missing metadata columns: ", paste(missing, collapse = ", "))

meta$group <- map_group(meta$disease_group)
meta$cell_identity <- map_identity(meta$major_compartment_primary, meta$display_label)

sample_group <- unique(meta[, c("patient_id", "group")])
if (any(duplicated(sample_group$patient_id))) stop("A patient_id maps to multiple groups.")

totals <- as.data.frame(table(meta$patient_id), stringsAsFactors = FALSE)
names(totals) <- c("patient_id", "total_cells")

counts <- as.data.frame(table(meta$patient_id, meta$cell_identity), stringsAsFactors = FALSE)
names(counts) <- c("patient_id", "cell_identity", "cell_count")

all_grid <- expand.grid(patient_id = sample_group$patient_id, cell_identity = cell_order, stringsAsFactors = FALSE)
comp <- merge(all_grid, counts, by = c("patient_id", "cell_identity"), all.x = TRUE)
comp$cell_count[is.na(comp$cell_count)] <- 0L
comp <- merge(comp, totals, by = "patient_id", all.x = TRUE)
comp <- merge(comp, sample_group, by = "patient_id", all.x = TRUE)
comp$composition_percent <- 100 * comp$cell_count / comp$total_cells
comp$group <- factor(comp$group, levels = c("Normal", "HGSOC"))
comp$cell_identity <- factor(comp$cell_identity, levels = cell_order)
comp <- comp[order(comp$cell_identity, comp$group, comp$patient_id), ]

stats <- do.call(rbind, lapply(cell_order, function(ct) {
  d <- comp[comp$cell_identity == ct, ]
  normal <- d$composition_percent[d$group == "Normal"]
  hgsoc <- d$composition_percent[d$group == "HGSOC"]
  wt <- suppressWarnings(wilcox.test(hgsoc, normal, alternative = "two.sided", exact = FALSE))
  med_normal <- median(normal, na.rm = TRUE)
  med_hgsoc <- median(hgsoc, na.rm = TRUE)
  data.frame(
    cell_identity = ct,
    test_method = "Wilcoxon rank-sum test (Mann-Whitney), two-sided, normal approximation",
    P = wt$p.value,
    P_type = "approximate",
    n_normal = length(normal),
    n_hgsoc = length(hgsoc),
    median_normal = med_normal,
    median_hgsoc = med_hgsoc,
    hgsoc_minus_normal_median = med_hgsoc - med_normal,
    direction = ifelse(med_hgsoc > med_normal, "increased_in_HGSOC",
                       ifelse(med_hgsoc < med_normal, "decreased_in_HGSOC", "no_median_difference")),
    stringsAsFactors = FALSE
  )
}))
stats$BH_q <- p.adjust(stats$P, method = "BH")
stats <- stats[match(cell_order, stats$cell_identity), ]

comp_out <- merge(comp, stats[, c("cell_identity", "test_method", "P", "BH_q", "n_normal", "n_hgsoc")],
                  by = "cell_identity", all.x = TRUE)
comp_out <- comp_out[order(factor(comp_out$cell_identity, levels = cell_order), comp_out$group, comp_out$patient_id), ]

write.csv(comp_out, file.path(out_dir, "Fig4_cell_composition_violin_source.csv"), row.names = FALSE)
write.csv(stats, file.path(out_dir, "Fig4_cell_composition_violin_statistics.csv"), row.names = FALSE)

ann <- stats
ann$cell_identity <- factor(ann$cell_identity, levels = cell_order)
yr <- do.call(rbind, lapply(cell_order, function(ct) {
  y <- comp$composition_percent[comp$cell_identity == ct]
  ymax <- max(y, na.rm = TRUE)
  pad <- max(2, ymax * 0.12)
  data.frame(cell_identity = ct, y_line = ymax + pad * 0.45, y_text = ymax + pad * 0.80, y_max = ymax + pad * 1.05)
}))
yr$cell_identity <- factor(yr$cell_identity, levels = cell_order)
ann <- merge(ann, yr, by = "cell_identity", all.x = TRUE)
ann$label <- paste0("q = ", fmt_p(ann$BH_q))

pal <- c(Normal = "#4F9A9A", HGSOC = "#C56A5D")

p <- ggplot(comp, aes(x = group, y = composition_percent, fill = group, colour = group)) +
  geom_violin(width = 0.78, alpha = 0.22, linewidth = 0.30, trim = TRUE, scale = "width", na.rm = TRUE) +
  geom_boxplot(width = 0.22, alpha = 0.78, outlier.shape = NA, linewidth = 0.32, colour = "#333333", na.rm = TRUE) +
  geom_point(
    position = position_jitter(width = 0.075, height = 0, seed = 42),
    size = 1.45, alpha = 0.88, stroke = 0.20, shape = 21, colour = "#333333", na.rm = TRUE
  ) +
  geom_segment(data = ann, aes(x = 1, xend = 2, y = y_line, yend = y_line),
               inherit.aes = FALSE, linewidth = 0.28, colour = "#2B2B2B") +
  geom_segment(data = ann, aes(x = 1, xend = 1, y = y_line, yend = y_line * 0.985),
               inherit.aes = FALSE, linewidth = 0.28, colour = "#2B2B2B") +
  geom_segment(data = ann, aes(x = 2, xend = 2, y = y_line, yend = y_line * 0.985),
               inherit.aes = FALSE, linewidth = 0.28, colour = "#2B2B2B") +
  geom_text(data = ann, aes(x = 1.5, y = y_text, label = label),
            inherit.aes = FALSE, size = 2.25, colour = "#222222") +
  facet_wrap(~ cell_identity, ncol = 4, scales = "free_y") +
  scale_fill_manual(values = pal, name = NULL) +
  scale_colour_manual(values = pal, guide = "none") +
  labs(x = NULL, y = "Cell composition (%)", tag = "B") +
  theme_classic(base_size = 7.2, base_family = "Arial") +
  theme(
    plot.background = element_rect(fill = "white", colour = NA),
    panel.background = element_rect(fill = "white", colour = NA),
    panel.grid.major.y = element_line(colour = "#ECEFF1", linewidth = 0.22),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(linewidth = 0.30, colour = "#222222"),
    axis.ticks = element_line(linewidth = 0.25, colour = "#222222"),
    axis.text.x = element_text(size = 6.4, colour = "#222222", angle = 25, hjust = 1),
    axis.text.y = element_text(size = 6.2, colour = "#222222"),
    axis.title.y = element_text(size = 7.4, margin = margin(r = 5)),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 7.2, colour = "#202020"),
    legend.position = "bottom",
    legend.text = element_text(size = 6.7),
    legend.key.width = unit(5, "mm"),
    panel.spacing.x = unit(4.5, "mm"),
    panel.spacing.y = unit(5.2, "mm"),
    plot.tag = element_text(face = "bold", size = 12, colour = "#111111"),
    plot.tag.position = c(0.004, 0.996),
    plot.margin = margin(5.5, 5.5, 5.5, 5.5)
  )

base <- file.path(out_dir, "Fig4_cell_composition_violin")
width_mm <- 178
height_mm <- 128
svglite::svglite(paste0(base, ".svg"), width = width_mm / 25.4, height = height_mm / 25.4)
print(p)
dev.off()
grDevices::cairo_pdf(paste0(base, ".pdf"), width = width_mm / 25.4, height = height_mm / 25.4, family = "Arial")
print(p)
dev.off()
ragg::agg_png(paste0(base, ".png"), width = width_mm, height = height_mm, units = "mm", res = 600, background = "white")
print(p)
dev.off()

sample_counts <- as.data.frame(table(sample_group$group), stringsAsFactors = FALSE)
names(sample_counts) <- c("group", "patient_sample_n")

audit <- c(
  "# Figure 4 cell composition violin audit",
  "",
  "Scope: cell-composition violin only. Cell annotations, patient labels, Normal/HGSOC grouping, UMAP, Figure 3, myeloid substates, and S01/S04 analyses were not modified.",
  paste0("Frozen cell metadata: ", meta_path),
  "",
  "## Input summary",
  paste0("- Cells: ", nrow(meta)),
  paste0("- Patient/sample N: ", nrow(sample_group)),
  paste0("- Normal patient/sample N: ", sample_counts$patient_sample_n[sample_counts$group == "Normal"]),
  paste0("- HGSOC patient/sample N: ", sample_counts$patient_sample_n[sample_counts$group == "HGSOC"]),
  "- Fixed cell identities: Malignant, Normal epithelial, Myeloid, T/NK, B/plasma, Stromal, Endothelial, unresolved.",
  "- Mapping: NON_HGSOC_NORMAL -> Normal; Fibroblast/stromal -> Stromal; Unresolved/Other_unresolved -> unresolved.",
  "",
  "## Statistics",
  "- Unit of analysis: patient/sample.",
  "- Test: Wilcoxon rank-sum test, two-sided, normal approximation.",
  "- Multiple testing: BH-FDR across the 8 fixed cell identities.",
  "",
  "## Outputs",
  paste0("- ", base, ".svg"),
  paste0("- ", base, ".pdf"),
  paste0("- ", base, ".png"),
  paste0("- ", file.path(out_dir, "Fig4_cell_composition_violin_source.csv")),
  paste0("- ", file.path(out_dir, "Fig4_cell_composition_violin_statistics.csv")),
  paste0("- ", file.path(out_dir, "Fig4_cell_composition_violin_audit.md"))
)
writeLines(audit, file.path(out_dir, "Fig4_cell_composition_violin_audit.md"))

message("Figure 4 cell composition violin complete.")
message("Output: ", out_dir)
