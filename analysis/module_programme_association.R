suppressPackageStartupMessages({
  library(ggplot2)
  library(svglite)
  library(ragg)
})

ROOT <- normalizePath(Sys.getenv("OC_PROJECT_ROOT", unset = "."), winslash = "/", mustWork = TRUE)
OUT <- file.path(ROOT, "outputs/Figure2C_WGCNA_6programme_20260916")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

defs_file <- file.path(ROOT, "outputs/Figure2A_6stage_myeloid_programmes_20260916/02_programme_definitions_frozen.csv")
expr_file <- file.path(ROOT, "data/TCGA_OV/processed/TCGA_OV_log2TPM_primary_tumor_for_score.tsv")
me_file <- file.path(ROOT, "results/15_TCGA_WGCNA/WGCNA_module_eigengenes.csv")
assign_file <- file.path(ROOT, "results/15_TCGA_WGCNA/WGCNA_module_assignment.csv")

defs <- read.csv(defs_file, stringsAsFactors = FALSE, check.names = FALSE)
programmes <- c("Myelopoietic Signalling", "Inflammatory Priming", "Chemokine Recruitment",
                "Endothelial Transmigration", "Myeloid State Differentiation", "Tumour Myeloid Remodelling")
if (!all(programmes %in% defs$programme)) stop("HARD STOP: missing one or more fixed programmes in Figure2A definitions.")
if (!all(grepl("FROZEN|PASS", defs$freeze_status[match(programmes, defs$programme)], ignore.case = TRUE))) {
  stop("HARD STOP: one or more Figure2A programme definitions are not FROZEN/PASS.")
}
program_genes <- setNames(lapply(programmes, function(p) {
  unique(strsplit(defs$exact_gene_list[defs$programme == p][1], ";", fixed = TRUE)[[1]])
}), programmes)
program_ids <- setNames(paste0("P", seq_along(programmes)), programmes)

message("Reading WGCNA eigengenes...")
mes <- read.csv(me_file, stringsAsFactors = FALSE, check.names = FALSE)
module_cols <- setdiff(colnames(mes), c("sample_barcode", "MEgrey"))
if (length(module_cols) != 13) stop(sprintf("HARD STOP: expected 13 non-grey modules, found %d.", length(module_cols)))
rownames(mes) <- mes$sample_barcode
me_mat <- as.matrix(mes[, module_cols, drop = FALSE])
storage.mode(me_mat) <- "double"

message("Reading WGCNA assignment...")
assign <- read.csv(assign_file, stringsAsFactors = FALSE, check.names = FALSE)
gene_col <- if ("gene_symbol" %in% names(assign)) "gene_symbol" else names(assign)[1]
module_col <- if ("module_eigengene" %in% names(assign)) "module_eigengene" else if ("module" %in% names(assign)) "module" else names(assign)[2]
assign[[module_col]] <- ifelse(grepl("^ME", assign[[module_col]]), assign[[module_col]], paste0("ME", assign[[module_col]]))
module_genes <- setNames(lapply(module_cols, function(m) unique(assign[[gene_col]][assign[[module_col]] == m])), module_cols)

message("Reading TCGA expression matrix...")
expr <- read.delim(expr_file, check.names = FALSE, stringsAsFactors = FALSE)
gene_col_expr <- colnames(expr)[1]
genes <- expr[[gene_col_expr]]
expr_mat <- as.matrix(expr[, -1, drop = FALSE])
storage.mode(expr_mat) <- "double"
rownames(expr_mat) <- make.unique(genes)
colnames(expr_mat) <- colnames(expr)[-1]
rm(expr); gc()

common_samples <- intersect(rownames(me_mat), colnames(expr_mat))
if (length(common_samples) != 426) stop(sprintf("HARD STOP: expected 426 matched TCGA samples, found %d.", length(common_samples)))
me_mat <- me_mat[common_samples, , drop = FALSE]
expr_mat <- expr_mat[, common_samples, drop = FALSE]

score_one <- function(mat, genes) {
  present <- intersect(genes, rownames(mat))
  if (!length(present)) return(rep(NA_real_, ncol(mat)))
  sub <- mat[present, , drop = FALSE]
  sds <- apply(sub, 1, sd, na.rm = TRUE)
  keep <- is.finite(sds) & sds > 0
  sub <- sub[keep, , drop = FALSE]
  if (!nrow(sub)) return(rep(NA_real_, ncol(mat)))
  z <- t(scale(t(sub)))
  colMeans(z, na.rm = TRUE)
}

coverage <- do.call(rbind, lapply(programmes, function(p) {
  frozen <- length(program_genes[[p]])
  found <- length(intersect(program_genes[[p]], rownames(expr_mat)))
  data.frame(programme_id = program_ids[[p]], programme = p, frozen_n = frozen,
             found_n = found, missing_n = frozen - found, coverage = found / frozen,
             status = ifelse(found / frozen >= 0.50, "ESTIMABLE", "NOT_ESTIMABLE"),
             missing_genes = paste(setdiff(program_genes[[p]], rownames(expr_mat)), collapse = ";"),
             stringsAsFactors = FALSE)
}))
write.csv(coverage, file.path(OUT, "02_programme_coverage.csv"), row.names = FALSE)
if (any(coverage$status == "NOT_ESTIMABLE")) stop("HARD STOP: one or more programmes has coverage <50%.")

message("Scoring six frozen programmes...")
score_mat <- sapply(programmes, function(p) score_one(expr_mat, program_genes[[p]]))
score_df <- data.frame(sample_barcode = common_samples, score_mat, check.names = FALSE)
write.csv(score_df, file.path(OUT, "01_TCGA_6programme_scores_426x6.csv"), row.names = FALSE)

cor_test <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  ct <- suppressWarnings(cor.test(x[ok], y[ok], method = "pearson"))
  c(N = sum(ok), r = unname(ct$estimate), P = ct$p.value)
}

message("Primary module-programme correlations...")
primary <- list()
for (m in module_cols) {
  for (p in programmes) {
    v <- cor_test(me_mat[, m], score_mat[, p])
    primary[[length(primary) + 1]] <- data.frame(module = m, programme_id = program_ids[[p]], programme = p,
                                                 N = v["N"], r = v["r"], P = v["P"], stringsAsFactors = FALSE)
  }
}
primary <- do.call(rbind, primary)
primary$q78 <- p.adjust(primary$P, method = "BH")
write.csv(primary, file.path(OUT, "03_module_programme_primary_78.csv"), row.names = FALSE)

message("Overlap-excluded sensitivity...")
overlap_rows <- list()
sens_rows <- list()
for (m in module_cols) {
  mg <- unique(module_genes[[m]])
  mg <- mg[!is.na(mg) & nzchar(mg)]
  for (p in programmes) {
    pg <- program_genes[[p]]
    overlap <- intersect(mg, pg)
    union_n <- length(union(mg, pg))
    jacc <- ifelse(union_n > 0, length(overlap) / union_n, NA_real_)
    primary_row <- primary[primary$module == m & primary$programme == p, ]
    overlap_rows[[length(overlap_rows) + 1]] <- data.frame(
      module = m, programme_id = program_ids[[p]], programme = p,
      module_gene_n = length(mg), programme_gene_n = length(pg),
      overlap_n = length(overlap), Jaccard = jacc,
      overlap_genes = paste(overlap, collapse = ";"),
      primary_r = primary_row$r, primary_P = primary_row$P, primary_q78 = primary_row$q78,
      stringsAsFactors = FALSE
    )
    remaining <- setdiff(pg, mg)
    found_remaining <- intersect(remaining, rownames(expr_mat))
    coverage_remaining <- length(found_remaining) / length(pg)
    if (length(remaining) >= 10 && coverage_remaining >= 0.50 && length(found_remaining) >= 10) {
      score <- score_one(expr_mat, remaining)
      v <- cor_test(me_mat[, m], score)
      status <- "ESTIMABLE"
      rr <- data.frame(module = m, programme_id = program_ids[[p]], programme = p,
                       module_gene_n = length(mg), programme_gene_n = length(pg),
                       overlap_n = length(overlap), Jaccard = jacc,
                       remaining_gene_n = length(remaining), remaining_found_n = length(found_remaining),
                       remaining_coverage = coverage_remaining,
                       primary_r = primary_row$r, primary_P = primary_row$P, primary_q78 = primary_row$q78,
                       overlap_excluded_r = v["r"], overlap_excluded_P = v["P"],
                       overlap_excluded_status = status, stringsAsFactors = FALSE)
    } else {
      rr <- data.frame(module = m, programme_id = program_ids[[p]], programme = p,
                       module_gene_n = length(mg), programme_gene_n = length(pg),
                       overlap_n = length(overlap), Jaccard = jacc,
                       remaining_gene_n = length(remaining), remaining_found_n = length(found_remaining),
                       remaining_coverage = coverage_remaining,
                       primary_r = primary_row$r, primary_P = primary_row$P, primary_q78 = primary_row$q78,
                       overlap_excluded_r = NA_real_, overlap_excluded_P = NA_real_,
                       overlap_excluded_status = "NOT_ESTIMABLE", stringsAsFactors = FALSE)
    }
    sens_rows[[length(sens_rows) + 1]] <- rr
  }
}
overlap_audit <- do.call(rbind, overlap_rows)
sens <- do.call(rbind, sens_rows)
sens$overlap_excluded_q <- NA_real_
idx <- which(is.finite(sens$overlap_excluded_P))
sens$overlap_excluded_q[idx] <- p.adjust(sens$overlap_excluded_P[idx], method = "BH")
sens$direction_preserved <- ifelse(is.finite(sens$overlap_excluded_r),
                                   sign(sens$primary_r) == sign(sens$overlap_excluded_r), NA)
write.csv(overlap_audit, file.path(OUT, "04_module_programme_overlap_audit.csv"), row.names = FALSE)
write.csv(sens, file.path(OUT, "05_module_programme_overlap_excluded.csv"), row.names = FALSE)

top_rows <- do.call(rbind, lapply(programmes, function(p) {
  d <- primary[primary$programme == p, ]
  pos <- d[which.max(d$r), ]
  neg <- d[which.min(d$r), ]
  spos <- sens[sens$module == pos$module & sens$programme == p, ]
  sneg <- sens[sens$module == neg$module & sens$programme == p, ]
  data.frame(programme_id = program_ids[[p]], programme = p,
             strongest_positive_module = pos$module, strongest_positive_r = pos$r,
             strongest_positive_q78 = pos$q78,
             strongest_positive_FDR_lt_0_05 = pos$q78 < 0.05,
             positive_overlap_excluded_r = spos$overlap_excluded_r,
             positive_overlap_excluded_q = spos$overlap_excluded_q,
             positive_direction_preserved = spos$direction_preserved,
             strongest_negative_module = neg$module, strongest_negative_r = neg$r,
             strongest_negative_q78 = neg$q78,
             strongest_negative_FDR_lt_0_05 = neg$q78 < 0.05,
             negative_overlap_excluded_r = sneg$overlap_excluded_r,
             negative_overlap_excluded_q = sneg$overlap_excluded_q,
             negative_direction_preserved = sneg$direction_preserved,
             stringsAsFactors = FALSE)
}))
write.csv(top_rows, file.path(OUT, "06_top_module_per_programme.csv"), row.names = FALSE)

key <- primary[primary$module %in% c("MEblue", "MEbrown", "MEturquoise"), ]
key <- merge(key, sens[, c("module","programme","overlap_n","Jaccard","overlap_excluded_r","overlap_excluded_P","overlap_excluded_q","direction_preserved","overlap_excluded_status")],
             by = c("module", "programme"), all.x = TRUE)
key$programme_id <- program_ids[key$programme]
key <- key[order(factor(key$module, c("MEblue","MEbrown","MEturquoise")), factor(key$programme, programmes)), ]
write.csv(key, file.path(OUT, "07_key_modules_MEblue_MEbrown_MEturquoise.csv"), row.names = FALSE)

plot_df <- primary
plot_df$module <- factor(plot_df$module, levels = rev(module_cols))
plot_df$programme <- factor(plot_df$programme, levels = programmes,
                            labels = paste0(program_ids[programmes], "\n", programmes))
plot_df$text_col <- ifelse(plot_df$q78 < 0.05, "#111111", "#A8A8A8")
plot_df$label <- sprintf("%.2f", plot_df$r)
plot_df$is_col_top <- FALSE
for (p in programmes) {
  d <- primary[primary$programme == p & primary$q78 < 0.05, ]
  if (nrow(d)) {
    mtop <- d$module[which.max(d$r)]
    plot_df$is_col_top[plot_df$programme == paste0(program_ids[[p]], "\n", p) & plot_df$module == mtop] <- TRUE
  }
}

p_heat <- ggplot(plot_df, aes(x = programme, y = module, fill = r)) +
  geom_tile(color = "#FFFFFF", linewidth = 0.6) +
  geom_tile(data = plot_df[plot_df$is_col_top, ], fill = NA, color = "#222222", linewidth = 0.75) +
  geom_text(aes(label = label, color = text_col), size = 2.55, family = "Arial") +
  scale_color_identity() +
  scale_fill_gradient2(low = "#3f78a8", mid = "#f7f7f7", high = "#c56e56",
                       midpoint = 0, limits = c(-1, 1), name = "Pearson r") +
  labs(x = NULL, y = NULL, tag = "C") +
  theme_minimal(base_family = "Arial", base_size = 8) +
  theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 35, hjust = 1, vjust = 1, color = "#222222"),
        axis.text.y = element_text(color = "#222222"), legend.position = "right",
        plot.margin = margin(8, 8, 8, 18), plot.tag = element_text(face = "bold", size = 16, color = "#111111"), plot.tag.position = c(0.01, 0.99))

ggsave(file.path(OUT, "Fig2C_module_6programme_heatmap.svg"), p_heat, width = 183/25.4, height = 105/25.4, device = svglite)
ggsave(file.path(OUT, "Fig2C_module_6programme_heatmap.pdf"), p_heat, width = 183/25.4, height = 105/25.4, device = cairo_pdf)
ragg::agg_png(file.path(OUT, "Fig2C_module_6programme_heatmap.png"), width = 183, height = 105, units = "mm", res = 600, background = "white")
print(p_heat)
dev.off()

serious_overlap <- sens[sens$overlap_n > 0 & is.finite(sens$overlap_excluded_r) &
                          sign(sens$primary_r) != sign(sens$overlap_excluded_r) &
                          sens$primary_q78 < 0.05, ]
log_lines <- c(
  "# Figure 2C WGCNA x six frozen myeloid programmes audit log",
  "",
  paste0("Created: ", Sys.time()),
  paste0("Definitions: ", defs_file),
  paste0("TCGA expression: ", expr_file),
  paste0("WGCNA eigengenes: ", me_file),
  paste0("WGCNA assignment: ", assign_file),
  paste0("Matched samples: ", length(common_samples), " (expected 426)."),
  paste0("Non-grey modules reused: ", paste(module_cols, collapse = "; ")),
  "No WGCNA reconstruction, dynamicTreeCut, merge, core-gene screening, ssGSEA or GSVA was run.",
  "Programme score method: gene-wise z-score across matched TCGA samples, then mean across available frozen programme genes.",
  paste0("Serious overlap-driven FDR-supported direction flips after overlap exclusion: ", nrow(serious_overlap))
)
writeLines(log_lines, file.path(OUT, "audit_log.md"), useBytes = TRUE)
message("DONE Figure2C outputs: ", OUT)


