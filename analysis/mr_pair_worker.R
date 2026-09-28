options(stringsAsFactors = FALSE, warn = 1, timeout = 1800)

for (loc in c(".UTF-8", "Chinese_China.UTF-8", "English_United States.utf8", "zh_CN.UTF-8")) {
  ok <- suppressWarnings(try(Sys.setlocale("LC_CTYPE", loc), silent = TRUE))
  if (!inherits(ok, "try-error") && nzchar(ok)) break
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1 || !file.exists(args[1])) stop("Missing worker config CSV")
cfg <- read.csv(args[1], check.names = FALSE, stringsAsFactors = FALSE, fileEncoding = "UTF-8-BOM")

PROJECT_DIR <- cfg$project_dir[1]
setwd(PROJECT_DIR)
OUT_ROOT <- cfg$out_root[1]
PAIR_DIR <- file.path(OUT_ROOT, "01_pair_results")
HARM_DIR <- file.path(OUT_ROOT, "02_harmonised")
INST_DIR <- file.path(OUT_ROOT, "03_instruments")
STATUS_FILE <- cfg$worker_status[1]
for (d in c(PAIR_DIR, HARM_DIR, INST_DIR, dirname(STATUS_FILE))) dir.create(d, showWarnings = FALSE, recursive = TRUE)

timestamp <- function() format(Sys.time(), "%Y-%m-%d %H:%M:%S")
clean_msg <- function(x) {
  if (length(x) == 0 || is.null(x) || all(is.na(x))) return("")
  gsub("[\r\n\t]+", " ", as.character(x[1]))
}
safe_id <- function(x) gsub("[^A-Za-z0-9._-]+", "_", x)

status_cols <- c(
  "pair_id", "exposure_id", "exposure_trait", "exposure_name_zh",
  "outcome_id", "outcome_trait", "outcome_curation_class",
  "status", "n_instruments", "n_outcome_snps", "n_harmonised_snps",
  "mr_methods_n", "primary_method", "primary_beta", "primary_se", "primary_p",
  "primary_or", "primary_or_lci95", "primary_or_uci95",
  "error_message", "result_source", "start_time", "end_time"
)

make_row <- function(status, n_instruments = NA, n_outcome_snps = NA, n_harmonised_snps = NA,
                     mr_methods_n = NA, primary_method = NA, primary_beta = NA,
                     primary_se = NA, primary_p = NA, error_message = "",
                     result_source = "OMEGA3_PROTECTIVE_V1_worker", start_time = NA, end_time = timestamp()) {
  primary_or <- primary_or_lci95 <- primary_or_uci95 <- NA_real_
  if (!is.na(suppressWarnings(as.numeric(primary_beta))) && !is.na(suppressWarnings(as.numeric(primary_se)))) {
    b <- as.numeric(primary_beta)
    se <- as.numeric(primary_se)
    primary_or <- exp(b)
    primary_or_lci95 <- exp(b - 1.96 * se)
    primary_or_uci95 <- exp(b + 1.96 * se)
  }
  out <- data.frame(
    pair_id = cfg$pair_id[1],
    exposure_id = cfg$exposure_id[1],
    exposure_trait = cfg$exposure_trait[1],
    exposure_name_zh = cfg$exposure_name_zh[1],
    outcome_id = cfg$outcome_id[1],
    outcome_trait = cfg$outcome_trait[1],
    outcome_curation_class = cfg$outcome_curation_class[1],
    status = status,
    n_instruments = n_instruments,
    n_outcome_snps = n_outcome_snps,
    n_harmonised_snps = n_harmonised_snps,
    mr_methods_n = mr_methods_n,
    primary_method = primary_method,
    primary_beta = primary_beta,
    primary_se = primary_se,
    primary_p = primary_p,
    primary_or = primary_or,
    primary_or_lci95 = primary_or_lci95,
    primary_or_uci95 = primary_or_uci95,
    error_message = clean_msg(error_message),
    result_source = result_source,
    start_time = start_time,
    end_time = end_time,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  out[, status_cols, drop = FALSE]
}

write_status <- function(row) write.csv(row, STATUS_FILE, row.names = FALSE, na = "", fileEncoding = "UTF-8")

with_retry <- function(expr, label, max_retries = 2, wait_seconds = 10) {
  expr_sub <- substitute(expr)
  expr_env <- parent.frame()
  last_error <- NULL
  for (attempt in seq_len(max_retries + 1)) {
    ok_env <- new.env(parent = emptyenv())
    ok_env$ok <- FALSE
    value <- tryCatch({
      z <- eval(expr_sub, envir = expr_env)
      ok_env$ok <- TRUE
      z
    }, error = function(e) {
      last_error <<- e
      NULL
    })
    if (ok_env$ok) return(value)
    if (attempt <= max_retries) {
      cat(timestamp(), "|", label, "failed attempt", attempt, "-", clean_msg(conditionMessage(last_error)), "\n")
      Sys.sleep(wait_seconds)
    }
  }
  stop(label, " failed: ", clean_msg(conditionMessage(last_error)))
}

start_time <- timestamp()
row <- tryCatch({
  suppressPackageStartupMessages({
    library(ieugwasr)
    library(TwoSampleMR)
  })
  if (requireNamespace("httr", quietly = TRUE)) httr::set_config(httr::timeout(1750))
  jwt <- ieugwasr::get_opengwas_jwt()
  if (is.null(jwt) || !nzchar(jwt)) stop("No OpenGWAS JWT found")

  exposure_id <- cfg$exposure_id[1]
  outcome_id <- cfg$outcome_id[1]
  pair_id <- cfg$pair_id[1]

  cached_inst <- file.path(OUT_ROOT, "00_exposure_audit", paste0(safe_id(exposure_id), "_instruments.csv"))
  inst <- if (file.exists(cached_inst)) {
    read.csv(cached_inst, check.names = FALSE, stringsAsFactors = FALSE, fileEncoding = "UTF-8-BOM")
  } else {
    with_retry(TwoSampleMR::extract_instruments(outcomes = exposure_id, p1 = 5e-8, opengwas_jwt = jwt),
               "extract_instruments")
  }
  if (is.null(inst) || nrow(inst) == 0 || !"SNP" %in% names(inst)) {
    make_row("no_instruments", n_instruments = 0, error_message = "extract_instruments returned no SNPs",
             start_time = start_time)
  } else {
    inst_file <- file.path(INST_DIR, paste0(safe_id(pair_id), "_instruments.csv"))
    write.csv(inst, inst_file, row.names = FALSE, na = "", fileEncoding = "UTF-8")
    out <- with_retry(
      TwoSampleMR::extract_outcome_data(
        snps = unique(inst$SNP),
        outcomes = outcome_id,
        proxies = FALSE,
        opengwas_jwt = jwt,
        splitsize = 25,
        proxy_splitsize = 25
      ),
      "extract_outcome_data"
    )
    if (is.null(out) || nrow(out) == 0 || !"SNP" %in% names(out)) {
      make_row("no_outcome_data", n_instruments = length(unique(inst$SNP)), n_outcome_snps = 0,
               error_message = "extract_outcome_data returned no SNPs", start_time = start_time)
    } else {
      harm <- TwoSampleMR::harmonise_data(inst, out)
      if ("mr_keep" %in% names(harm)) harm_keep <- harm[harm$mr_keep, , drop = FALSE] else harm_keep <- harm
      n_harm <- if ("SNP" %in% names(harm_keep)) length(unique(harm_keep$SNP)) else nrow(harm_keep)
      if (nrow(harm_keep) < 1) {
        make_row("too_few_harmonised_snps",
                 n_instruments = length(unique(inst$SNP)),
                 n_outcome_snps = length(unique(out$SNP)),
                 n_harmonised_snps = 0,
                 error_message = "0 harmonised SNPs kept", start_time = start_time)
      } else {
        res <- TwoSampleMR::mr(harm_keep)
        res$pair_id <- pair_id
        res$exposure_id <- exposure_id
        res$exposure_trait_clean <- cfg$exposure_trait[1]
        res$outcome_id_clean <- outcome_id
        res$outcome_trait_clean <- cfg$outcome_trait[1]
        res$outcome_curation_class <- cfg$outcome_curation_class[1]
        if ("b" %in% names(res) && "se" %in% names(res)) {
          res$or <- exp(res$b)
          res$or_lci95 <- exp(res$b - 1.96 * res$se)
          res$or_uci95 <- exp(res$b + 1.96 * res$se)
        }
        write.csv(harm_keep, file.path(HARM_DIR, paste0(safe_id(pair_id), "_harmonised.csv")), row.names = FALSE, na = "", fileEncoding = "UTF-8")
        write.csv(res, file.path(PAIR_DIR, paste0(safe_id(pair_id), "_MR_results.csv")), row.names = FALSE, na = "", fileEncoding = "UTF-8")
        primary_idx <- which(res$method == "Inverse variance weighted")
        if (length(primary_idx) == 0) primary_idx <- which(res$method == "Wald ratio")
        if (length(primary_idx) == 0) primary_idx <- 1
        rr <- res[primary_idx[1], , drop = FALSE]
        make_row("success",
                 n_instruments = length(unique(inst$SNP)),
                 n_outcome_snps = length(unique(out$SNP)),
                 n_harmonised_snps = n_harm,
                 mr_methods_n = nrow(res),
                 primary_method = rr$method[1],
                 primary_beta = rr$b[1],
                 primary_se = rr$se[1],
                 primary_p = rr$pval[1],
                 error_message = "",
                 start_time = start_time)
      }
    }
  }
}, error = function(e) {
  make_row("api_error", error_message = clean_msg(conditionMessage(e)), start_time = start_time)
})

write_status(row)
quit(save = "no", status = 0)
