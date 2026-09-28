# Portable entry point. Run from the repository root with Rscript --vanilla.
args <- commandArgs(trailingOnly = TRUE)
if (!length(args) || args[1] %in% c('--help', '-h')) {
  cat('Usage: Rscript --vanilla run_analysis.R <module-association|core-genes|enrichment|composition> <project-dir> [--check-only]\n')
  cat('       Rscript --vanilla run_analysis.R mr <project-dir> <exposure-id> <outcome-id> [--check-only]\n')
  quit(status = 0)
}
if (length(args) < 2) stop('Supply a project directory. See --help.')
repo <- normalizePath('.', winslash = '/', mustWork = TRUE)
task <- args[1]
check_only <- '--check-only' %in% args
scripts <- c('module-association'='module_programme_association.R',
             'core-genes'='core_gene_prioritization.R', 'enrichment'='module_enrichment.R',
             'composition'='donor_cell_composition.R', 'mr'='mr_pair_worker.R')
if (!task %in% names(scripts)) stop('Unknown workflow. See --help.')
script <- file.path(repo, 'analysis', scripts[[task]])
if (!file.exists(script)) stop('Run this entry point from the repository root.')
if (task == 'mr' && !check_only) dir.create(args[2], recursive=TRUE, showWarnings=FALSE)
project <- normalizePath(args[2], winslash='/', mustWork=TRUE)
expression <- 'data/TCGA_OV/processed/TCGA_OV_log2TPM_primary_tumor_for_score.tsv'
assignment <- 'results/15_TCGA_WGCNA/WGCNA_module_assignment.csv'
eigengenes <- 'results/15_TCGA_WGCNA/WGCNA_module_eigengenes.csv'
definitions <- 'outputs/Figure2A_6stage_myeloid_programmes_20260916/02_programme_definitions_frozen.csv'
scores <- 'outputs/Figure2C_WGCNA_6programme_20260916/01_TCGA_6programme_scores_426x6.csv'
inputs <- switch(task,
  'module-association'=c(expression, assignment, eigengenes, definitions),
  'core-genes'=c(expression, assignment, eigengenes, definitions, scores),
  'enrichment'=c(assignment, paste0('data/gene_sets/', c('h.all.v2025.1.Hs.symbols.gmt', 'c2.cp.reactome.v2025.1.Hs.symbols.gmt', 'c5.go.bp.v2025.1.Hs.symbols.gmt'))),
  'composition'='outputs/FIG4A_SINGLE_CELL_ECOSYSTEM_ATLAS_PUBLICATION_REDRAW_V1/02_PLOT_DATA/FIG4A_umap_plot_ready.csv',
  'mr'=character())
if (length(inputs)) {
  absent <- inputs[!file.exists(file.path(project, inputs))]
  if (length(absent)) stop('Required input files are missing:\n', paste(absent, collapse='\n'), '\nSee docs/INPUTS.md; no input data are bundled.')
}
if (task == 'mr') {
  if (length(args) < 4 || any(startsWith(args[3:4], '--'))) stop('MR requires explicit exposure and outcome IDs.')
  if (!all(grepl('^[A-Za-z0-9_.-]+$', args[3:4]))) stop('Invalid OpenGWAS ID format.')
}
if (check_only) {cat('Input file preflight passed. No analysis executed.\n'); quit(status=0)}
Sys.setenv(OC_PROJECT_ROOT=project)
if (task == 'mr') {
  if (!nzchar(Sys.getenv('OPENGWAS_JWT'))) stop('Set OPENGWAS_JWT in your private environment; never commit it.')
  out <- file.path(project, 'results', 'mr')
  dir.create(out, recursive=TRUE, showWarnings=FALSE)
  pair <- paste(args[3], args[4], sep='__')
  status_file <- file.path(out, paste0(pair, '_status.csv'))
  config <- data.frame(project_dir=project, out_root=out, worker_status=status_file,
    pair_id=pair, exposure_id=args[3], exposure_trait=args[3], exposure_name_zh=args[3],
    outcome_id=args[4], outcome_trait=args[4], outcome_curation_class='user-selected', stringsAsFactors=FALSE)
  tmp <- tempfile(fileext='.csv')
  write.csv(config, tmp, row.names=FALSE, fileEncoding='UTF-8')
  rscript <- file.path(R.home('bin'), if (.Platform$OS.type == 'windows') 'Rscript.exe' else 'Rscript')
  status <- system2(rscript, c('--vanilla', shQuote(script), shQuote(tmp)))
  unlink(tmp)
  if (status != 0) stop('MR worker failed with exit status ', status)
  if (!file.exists(status_file)) stop('MR worker did not create a status record.')
  result <- read.csv(status_file, stringsAsFactors=FALSE)
  if (!all(result$status == 'success')) stop('MR did not complete successfully; inspect its status CSV.')
} else {
  source(script, chdir=FALSE)
}
