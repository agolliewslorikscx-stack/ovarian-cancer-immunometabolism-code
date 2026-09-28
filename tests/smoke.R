# Offline software fixtures only; no manuscript observations or outputs.
files <- c(list.files('analysis', pattern='\\.R$', full.names=TRUE),
           'run_analysis.R', 'check_dependencies.R', 'tests/smoke.R')
for (f in files) {parse(f); cat('PARSE OK:', f, '\n')}
load_functions <- function(file, wanted, env) {
  for (expr in parse(file)) {
    if (is.call(expr) && identical(expr[[1]], as.name('<-')) &&
        is.symbol(expr[[2]]) && as.character(expr[[2]]) %in% wanted) eval(expr, envir=env)
  }
  stopifnot(all(vapply(wanted, exists, logical(1), envir=env, inherits=FALSE)))
}
e <- new.env(parent=globalenv())
load_functions('analysis/module_programme_association.R', 'score_one', e)
m <- rbind(A=c(1,2,3), B=c(10,20,30), C=c(4,4,4))
stopifnot(isTRUE(all.equal(unname(e$score_one(m, c('A','B','C'))), c(-1,0,1))))
stopifnot(all(is.na(e$score_one(m, 'missing'))))
load_functions('analysis/module_enrichment.R', c('jaccard','ora_one'), e)
stopifnot(e$jaccard(c('A','B'),c('B','C')) == 1/3)
g <- data.frame(term_id=c(rep('T1',5),rep('T2',5)), gene=c(LETTERS[1:5],LETTERS[4:8]))
r <- e$ora_one(g, LETTERS[1:3], LETTERS[1:10], 'test', 'test')
stopifnot(nrow(r)==1, r$P==phyper(2,5,5,3,lower.tail=FALSE), r$overlap==3)
load_functions('analysis/donor_cell_composition.R', 'map_group', e)
stopifnot(identical(e$map_group(c('HGSOC','NON_HGSOC_NORMAL')),c('HGSOC','Normal')))
stopifnot(inherits(try(e$map_group('unexpected'),silent=TRUE),'try-error'))
cat('PASS: offline function checks; no research analysis executed.\n')
