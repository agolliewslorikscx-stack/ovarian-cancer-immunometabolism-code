pkgs <- c('ieugwasr','TwoSampleMR','httr','ggplot2','svglite','ragg',
          'readr','dplyr','tidyr','stringr','gridExtra','scales')
versions <- vapply(pkgs, function(p) {
  if (requireNamespace(p, quietly=TRUE)) as.character(packageVersion(p)) else 'NOT INSTALLED'
}, character(1))
cat(R.version.string, '\n')
print(data.frame(package=pkgs, version=unname(versions)), row.names=FALSE)
cat('Install only dependencies for the workflow you intend to run. See docs/DEPENDENCIES.md.\n')
