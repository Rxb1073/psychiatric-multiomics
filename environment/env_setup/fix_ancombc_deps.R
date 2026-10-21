personal <- Sys.getenv("R_LIBS_USER"); if (nchar(personal)==0) personal <- path.expand("~/R/mmuphin_lib")
.libPaths(c(personal, .libPaths()))
if (!requireNamespace("BiocManager", quietly=TRUE))
  install.packages("BiocManager", repos="https://cloud.r-project.org", lib=personal)

# reinstall the tree stack + ANCOMBC coherently, forcing rebuild so versions align
pkgs <- c("tidytree","treeio","ape","TreeSummarizedExperiment","ANCOMBC")
BiocManager::install(pkgs, update=TRUE, ask=FALSE, lib=personal, force=TRUE)

for (p in pkgs) {
  ok <- requireNamespace(p, quietly=TRUE)
  cat(sprintf("%-28s %s\n", p, if (ok) "OK" else "FAILED"))
}
suppressMessages({ library(TreeSummarizedExperiment); library(ANCOMBC) })
cat("BOTH LOADED: TRUE\n")
