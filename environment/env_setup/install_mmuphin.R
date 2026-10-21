personal <- Sys.getenv("R_LIBS_USER")
if (nchar(personal)==0) personal <- path.expand("~/R/mmuphin_lib")
dir.create(personal, recursive=TRUE, showWarnings=FALSE)
.libPaths(c(personal, .libPaths()))
if (!requireNamespace("BiocManager", quietly=TRUE))
  install.packages("BiocManager", repos="https://cloud.r-project.org")
for (pkg in c("MMUPHin","vegan")) {
  if (!requireNamespace(pkg, quietly=TRUE)) BiocManager::install(pkg, update=FALSE, ask=FALSE)
}
cat("MMUPHin available:", requireNamespace("MMUPHin", quietly=TRUE), "\n")
cat("vegan available:", requireNamespace("vegan", quietly=TRUE), "\n")
