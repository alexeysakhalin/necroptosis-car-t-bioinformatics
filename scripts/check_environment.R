expected <- c(Seurat="5.0.1", SeuratObject="5.0.1", Matrix="1.6-5",
              readxl="1.4.3", jsonlite="1.8.8", AnnotationDbi="1.64.1",
              org.Hs.eg.db="3.18.0", GO.db="3.18.0", fgsea="1.32.2",
              BiocParallel="1.36.0")
failures <- character()
if (as.character(getRversion()) != "4.3.3") {
  failures <- c(failures, paste("R", getRversion(), "does not match R 4.3.3"))
}
for (package in names(expected)) {
  actual <- if (requireNamespace(package, quietly=TRUE)) as.character(packageVersion(package)) else "missing"
  if (actual != as.character(package_version(expected[[package]]))) {
    failures <- c(failures, paste(package, actual, "does not match", expected[[package]]))
  }
}
if (length(failures)) stop(paste(c("Reference environment mismatch", failures), collapse="\n"))
cat("R and the ten required package versions match the reference environment.\n")
