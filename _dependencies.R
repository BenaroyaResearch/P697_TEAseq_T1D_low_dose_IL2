# renv implicit-dependency manifest.
#
# Packages listed here are used by the project but are never referenced via
# library()/require()/`::` in analysis code, so renv's static dependency scan
# cannot see them. Declaring them here keeps renv::snapshot() from silently
# dropping them from renv.lock.
#
# presto: used internally by Seurat::FindMarkers() (fast Wilcoxon / AUC) when
#   installed -- nothing in this project calls presto:: directly.
library(presto)
