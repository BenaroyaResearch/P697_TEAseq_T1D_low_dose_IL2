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

# gitcreds: not used by any analysis code -- it is what lets renv authenticate to
#   the PRIVATE BenaroyaResearch/theTools repo when restoring. Without it,
#   renv::restore() fails on theTools with a bare "error code 56" download error.
#   It is on CRAN, so renv installs it before it needs to reach GitHub and the
#   restore becomes self-sufficient.
library(gitcreds)
