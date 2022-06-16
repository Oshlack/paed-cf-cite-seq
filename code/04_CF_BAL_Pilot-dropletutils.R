# Create SingleCellExperiment object from 10X data and identify empty droplets
# using DropletUtils.
# Jovana Maksimovic modified code from Peter Hickey
# 2022-03-07

# Setup ------------------------------------------------------------------------

library(DropletUtils)
library(here)

dir.create(here("data", "SCEs"), recursive = TRUE)
dir.create(here("data", "emptyDrops"))

# Construct SingleCellExperiment object ----------------------------------------
folder <- here()
capture_names <- c("A","B","C","D")
capture_names <- setNames(capture_names, capture_names)
captures <- setNames(
  file.path(
    folder,
    "data/190930_A00152_0150_BHTYCMDSXX/GE",
    capture_names,
    capture_names,
    "outs/multi/count/raw_feature_bc_matrix"),
  capture_names)
sce <- read10xCounts(samples = captures, col.names = TRUE)
stopifnot(!anyDuplicated(colnames(sce)))
sce <- splitAltExps(
  sce,
  rowData(sce)$Type,
  "Gene Expression")

# Identify empty droplets ------------------------------------------------------

set.seed(100)
list_of_empties <- lapply(capture_names, function(cn) {
  message(cn)
  keep <- !(grepl("^MT-", rowData(sce)$Symbol) | grepl("^RP[SL]", rowData(sce)$Symbol))
  emptyDrops(counts(sce)[keep, sce$Sample == cn],
             BPPARAM = BiocParallel::MulticoreParam(workers = BiocParallel::multicoreWorkers() - 1))
})

# Check if more permutations are needed; see
# https://osca.bioconductor.org/quality-control.html#testing-for-empty-droplets
more_permutations_needed <- sapply(list_of_empties, function(e) {
  table(
    Sig = e$FDR <= 0.001,
    Limited = e$Limited)[1, 2] > 0
})
stopifnot(all(!more_permutations_needed))

saveRDS(
  object = sce,
  file = here(
    "data",
    "SCEs",
    "04_CF_BAL_Pilot.CellRanger_v6.SCE.rds"),
  compress = "xz")

for (cn in capture_names) {
  message(cn)
  empties <- list_of_empties[[cn]]
  saveRDS(
    object = empties,
    file = here(
      "data",
      "emptyDrops",
      paste0(cn, ".emptyDrops.rds")),
    compress = "xz")
}
