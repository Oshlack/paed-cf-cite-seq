# Create SingleCellExperiment object from 10X data and identify empty droplets
# using DropletUtils.
# Peter Hickey
# 2021-05-10

# Setup ------------------------------------------------------------------------

library(DropletUtils)
library(here)

dir.create(here("data", "SCEs"), recursive = TRUE)
dir.create(here("data", "emptyDrops"))

# Construct SingleCellExperiment object ----------------------------------------

capture_names <- paste0("C133_", 1:2)
capture_names <- setNames(capture_names, capture_names)
captures <- setNames(
  here(
    "extdata",
    "210312_A01221_0027_BHJ7TKDSXY",
    "CellRanger",
    capture_names,
    "outs",
    "raw_feature_bc_matrix"),
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
  emptyDrops(counts(sce)[, sce$Sample == cn])
})

# Check if more permutations are needed; see
# https://osca.bioconductor.org/quality-control.html#testing-for-empty-droplets
more_permutations_needed <- sapply(list_of_empties, function(e) {
  table(
    Sig = e$FDR <= 0.001,
    Limited = e$Limited)[1, 2] > 0
})
stopifnot(all(!more_permutations_needed))

# Compare `lower` cutoffs ------------------------------------------------------

set.seed(100)
lower <- lower <- c(100, 50, 20)
names(lower) <- lower
tmp <- lapply(c(100, 50, 20), function(lower) {
  message(lower)
  lapply(capture_names, function(cn) {
    message(cn)
    emptyDrops(counts(sce)[, sce$Sample == cn], lower = lower)
  })
})
tmp2 <- data.frame(
  C133_1 = sapply(
    lapply(tmp, "[[", "C133_1"), function(x) {
      sum(x$FDR < 0.001, na.rm = TRUE)
    }),
  C133_2 = sapply(
    lapply(tmp, "[[", "C133_2"), function(x) {
      sum(x$FDR < 0.001, na.rm = TRUE)
    }))
rownames(tmp2) <- paste0("emptyDrops(..., lower = ", lower, ")")
cr <- read10xCounts(sub("raw", "filtered", captures))
tmp2 <- rbind(
  tmp2,
  data.frame(
    C133_1 = sum(cr$Sample == "C133_1"),
    C133_2 = sum(cr$Sample == "C133_2"),
    row.names = "CellRanger"))

# Save outputs -----------------------------------------------------------------

saveRDS(
  object = sce,
  file = here(
    "data",
    "SCEs",
    "C133_Neeland.CellRanger.SCE.rds"),
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
  writeLines(
    text = sce[["Barcode"]][sce$Sample == cn][which(empties$FDR <= 0.001)],
    con = here(
      "data",
      "emptyDrops",
      paste0(cn, ".barcodes.txt")))
}

saveRDS(tmp2, file = here("data/emptyDrops/varying_lower.emptyDrops.rds"))
