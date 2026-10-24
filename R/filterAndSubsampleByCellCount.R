#' Filter donors with few cells and subsample donors with very many cells
#'
#' @param seurat_obj Seurat object.
#' @param donor_col Metadata column identifying donors.
#' @param seed Random seed for reproducibility.
#' @return Filtered Seurat object retaining original metadata and assays.
#' @export
filterAndSubsampleByCellCount <- function(seurat_obj, donor_col = "projid", seed = 0) {
  meta <- seurat_obj[[]]
  if (!donor_col %in% names(meta)) stop("Missing donor column")
  ids <- as.character(meta[[donor_col]])
  if (anyNA(ids)) stop("Missing donor IDs")
  counts <- table(ids)
  lower <- as.numeric(stats::quantile(counts, 0.25))
  upper <- stats::median(counts) + 2 * stats::IQR(counts)
  eligible <- as.integer(counts[counts >= lower & counts <= upper])
  if (!length(eligible)) stop("No donors within cell-count thresholds")
  set.seed(seed)
  keep <- unlist(lapply(names(counts), function(id) {
    cells <- rownames(meta)[ids == id]
    n <- length(cells)
    if (n < lower) return(character(0))
    if (n > upper) return(sample(cells, min(n, sample(eligible, 1))))
    cells
  }), use.names = FALSE)
  subset(seurat_obj, cells = keep)
}
