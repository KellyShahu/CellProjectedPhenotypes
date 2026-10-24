#' Remove cells with low local mixing (LISI) in diffusion space
#'
#' @param seurat_obj Seurat object with metadata.
#' @param diffusion_coords Named cells x diffusion-dimensions matrix.
#' @param col_name Metadata column defining the groups to mix.
#' @param threshold Minimum LISI value for retention.
#' @return Filtered Seurat object.
#' @export
removePoorlyMixedCells <- function(seurat_obj, diffusion_coords, col_name, threshold = 3) {
  meta <- seurat_obj[[]]
  if (!col_name %in% names(meta)) stop("Missing metadata column")
  if (!all(colnames(seurat_obj) %in% rownames(diffusion_coords))) stop("Missing cell coordinates")
  X <- as.matrix(diffusion_coords[colnames(seurat_obj), , drop = FALSE])
  lisi <- lisi::compute_lisi(X, meta[colnames(seurat_obj), , drop = FALSE], label_colnames = col_name)
  keep <- colnames(seurat_obj)[!is.na(lisi[, col_name]) & lisi[, col_name] >= threshold]
  subset(seurat_obj, cells = keep)
}
