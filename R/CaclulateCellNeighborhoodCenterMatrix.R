#' Calculate neighborhood centers in diffusion space
#'
#' @param milo_obj Milo object with neighborhoods already constructed.
#' @param diffusion_coords Named matrix: cells x diffusion coordinates.
#' @param d Number of diffusion coordinates to use.
#' @return Neighborhood x diffusion-coordinate matrix, ordered as milo_obj@nhoods.
#' @export
calculateNeighborhoodCenterMatrix <- function(milo_obj, diffusion_coords, d = 12) {
  M <- as(milo_obj@nhoods, "CsparseMatrix")
  if (is.null(rownames(M)) || is.null(rownames(diffusion_coords))) stop("Cell names are required")
  if (!all(rownames(M) %in% rownames(diffusion_coords))) stop("Missing diffusion coordinates for cells")
  if (ncol(diffusion_coords) < d) stop("Not enough diffusion coordinates")
  X <- as.matrix(diffusion_coords[rownames(M), seq_len(d), drop = FALSE])
  centers <- matrix(NA_real_, nrow = ncol(M), ncol = d)
  for (j in seq_len(ncol(M))) {
    ix <- seq.int(M@p[j] + 1L, M@p[j + 1L])
    if (length(ix) == 0L || M@p[j] == M@p[j + 1L]) stop("Empty neighborhood")
    centers[j, ] <- colMeans(X[M@i[ix] + 1L, , drop = FALSE])
  }
  rownames(centers) <- if (is.null(colnames(M))) as.character(seq_len(ncol(M))) else colnames(M)
  colnames(centers) <- colnames(X)
  centers
}
