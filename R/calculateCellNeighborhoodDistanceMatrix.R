#' Sparse cell-to-neighborhood distances in diffusion space
#'
#' Only membership pairs are evaluated; this avoids a dense cells-by-neighborhood matrix.
#' @param NeighborhoodCenters Matrix of neighborhood centers in Milo neighborhood order.
#' @param diffusion_coords Named matrix of cell diffusion coordinates.
#' @param milo_obj Milo object containing @nhoods.
#' @param d Number of coordinates.
#' @param chunk_size Number of membership pairs per chunk.
#' @return Sparse cells x neighborhoods distance matrix, same sparsity as @nhoods.
#' @export
calculateCellNeighborhoodDistanceMatrix <- function(NeighborhoodCenters, diffusion_coords, milo_obj, d = 12, chunk_size = 500000L) {
  M <- as(milo_obj@nhoods, "CsparseMatrix")
  if (!all(rownames(M) %in% rownames(diffusion_coords))) stop("Unmatched cell names")
  if (nrow(NeighborhoodCenters) != ncol(M)) stop("Centers must follow Milo neighborhood order")
  if (ncol(diffusion_coords) < d || ncol(NeighborhoodCenters) < d) stop("Insufficient dimensions")
  if (chunk_size < 1) stop("chunk_size must be positive")
  X <- as.matrix(diffusion_coords[rownames(M), seq_len(d), drop = FALSE])
  C <- as.matrix(NeighborhoodCenters[, seq_len(d), drop = FALSE])
  rows <- M@i + 1L
  cols <- rep(seq_len(ncol(M)), diff(M@p))
  vals <- numeric(length(rows))
  if (length(rows)) for (s in seq.int(1L, length(rows), by = chunk_size)) {
    ix <- s:min(s + chunk_size - 1L, length(rows))
    delta <- X[rows[ix], , drop = FALSE] - C[cols[ix], , drop = FALSE]
    vals[ix] <- sqrt(rowSums(delta^2))
  }
  D <- M
  D@x <- vals
  D
}
