
#' Compute Diffusion CPP Coordinates
#'
#' Computes diffusion coordinates using cosine kNN, Gaussian connectivity,
#' and diffusion-map eigenvectors scaled by eigenvalues raised to diffusion time.
#'
#' @param adata AnnData object containing PCA coordinates.
#' @param pca_key Name of the PCA representation.
#' @param k Number of nearest neighbors.
#' @param n_components Number of diffusion components.
#' @param diffusion_time Diffusion time.
#' @param random_state Random seed.
#'
#' @return Matrix of diffusion coordinates in the original cell order.
#' @export

computeDiffusionCoordinates <- function(adata, pca_key = "X_pca", k = 30, n_components = 12, diffusion_time = 3, random_state = 0) {

  sc <- reticulate::import("scanpy", convert = FALSE)
  np <- reticulate::import("numpy", convert = FALSE)

  tmp <- sc$AnnData(X = np$zeros(reticulate::tuple(adata$n_obs, 1L), dtype = np$float32))
  tmp$obs_names <- adata$obs_names$copy()
  tmp$obsm[[pca_key]] <- np$asarray(adata$obsm[[pca_key]], dtype = np$float32)

  sc$pp$neighbors(tmp, n_neighbors = as.integer(k), use_rep = pca_key, metric = "cosine", method = "gauss", knn = TRUE, random_state = as.integer(random_state))
  sc$tl$diffmap(tmp, n_comps = as.integer(n_components + 1), random_state = as.integer(random_state))

  psi <- reticulate::py_to_r(tmp$obsm[["X_diffmap"]])[, 2:(n_components + 1), drop = FALSE]
  eigenvalues <- reticulate::py_to_r(tmp$uns[["diffmap_evals"]])[2:(n_components + 1)]
  eigenvalues <- pmin(pmax(eigenvalues, 0), 1)

  if (ncol(psi) != length(eigenvalues)) stop("Too few diffusion components")

  coordinates <- sweep(psi, 2, eigenvalues^diffusion_time, FUN = "*")
  rownames(coordinates) <- reticulate::py_to_r(adata$obs_names$to_list())

  return(coordinates)
}
