#' Softmax distance-weighted CPP scores for cells
#'
#' @param milo_obj Milo object with neighborhood membership.
#' @param neighborhoodPhenotypeAssociations Named list of testNhoods tables, each with logFC.
#' @param cellNeighborhoodDistanceMatrix Sparse cells x neighborhoods diffusion distances.
#' @param col_name Cell metadata identifier column; NULL uses Milo cell names.
#' @param varsToTest Phenotypes to score; defaults to list names.
#' @param temperature Positive softmax temperature (notebook default: 0.3).
#' @return Data frame with one row per cell and one column per phenotype.
#' @export
calculateCellScores <- function(milo_obj, neighborhoodPhenotypeAssociations, cellNeighborhoodDistanceMatrix,
                                col_name = NULL, varsToTest = names(neighborhoodPhenotypeAssociations), temperature = 0.3) {
  if (temperature <= 0) stop("temperature must be positive")
  M <- as(milo_obj@nhoods, "CsparseMatrix")
  D <- as(cellNeighborhoodDistanceMatrix, "CsparseMatrix")
  if (!identical(dim(M), dim(D)) || !identical(M@p, D@p) || !identical(M@i, D@i))
    stop("Distance matrix must have the same sparse structure as Milo neighborhoods")
  if (!all(varsToTest %in% names(neighborhoodPhenotypeAssociations))) stop("Missing phenotype associations")
  score_one <- function(result) {
    logfc <- as.numeric(result$logFC)
    if (length(logfc) != ncol(M)) stop("logFC must be ordered by Milo neighborhoods")
    if (anyNA(logfc)) stop("Missing neighborhood logFC")
    out <- rep(NA_real_, nrow(M))
    A <- M
    A@x <- logfc[rep(seq_len(ncol(M)), diff(M@p))]
    At <- as(Matrix::t(A), "CsparseMatrix")
    Dt <- as(Matrix::t(D), "CsparseMatrix")
    for (cell in seq_len(ncol(At))) {
      start <- At@p[cell] + 1L; end <- At@p[cell + 1L]
      if (start > end) next
      ix <- start:end
      z <- -Dt@x[ix] / temperature
      z <- z - max(z)
      w <- exp(z); w <- w / sum(w)
      out[cell] <- sum(At@x[ix] * w)
    }
    out
  }
  values <- lapply(varsToTest, function(v) score_one(neighborhoodPhenotypeAssociations[[v]]))
  scores <- as.data.frame(setNames(values, varsToTest), check.names = FALSE)
  rownames(scores) <- rownames(M)
  if (!is.null(col_name)) {
    meta <- as.data.frame(SummarizedExperiment::colData(milo_obj))
    if (!col_name %in% colnames(meta)) stop("col_name not found in Milo metadata")
    scores[[col_name]] <- meta[rownames(M), col_name]
  }
  scores
}
