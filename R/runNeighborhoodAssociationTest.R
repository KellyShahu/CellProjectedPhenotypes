#' Test phenotype associations with Milo neighborhoods
#'
#' @param milo_obj Milo object after countCells().
#' @param donor_metadata Data frame with one row per donor/sample.
#' @param varsToTest Phenotype column names.
#' @param varsToCorrect Covariate column names, or character(0).
#' @param col_name Donor/sample ID column in donor_metadata matching nhoodCounts columns.
#' @return Named list of Milo differential-abundance tables in neighborhood order.
#' @export
runNeighborhoodAssociationTest <- function(milo_obj, donor_metadata, varsToTest, varsToCorrect = character(0), col_name = "subject") {
  ids <- colnames(miloR::nhoodCounts(milo_obj))
  if (is.null(ids)) stop("Neighborhood counts must have sample IDs")
  if (!all(c(col_name, varsToTest, varsToCorrect) %in% names(donor_metadata))) stop("Metadata columns missing")
  if (anyDuplicated(as.character(donor_metadata[[col_name]]))) stop("Donor IDs must be unique")
  idx <- match(ids, as.character(donor_metadata[[col_name]]))
  if (anyNA(idx)) stop("Missing donor metadata")
  design <- donor_metadata[idx, unique(c(col_name, varsToTest, varsToCorrect)), drop = FALSE]
  rownames(design) <- ids
  if (anyNA(design)) stop("Missing phenotype or covariate values")
  out <- lapply(varsToTest, function(phenotype) {
    form <- stats::reformulate(c(varsToCorrect, phenotype))
    res <- miloR::testNhoods(milo_obj, design = form, design.df = design, fdr.weighting = "none")
    # Milo returns rows in neighborhood order; verify before using for cell scores.
    if (nrow(res) != ncol(milo_obj@nhoods)) stop("Association rows do not match neighborhood count")
    res
  })
  names(out) <- varsToTest
  out
}
