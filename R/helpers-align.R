## ----
## align-helpers.R
## Internal helper functions for the alignment system.
## None of these are exported.
## ----

#' Harmonize metabolite space between two CMs
#'
#' Builds a shared metabolite index from two CM objects so their
#' binary matrices can be compared element-wise. Returns
#' pre-expanded matrices in the union metabolite space.
#'
#' @param cm1 A [ConsortiumMetabolism] object (query).
#' @param cm2 A [ConsortiumMetabolism] object (reference).
#'
#' @return A named list with elements:
#'   \describe{
#'     \item{`X`}{Sparse binary matrix for `cm1` in union space.}
#'     \item{`Y`}{Sparse binary matrix for `cm2` in union space.}
#'     \item{`metabolites`}{Character vector of union
#'       metabolites (sorted).}
#'   }
#'
#' @noRd
.harmonizeMetaboliteSpace <- function(cm1, cm2) {
    cli::cli_abort(
        "Not yet implemented (Phase 1).",
        .internal = TRUE
    )
}

#' Expand an assay matrix to a target metabolite space
#'
#' Generic expansion of any assay matrix (not just binary) to a
#' larger metabolite space. Used when comparing weighted assays
#' (e.g., Consumption, Production) between CMs with different
#' metabolite sets.
#'
#' @param mat A sparse or dense matrix with metabolite
#'   row/colnames.
#' @param target_mets Character vector of target metabolite names.
#'
#' @return A matrix of dimensions
#'   `length(target_mets) x length(target_mets)` with the
#'   original values placed at the correct positions and zeros
#'   elsewhere.
#'
#' @noRd
.expandMatrix <- function(mat, target_mets) {
    cli::cli_abort(
        "Not yet implemented (Phase 1).",
        .internal = TRUE
    )
}

#' Compute all similarity scores between two matrices
#'
#' Given pre-expanded binary and weighted matrices for a query and
#' reference, computes all four metrics and returns them as a named
#' list.
#'
#' @param xBin Sparse binary matrix for query (union space).
#' @param yBin Sparse binary matrix for reference (union space).
#' @param xWeighted Named list of weighted assay matrices for
#'   query. Expected names: `"Consumption"`,
#'   `"Production"`.
#' @param yWeighted Named list of weighted assay matrices for
#'   reference. Same structure as `xWeighted`.
#'
#' @return A named list with elements:
#'   \describe{
#'     \item{`FOS`}{Functional Overlap Score (numeric).}
#'     \item{`jaccard`}{Jaccard index (numeric).}
#'     \item{`brayCurtis`}{Bray-Curtis similarity (numeric).}
#'     \item{`redundancyOverlap`}{Redundancy overlap (numeric).}
#'   }
#'
#' @noRd
.computeAllScores <- function(xBin, yBin,
                              xWeighted = NULL,
                              yWeighted = NULL) {
    cli::cli_abort(
        "Not yet implemented (Phase 1).",
        .internal = TRUE
    )
}

#' Functional Overlap Score (Szymkiewicz-Simpson)
#'
#' Computes FOS = |X AND Y| / min(|X|, |Y|), the
#' Szymkiewicz-Simpson overlap coefficient applied to binary
#' pathway matrices. Equivalent to the existing `.binMatOverlap()`
#' but operates on pre-expanded matrices in the same space.
#'
#' @param xBin Sparse binary matrix (query, union space).
#' @param yBin Sparse binary matrix (reference, union space).
#'
#' @return Numeric scalar in \[0, 1\].
#'
#' @noRd
.functionalOverlap <- function(xBin, yBin) {
    cli::cli_abort(
        "Not yet implemented (Phase 1).",
        .internal = TRUE
    )
}

#' Jaccard Index
#'
#' Computes Jaccard = |X AND Y| / |X OR Y| on binary pathway
#' matrices. Symmetric measure of set similarity.
#'
#' @param xBin Sparse binary matrix (query, union space).
#' @param yBin Sparse binary matrix (reference, union space).
#'
#' @return Numeric scalar in \[0, 1\].
#'
#' @noRd
.jaccardIndex <- function(xBin, yBin) {
    cli::cli_abort(
        "Not yet implemented (Phase 1).",
        .internal = TRUE
    )
}

#' Bray-Curtis Similarity
#'
#' Computes Bray-Curtis similarity = 1 - BC dissimilarity on
#' weighted pathway matrices. Uses flux magnitudes from the
#' Consumption and Production assays.
#'
#' @param xWeighted Named list with `Consumption` and
#'   `Production` matrices (query).
#' @param yWeighted Named list with `Consumption` and
#'   `Production` matrices (reference).
#'
#' @return Numeric scalar in \[0, 1\].
#'
#' @noRd
.brayCurtisSimilarity <- function(xWeighted, yWeighted) {
    cli::cli_abort(
        "Not yet implemented (Phase 1).",
        .internal = TRUE
    )
}

#' Redundancy Overlap
#'
#' Computes weighted Jaccard on species-count (nEdges) matrices.
#' Measures how similarly two consortia distribute metabolic
#' labor across species for shared pathways.
#'
#' @param xEdges nEdges assay matrix for query (union space).
#' @param yEdges nEdges assay matrix for reference (union space).
#'
#' @return Numeric scalar in \[0, 1\].
#'
#' @noRd
.redundancyOverlap <- function(xEdges, yEdges) {
    cli::cli_abort(
        "Not yet implemented (Phase 1).",
        .internal = TRUE
    )
}

#' Compute MAAS composite score
#'
#' Metabolic Alignment Aggregate Score. Weighted combination of
#' the four individual metrics. Default weights give FOS 0.4 and
#' each other metric 0.2, but users can supply custom weights.
#'
#' @param scores Named list of individual metric scores (as
#'   returned by `.computeAllScores()`).
#' @param weights Named numeric vector of weights. Must sum to 1.
#'   Names must match elements of `scores`. Default:
#'   `c(FOS = 0.4, jaccard = 0.2, brayCurtis = 0.2,
#'   redundancyOverlap = 0.2)`.
#'
#' @return Numeric scalar in \[0, 1\].
#'
#' @noRd
.computeMAAS <- function(
    scores,
    weights = c(
        FOS = 0.4,
        jaccard = 0.2,
        brayCurtis = 0.2,
        redundancyOverlap = 0.2
    )
) {
    cli::cli_abort(
        "Not yet implemented (Phase 1).",
        .internal = TRUE
    )
}

#' Identify pathway correspondences
#'
#' Given two binary matrices in the same metabolite space,
#' identifies shared pathways (edges present in both), pathways
#' unique to the query, and pathways unique to the reference.
#' Returns results as tibbles with columns: consumed, produced,
#' source (species involved).
#'
#' @param xBin Sparse binary matrix (query, union space).
#' @param yBin Sparse binary matrix (reference, union space).
#' @param xEdges Edge data.frame from the query CM.
#' @param yEdges Edge data.frame from the reference CM.
#'
#' @return A named list with elements:
#'   \describe{
#'     \item{`shared`}{Tibble of pathways in both CMs.}
#'     \item{`uniqueQuery`}{Tibble of pathways only in query.}
#'     \item{`uniqueReference`}{Tibble of pathways only in
#'       reference.}
#'   }
#'
#' @noRd
.identifyPathwayCorrespondences <- function(xBin, yBin,
                                            xEdges, yEdges) {
    cli::cli_abort(
        "Not yet implemented (Phase 1).",
        .internal = TRUE
    )
}

#' Compute permutation p-value
#'
#' Estimates the probability of observing a score as extreme as
#' `observed` under a null model of degree-preserving network
#' rewiring. Uses `igraph::rewire(keeping_degseq())` to generate
#' permuted networks.
#'
#' @param xGraph An `igraph` graph object (query network).
#' @param yBin Sparse binary matrix (reference, union space).
#' @param observed Numeric scalar; the observed metric value.
#' @param metricFn A function that takes two sparse binary
#'   matrices and returns a numeric score (e.g.,
#'   `.functionalOverlap`).
#' @param nPerm Integer; number of permutations. Default `999L`.
#' @param metabolites Character vector of metabolite names in
#'   universal space.
#'
#' @return Numeric scalar: (sum(null >= observed) + 1) /
#'   (nPerm + 1).
#'
#' @noRd
.computePvalue <- function(xGraph, yBin, observed, metricFn,
                           nPerm = 999L, metabolites = NULL) {
    cli::cli_abort(
        "Not yet implemented (Phase 1).",
        .internal = TRUE
    )
}
