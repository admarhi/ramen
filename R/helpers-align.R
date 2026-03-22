## ----
## align-helpers.R
## Internal helper functions for the alignment system.
## None of these are exported.
## ----

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
#' @return A sparse matrix of dimensions
#'   `length(target_mets) x length(target_mets)` with the
#'   original values placed at the correct positions and zeros
#'   elsewhere.
#'
#' @noRd
.expandMatrix <- function(mat, target_mets) {
    local_mets <- rownames(mat)
    n <- length(target_mets)

    ## Coerce dense to sparse if needed
    if (!is(mat, "dgCMatrix")) {
        mat <- as(mat, "dgCMatrix")
    }

    ## Extract triplet form
    triplet <- Matrix::summary(mat)
    if (nrow(triplet) == 0L) {
        return(
            Matrix::sparseMatrix(
                i = integer(0L),
                j = integer(0L),
                x = numeric(0L),
                dims = c(n, n),
                dimnames = list(target_mets, target_mets)
            )
        )
    }

    ## Map local indices to target indices
    idx <- match(local_mets, target_mets)

    Matrix::sparseMatrix(
        i = idx[triplet$i],
        j = idx[triplet$j],
        x = triplet$x,
        dims = c(n, n),
        dimnames = list(target_mets, target_mets)
    )
}

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
    mets1 <- rownames(assays(cm1)$Binary)
    mets2 <- rownames(assays(cm2)$Binary)
    union_mets <- sort(unique(c(mets1, mets2)))

    X <- .expandMatrix(assays(cm1)$Binary, union_mets)
    Y <- .expandMatrix(assays(cm2)$Binary, union_mets)

    list(X = X, Y = Y, metabolites = union_mets)
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
    intersection <- xBin * yBin
    denom <- min(sum(xBin), sum(yBin))
    if (denom == 0) return(0)
    sum(intersection) / denom
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
    intersection <- sum(xBin * yBin)
    union_size <- sum(xBin) + sum(yBin) - intersection
    if (union_size == 0) return(0)
    intersection / union_size
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
    xC <- xWeighted$Consumption
    yC <- yWeighted$Consumption
    xP <- xWeighted$Production
    yP <- yWeighted$Production

    diff_C <- abs(xC - yC)
    diff_P <- abs(xP - yP)
    numerator <- sum(diff_C) + sum(diff_P)

    total <- sum(xC) + sum(yC) + sum(xP) + sum(yP)
    if (total == 0) return(0)
    1 - (numerator / total)
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
    diff_mat <- abs(xEdges - yEdges)
    total <- xEdges + yEdges
    numer <- sum(total - diff_mat) / 2
    denom <- sum(total + diff_mat) / 2
    if (denom == 0) return(0)
    numer / denom
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
#'   query. Expected names: `"Consumption"`, `"Production"`,
#'   `"nEdges"`. `NULL` if unweighted.
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
    scores <- list(
        FOS = .functionalOverlap(xBin, yBin),
        jaccard = .jaccardIndex(xBin, yBin)
    )

    if (!is.null(xWeighted) && !is.null(yWeighted)) {
        scores$brayCurtis <- .brayCurtisSimilarity(
            list(
                Consumption = xWeighted$Consumption,
                Production = xWeighted$Production
            ),
            list(
                Consumption = yWeighted$Consumption,
                Production = yWeighted$Production
            )
        )
        scores$redundancyOverlap <- .redundancyOverlap(
            xWeighted$nEdges, yWeighted$nEdges
        )
    } else {
        scores$brayCurtis <- NA_real_
        scores$redundancyOverlap <- NA_real_
    }

    scores
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
    if (abs(sum(weights) - 1) > 1e-10) {
        cli::cli_abort("Weights must sum to 1.")
    }

    available <- names(scores)[
        !is.na(vapply(scores, identity, numeric(1L)))
    ]
    if (length(available) == 0L) return(NA_real_)

    w <- weights[available]
    w <- w / sum(w)
    s <- vapply(scores[available], identity, numeric(1L))
    sum(w * s)
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
    metabolites <- rownames(xBin)

    shared_mat <- Matrix::drop0(xBin * yBin)
    unique_x <- Matrix::drop0(xBin - shared_mat)
    unique_y <- Matrix::drop0(yBin - shared_mat)

    .matToEdgeList <- function(mat) {
        triplet <- Matrix::summary(mat)
        if (nrow(triplet) == 0L) {
            return(tibble::tibble(
                consumed = character(0L),
                produced = character(0L)
            ))
        }
        tibble::tibble(
            consumed = metabolites[triplet$i],
            produced = metabolites[triplet$j]
        )
    }

    shared_edges <- .matToEdgeList(shared_mat)
    unique_x_edges <- .matToEdgeList(unique_x)
    unique_y_edges <- .matToEdgeList(unique_y)

    ## Enrich shared edges with species info
    if (nrow(shared_edges) > 0L) {
        shared_edges <- shared_edges |>
            dplyr::left_join(
                dplyr::select(
                    xEdges, "consumed", "produced", "data"
                ),
                by = c("consumed", "produced")
            ) |>
            dplyr::rename(querySpecies = "data") |>
            dplyr::left_join(
                dplyr::select(
                    yEdges, "consumed", "produced", "data"
                ),
                by = c("consumed", "produced")
            ) |>
            dplyr::rename(referenceSpecies = "data")
    }

    list(
        shared = shared_edges,
        uniqueQuery = unique_x_edges,
        uniqueReference = unique_y_edges
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
    if (nPerm == 0L) return(NA_real_)

    n_iter <- igraph::ecount(xGraph) * 10L

    null_scores <- vapply(seq_len(nPerm), function(i) {
        rewired <- igraph::rewire(
            xGraph,
            igraph::keeping_degseq(niter = n_iter)
        )
        rewired_adj <- igraph::as_adjacency_matrix(
            rewired, sparse = TRUE
        )
        rewired_exp <- .expandMatrix(rewired_adj, metabolites)
        metricFn(rewired_exp, yBin)
    }, numeric(1L))

    (sum(null_scores >= observed) + 1L) / (nPerm + 1L)
}

## ---- Multiple alignment helpers -----------------------------------------

#' Pre-expand weighted assays to universal metabolite space
#'
#' Expands Consumption, Production, and nEdges assay matrices
#' for every consortium in a CMS to the universal metabolite
#' space. Called once at the start of non-FOS multiple alignment
#' to avoid redundant per-pair expansions.
#'
#' @param cms A [ConsortiumMetabolismSet] object.
#'
#' @return A named list (keyed by consortium name), each element
#'   a named list with `Consumption`, `Production`, and `nEdges`
#'   sparse matrices in universal space. Elements are `NULL` for
#'   unweighted consortia.
#'
#' @noRd
.expandWeightedAssays <- function(cms) {
    universal_mets <- rownames(cms@BinaryMatrices[[1L]])
    cm_names <- names(cms@BinaryMatrices)

    result <- lapply(cms@Consortia, function(cm) {
        if (!cm@Weighted) {
            return(NULL)
        }
        list(
            nEdges = .expandMatrix(
                assays(cm)$nEdges, universal_mets
            ),
            Consumption = .expandMatrix(
                assays(cm)$Consumption, universal_mets
            ),
            Production = .expandMatrix(
                assays(cm)$Production, universal_mets
            )
        )
    })
    names(result) <- cm_names
    result
}

#' Compute pairwise similarity matrix for a CMS
#'
#' Builds an n x n symmetric similarity matrix across all
#' consortia in a CMS using the requested metric. For FOS,
#' reuses the pre-computed OverlapMatrix. For other metrics,
#' computes all n*(n-1)/2 pairwise scores, optionally in
#' parallel via BiocParallel.
#'
#' @param cms A [ConsortiumMetabolismSet] object.
#' @param method Character scalar: metric name.
#' @param BPPARAM A [BiocParallel::BiocParallelParam] object.
#'
#' @return Numeric n x n matrix with dimnames = consortium
#'   names, diagonal = 1, off-diagonal = similarity scores.
#'
#' @noRd
.computePairwiseSimilarityMatrix <- function(cms, method,
                                             BPPARAM) {
    n <- length(cms@Consortia)
    cm_names <- names(cms@BinaryMatrices)
    bins <- cms@BinaryMatrices

    ## FOS shortcut: invert pre-computed distance matrix
    if (method == "FOS") {
        om <- cms@OverlapMatrix
        sim_mat <- 1 - (om + t(om))
        return(sim_mat)
    }

    ## Pre-expand weighted assays if needed
    weighted <- NULL
    needs_weighted <- method %in% c(
        "brayCurtis", "redundancyOverlap", "MAAS"
    )
    if (needs_weighted) {
        weighted <- .expandWeightedAssays(cms)
    }

    ## Generate pair indices
    pairs <- utils::combn(n, 2L)

    ## Compute pairwise scores
    scores <- BiocParallel::bplapply(
        seq_len(ncol(pairs)),
        function(k) {
            i <- pairs[1L, k]
            j <- pairs[2L, k]
            if (method == "jaccard") {
                .jaccardIndex(bins[[i]], bins[[j]])
            } else if (method == "brayCurtis") {
                .brayCurtisSimilarity(
                    weighted[[i]], weighted[[j]]
                )
            } else if (method == "redundancyOverlap") {
                .redundancyOverlap(
                    weighted[[i]]$nEdges,
                    weighted[[j]]$nEdges
                )
            } else {
                ## MAAS
                all_sc <- .computeAllScores(
                    bins[[i]], bins[[j]],
                    weighted[[i]], weighted[[j]]
                )
                .computeMAAS(all_sc)
            }
        },
        BPPARAM = BPPARAM
    )

    ## Fill symmetric matrix
    sim_mat <- matrix(
        0, n, n,
        dimnames = list(cm_names, cm_names)
    )
    diag(sim_mat) <- 1
    scores_vec <- vapply(
        scores, identity, numeric(1L)
    )
    for (k in seq_len(ncol(pairs))) {
        i <- pairs[1L, k]
        j <- pairs[2L, k]
        sim_mat[i, j] <- scores_vec[k]
        sim_mat[j, i] <- scores_vec[k]
    }
    sim_mat
}

#' Build edge prevalence from a CMS
#'
#' Extracts the Levels assay from a CMS (which counts how many
#' consortia have each metabolite-metabolite edge) and returns
#' a tidy data.frame of edge prevalence.
#'
#' @param cms A [ConsortiumMetabolismSet] object.
#'
#' @return A data.frame with columns:
#'   \describe{
#'     \item{`consumed`}{Character, consumed metabolite name.}
#'     \item{`produced`}{Character, produced metabolite name.}
#'     \item{`nConsortia`}{Integer, number of consortia with
#'       this edge.}
#'     \item{`proportion`}{Numeric, nConsortia / total
#'       consortia.}
#'   }
#'
#' @noRd
.buildPrevalence <- function(cms) {
    levels_mat <- assays(cms)$Levels
    mets <- rownames(levels_mat)
    n <- length(cms@Consortia)

    idx <- which(levels_mat > 0, arr.ind = TRUE)
    if (nrow(idx) == 0L) {
        return(data.frame(
            consumed = character(0L),
            produced = character(0L),
            nConsortia = integer(0L),
            proportion = numeric(0L),
            stringsAsFactors = FALSE
        ))
    }

    data.frame(
        consumed = mets[idx[, 1L]],
        produced = mets[idx[, 2L]],
        nConsortia = as.integer(levels_mat[idx]),
        proportion = levels_mat[idx] / n,
        stringsAsFactors = FALSE
    )
}
