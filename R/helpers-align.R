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
#' pathway matrices. Also used by CMS overlap matrix construction.
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
    if (denom == 0) {
        return(0)
    }
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
    if (union_size == 0) {
        return(0)
    }
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
    if (total == 0) {
        return(0)
    }
    max(0, 1 - (numerator / total))
}

#' Redundancy Overlap
#'
#' Computes weighted Jaccard on species-count (nSpecies)
#' matrices. Measures how similarly two consortia distribute
#' metabolic labor across species for shared pathways.
#'
#' @param xSpecies nSpecies assay matrix for query (union
#'   space).
#' @param ySpecies nSpecies assay matrix for reference (union
#'   space).
#'
#' @return Numeric scalar in \[0, 1\].
#'
#' @noRd
.redundancyOverlap <- function(xSpecies, ySpecies) {
    diff_mat <- abs(xSpecies - ySpecies)
    total <- xSpecies + ySpecies
    numer <- sum(total - diff_mat) / 2
    denom <- sum(total + diff_mat) / 2
    if (denom == 0) {
        return(0)
    }
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
#'   `"nSpecies"`. `NULL` if unweighted.
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
.computeAllScores <- function(xBin, yBin, xWeighted = NULL, yWeighted = NULL) {
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
            xWeighted$nSpecies,
            yWeighted$nSpecies
        )
    } else {
        scores$brayCurtis <- NA_real_
        scores$redundancyOverlap <- NA_real_
    }

    scores
}

#' Compute pathway coverage ratios
#'
#' Returns the fraction of each network's pathways that are
#' shared. Complements FOS by revealing size asymmetry: FOS = 1
#' with low reference coverage means the query is a strict
#' subset of the reference.
#'
#' @param xBin Sparse binary matrix (query, union space).
#' @param yBin Sparse binary matrix (reference, union space).
#'
#' @return Named list with `coverageQuery` and
#'   `coverageReference`, each a numeric scalar in
#'   \[0, 1\].
#'
#' @noRd
.computeCoverage <- function(xBin, yBin) {
    nShared <- sum(xBin * yBin)
    nQuery <- sum(xBin)
    nRef <- sum(yBin)
    list(
        coverageQuery = if (nQuery == 0) {
            0
        } else {
            nShared / nQuery
        },
        coverageReference = if (nRef == 0) {
            0
        } else {
            nShared / nRef
        }
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
    if (abs(sum(weights) - 1) > 1e-10) {
        cli::cli_abort("Weights must sum to 1.")
    }

    available <- names(scores)[
        !is.na(vapply(scores, identity, numeric(1L)))
    ]
    if (length(available) == 0L) {
        return(NA_real_)
    }

    w <- weights[available]
    w <- w / sum(w)
    s <- vapply(scores[available], identity, numeric(1L))
    sum(w * s)
}

#' Identify pathway correspondences
#'
#' Given two binary matrices in the same metabolite space,
#' identifies shared pathways (present in both), pathways
#' unique to the query, and pathways unique to the reference.
#' Returns results as tibbles with columns: consumed, produced,
#' source (species involved).
#'
#' @param xBin Sparse binary matrix (query, union space).
#' @param yBin Sparse binary matrix (reference, union space).
#' @param xPathways Pathway data.frame from the query CM.
#' @param yPathways Pathway data.frame from the reference CM.
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
.identifyPathwayCorrespondences <- function(xBin, yBin, xPathways, yPathways) {
    metabolites <- rownames(xBin)

    shared_mat <- Matrix::drop0(xBin * yBin)
    unique_x <- Matrix::drop0(xBin - shared_mat)
    unique_y <- Matrix::drop0(yBin - shared_mat)

    .matToPathwayList <- function(mat) {
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

    shared_pw <- .matToPathwayList(shared_mat)
    unique_x_pw <- .matToPathwayList(unique_x)
    unique_y_pw <- .matToPathwayList(unique_y)

    ## Enrich shared pathways with species info
    if (nrow(shared_pw) > 0L) {
        shared_pw <- shared_pw |>
            dplyr::left_join(
                dplyr::select(
                    xPathways,
                    "consumed",
                    "produced",
                    "data"
                ),
                by = c("consumed", "produced")
            ) |>
            dplyr::rename(querySpecies = "data") |>
            dplyr::left_join(
                dplyr::select(
                    yPathways,
                    "consumed",
                    "produced",
                    "data"
                ),
                by = c("consumed", "produced")
            ) |>
            dplyr::rename(referenceSpecies = "data")
    }

    list(
        shared = shared_pw,
        uniqueQuery = unique_x_pw,
        uniqueReference = unique_y_pw
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
.computePvalue <- function(
    xGraph,
    yBin,
    observed,
    metricFn,
    nPerm = 999L,
    metabolites = NULL
) {
    if (nPerm == 0L) {
        return(NA_real_)
    }

    n_iter <- igraph::ecount(xGraph) * 10L

    null_scores <- vapply(
        seq_len(nPerm),
        function(i) {
            rewired <- igraph::rewire(
                xGraph,
                igraph::keeping_degseq(niter = n_iter)
            )
            rewired_adj <- igraph::as_adjacency_matrix(
                rewired,
                sparse = TRUE
            )
            rewired_exp <- .expandMatrix(rewired_adj, metabolites)
            metricFn(rewired_exp, yBin)
        },
        numeric(1L)
    )

    (sum(null_scores >= observed) + 1L) / (nPerm + 1L)
}

## ---- Multiple alignment helpers -----------------------------------------

#' Pre-expand weighted assays to universal metabolite space
#'
#' Expands Consumption, Production, and nSpecies assay matrices
#' for every consortium in a CMS to the universal metabolite
#' space. Called once at the start of non-FOS multiple alignment
#' to avoid redundant per-pair expansions.
#'
#' @param cms A [ConsortiumMetabolismSet] object.
#'
#' @return A named list (keyed by consortium name), each element
#'   a named list with `Consumption`, `Production`, and
#'   `nSpecies` sparse matrices in universal space. Elements are
#'   `NULL` for unweighted consortia.
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
            nSpecies = .expandMatrix(
                assays(cm)$nSpecies,
                universal_mets
            ),
            Consumption = .expandMatrix(
                assays(cm)$Consumption,
                universal_mets
            ),
            Production = .expandMatrix(
                assays(cm)$Production,
                universal_mets
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
.computePairwiseSimilarityMatrix <- function(cms, method, BPPARAM) {
    n <- length(cms@Consortia)
    cm_names <- names(cms@BinaryMatrices)
    bins <- cms@BinaryMatrices

    ## FOS shortcut: invert pre-computed distance matrix
    if (method == "FOS") {
        om <- cms@OverlapMatrix
        sim_mat <- 1 - om
        return(sim_mat)
    }

    ## Pre-expand weighted assays if needed
    weighted <- NULL
    needs_weighted <- method %in%
        c(
            "brayCurtis",
            "redundancyOverlap",
            "MAAS"
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
                    weighted[[i]],
                    weighted[[j]]
                )
            } else if (method == "redundancyOverlap") {
                .redundancyOverlap(
                    weighted[[i]]$nSpecies,
                    weighted[[j]]$nSpecies
                )
            } else {
                ## MAAS
                all_sc <- .computeAllScores(
                    bins[[i]],
                    bins[[j]],
                    weighted[[i]],
                    weighted[[j]]
                )
                .computeMAAS(all_sc)
            }
        },
        BPPARAM = BPPARAM
    )

    ## Fill symmetric matrix
    sim_mat <- matrix(
        0,
        n,
        n,
        dimnames = list(cm_names, cm_names)
    )
    diag(sim_mat) <- 1
    scores_vec <- vapply(
        scores,
        identity,
        numeric(1L)
    )
    for (k in seq_len(ncol(pairs))) {
        i <- pairs[1L, k]
        j <- pairs[2L, k]
        sim_mat[i, j] <- scores_vec[k]
        sim_mat[j, i] <- scores_vec[k]
    }
    sim_mat
}

#' Summarise pairwise similarity scores
#'
#' Reduces the upper triangle of a square similarity matrix to a
#' named list of summary statistics used as the multi-alignment
#' `Scores` slot.
#'
#' @param simMat Numeric n x n similarity matrix (symmetric, with
#'   `diag = 1`).
#'
#' @return Named list with `mean`, `median`, `min`, `max`, `sd`,
#'   and `nPairs` (integer count of upper-triangle entries).
#'
#' @keywords internal
#' @noRd
.summariseSimilarityScores <- function(simMat) {
    pairwiseVals <- simMat[upper.tri(simMat)]
    list(
        mean = mean(pairwiseVals),
        median = stats::median(pairwiseVals),
        min = min(pairwiseVals),
        max = max(pairwiseVals),
        sd = stats::sd(pairwiseVals),
        nPairs = length(pairwiseVals)
    )
}

#' Build a dendrogram from a similarity matrix
#'
#' Converts a similarity matrix to a distance object via
#' `1 - simMat`, then runs `stats::hclust()` with the requested
#' linkage and coerces to a dendrogram.
#'
#' @param simMat Numeric n x n similarity matrix.
#' @param linkage Character scalar; agglomeration method passed to
#'   `stats::hclust()`.
#'
#' @return A `dendrogram` object.
#'
#' @keywords internal
#' @noRd
.dendrogramFromSimilarity <- function(simMat, linkage) {
    distMat <- stats::as.dist(1 - simMat)
    stats::hclust(distMat, method = linkage) |>
        stats::as.dendrogram()
}

#' Build pathway prevalence from a CMS
#'
#' Extracts the Levels assay from a CMS (which counts how many
#' consortia have each metabolite-metabolite pathway) and
#' returns a tidy data.frame of pathway prevalence.
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

## ---- Database search helpers --------------------------------------------

#' Expand the Consumption/Production/nSpecies triple for a CM
#'
#' Returns the three weighted assay matrices for `cm`, expanded
#' into the supplied `unionMets` space, or `NULL` for unweighted
#' consortia.
#'
#' @param cm A [ConsortiumMetabolism].
#' @param unionMets Character vector of target metabolite names.
#'
#' @return Named list with `Consumption`, `Production`, and
#'   `nSpecies` sparse matrices, or `NULL` if `cm@Weighted` is
#'   `FALSE`.
#'
#' @keywords internal
#' @noRd
.expandWeightedTriple <- function(cm, unionMets) {
    if (!cm@Weighted) {
        return(NULL)
    }
    list(
        Consumption = .expandMatrix(
            assays(cm)$Consumption,
            unionMets
        ),
        Production = .expandMatrix(
            assays(cm)$Production,
            unionMets
        ),
        nSpecies = .expandMatrix(
            assays(cm)$nSpecies,
            unionMets
        )
    )
}

#' Expand query and database matrices into a shared metabolite space
#'
#' Builds the union metabolite space across the query CM and a
#' database CMS, then expands the query binary matrix and every
#' database consortium's binary matrix into that shared space.
#' Optionally also expands the weighted assays (Consumption,
#' Production, nSpecies) when the requested metric set requires
#' them.
#'
#' @param x A [ConsortiumMetabolism] (query).
#' @param y A [ConsortiumMetabolismSet] (database).
#' @param metrics Character vector of metric names to be computed.
#' @param method Character scalar; primary metric (used to decide
#'   whether weighted assays are required even if not explicitly
#'   in `metrics`, e.g. `"MAAS"`).
#'
#' @return Named list with elements `xBin`, `xWeighted` (or
#'   `NULL`), `yBins` (named list of binary matrices keyed by
#'   consortium name), `yWeightedList` (named list, may contain
#'   `NULL` per element, or `NULL` overall), `unionMets`
#'   (character vector), and `cmNames` (character vector).
#'
#' @keywords internal
#' @noRd
.alignSearchExpandSpace <- function(x, y, metrics, method) {
    cmsUniversal <- rownames(y@BinaryMatrices[[1L]])
    queryMets <- rownames(assays(x)$Binary)
    unionMets <- sort(unique(c(cmsUniversal, queryMets)))

    xBin <- .expandMatrix(assays(x)$Binary, unionMets)

    needsWeighted <- any(
        c("brayCurtis", "redundancyOverlap") %in% metrics
    ) ||
        method == "MAAS"

    xWeighted <- if (needsWeighted) {
        .expandWeightedTriple(x, unionMets)
    } else {
        NULL
    }

    cmNames <- names(y@BinaryMatrices)
    yBins <- lapply(
        y@BinaryMatrices,
        function(m) .expandMatrix(m, unionMets)
    )

    yWeightedList <- NULL
    if (needsWeighted) {
        yWeightedList <- lapply(
            y@Consortia,
            function(cm) .expandWeightedTriple(cm, unionMets)
        )
        names(yWeightedList) <- cmNames
    }

    list(
        xBin = xBin,
        xWeighted = xWeighted,
        yBins = yBins,
        yWeightedList = yWeightedList,
        unionMets = unionMets,
        cmNames = cmNames
    )
}

#' Score a query against every consortium in a database
#'
#' Runs `.computeAllScores()` and `.computeCoverage()` for every
#' database consortium, optionally in parallel via BiocParallel.
#'
#' @param xBin Sparse binary matrix for the query in union space.
#' @param xWeighted Named list of expanded weighted assays for
#'   the query (or `NULL`).
#' @param yBins Named list of expanded binary matrices for each
#'   database consortium.
#' @param yWeightedList Named list of expanded weighted assay
#'   lists for each database consortium (or `NULL`).
#' @param BPPARAM A [BiocParallel::BiocParallelParam] object.
#'
#' @return A list (length = `length(yBins)`) of named scores
#'   lists. Each element contains `FOS`, `jaccard`, `brayCurtis`,
#'   `redundancyOverlap`, `coverageQuery`, and
#'   `coverageReference`.
#'
#' @keywords internal
#' @noRd
.alignSearchScoreAll <- function(
    xBin,
    xWeighted,
    yBins,
    yWeightedList,
    BPPARAM
) {
    BiocParallel::bplapply(
        seq_along(yBins),
        function(i) {
            yB <- yBins[[i]]
            yW <- if (is.null(yWeightedList)) {
                NULL
            } else {
                yWeightedList[[i]]
            }
            sc <- .computeAllScores(xBin, yB, xWeighted, yW)
            cov <- .computeCoverage(xBin, yB)
            sc$coverageQuery <- cov$coverageQuery
            sc$coverageReference <- cov$coverageReference
            sc
        },
        BPPARAM = BPPARAM
    )
}

#' Reduce per-consortium scores to a primary-score vector
#'
#' For each entry in `scoresList`, picks the named primary metric
#' (or, for `"MAAS"`, recomputes the composite from the four core
#' metrics).
#'
#' @param scoresList List of per-consortium score lists.
#' @param method Character scalar; primary metric.
#'
#' @return Numeric vector of primary scores, aligned with
#'   `scoresList`.
#'
#' @keywords internal
#' @noRd
.alignSearchPrimaryScores <- function(scoresList, method) {
    vapply(
        scoresList,
        function(s) {
            if (method == "MAAS") {
                core <- s[c(
                    "FOS",
                    "jaccard",
                    "brayCurtis",
                    "redundancyOverlap"
                )]
                .computeMAAS(core)
            } else {
                s[[method]]
            }
        },
        numeric(1L)
    )
}

#' Build the ranking tibble and similarity matrix for a search
#'
#' Assembles the per-consortium scores (with primary scores and
#' coverage ratios) into a ranking tibble sorted by primary score,
#' optionally truncated to `topK`. Also builds the 1-row
#' similarity matrix labelled with the query name and identifies
#' the top hit (always the overall best, regardless of `topK`).
#'
#' @param scoresList List of per-consortium score lists, as
#'   returned by `.alignSearchScoreAll()`.
#' @param primaryScores Numeric vector of primary scores aligned
#'   with `scoresList`, as returned by
#'   `.alignSearchPrimaryScores()`.
#' @param cmNames Character vector of consortium names (database).
#' @param topK Integer or `NULL`; if non-`NULL`, truncate the
#'   ranked results to the top `topK` hits.
#' @param queryLabel Character scalar; row name for the
#'   similarity matrix.
#'
#' @return Named list with `ranking` (tibble), `topIdx` (integer
#'   index into `scoresList` of the overall best hit, regardless
#'   of `topK`), `topName` (character), `topScore` (numeric), and
#'   `simMat` (1 x nrow(ranking) matrix).
#'
#' @keywords internal
#' @noRd
.alignSearchBuildRanking <- function(
    scoresList,
    primaryScores,
    cmNames,
    topK,
    queryLabel
) {
    pluck <- function(field) {
        vapply(scoresList, function(s) s[[field]], numeric(1L))
    }
    ranking <- tibble::tibble(
        reference = cmNames,
        score = primaryScores,
        FOS = pluck("FOS"),
        jaccard = pluck("jaccard"),
        brayCurtis = pluck("brayCurtis"),
        redundancyOverlap = pluck("redundancyOverlap"),
        coverageQuery = pluck("coverageQuery"),
        coverageReference = pluck("coverageReference")
    ) |>
        dplyr::arrange(dplyr::desc(.data$score))

    if (!is.null(topK)) {
        ranking <- utils::head(ranking, topK)
    }

    topIdx <- which.max(primaryScores)
    simMat <- matrix(
        ranking$score,
        nrow = 1L,
        ncol = nrow(ranking),
        dimnames = list(queryLabel, ranking$reference)
    )

    list(
        ranking = ranking,
        topIdx = topIdx,
        topName = cmNames[topIdx],
        topScore = primaryScores[topIdx],
        simMat = simMat
    )
}

#' Permutation p-value for the top hit of a database search
#'
#' Wraps `.computePvalue()` for the top-scoring database
#' consortium. Emits a warning and returns `NA_real_` for metrics
#' where degree-preserving permutation is not yet supported
#' (`"brayCurtis"`, `"redundancyOverlap"`).
#'
#' @param method Character scalar; primary metric.
#' @param x The query [ConsortiumMetabolism] object (needs
#'   `@Graphs[[1L]]` for permutation).
#' @param xWeighted Expanded weighted assays for the query (or
#'   `NULL`); used to recompute MAAS components on permuted nets.
#' @param yBin Sparse binary matrix for the top-hit database
#'   consortium in union space.
#' @param topYWeighted Expanded weighted assays for the top-hit
#'   database consortium (or `NULL`).
#' @param topScore Numeric; observed primary score for the top
#'   hit.
#' @param nPermutations Integer.
#' @param unionMets Character vector of union metabolite names.
#'
#' @return Numeric p-value in \[0, 1\], or `NA_real_` for
#'   unsupported metrics.
#'
#' @keywords internal
#' @noRd
.alignSearchTopPvalue <- function(
    method,
    x,
    xWeighted,
    yBin,
    topYWeighted,
    topScore,
    nPermutations,
    unionMets
) {
    if (method %in% c("brayCurtis", "redundancyOverlap")) {
        cli::cli_warn(
            paste0(
                "P-value computation for ",
                "{.val {method}} not yet ",
                "supported. Skipping."
            )
        )
        return(NA_real_)
    }
    metricFn <- switch(
        method,
        FOS = .functionalOverlap,
        jaccard = .jaccardIndex,
        MAAS = function(xBinPerm, yBinTop) {
            sc <- .computeAllScores(
                xBinPerm,
                yBinTop,
                xWeighted,
                topYWeighted
            )
            .computeMAAS(sc)
        }
    )
    .computePvalue(
        xGraph = x@Graphs[[1L]],
        yBin = yBin,
        observed = topScore,
        metricFn = metricFn,
        nPerm = nPermutations,
        metabolites = unionMets
    )
}
