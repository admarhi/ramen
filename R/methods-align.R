#' @include AllClasses.R AllGenerics.R helpers-align.R
NULL

#' @describeIn align Pairwise alignment of two
#'   [ConsortiumMetabolism] objects
#'
#' @param x A [ConsortiumMetabolism] object (query).
#' @param y A [ConsortiumMetabolism] object (reference).
#' @param method Character scalar specifying the metric.
#' @param computePvalue Logical; compute permutation p-value?
#'   Default `FALSE`.
#' @param nPermutations Integer; number of permutations.
#'   Default `999L`.
#' @param ... Additional arguments (currently unused).
#'
#' @return A [ConsortiumMetabolismAlignment] object of type
#'   `"pairwise"`.
#'
#' @export
setMethod(
    "align",
    signature(
        x = "ConsortiumMetabolism",
        y = "ConsortiumMetabolism"
    ),
    function(
        x,
        y,
        method = "FOS",
        computePvalue = FALSE,
        nPermutations = 999L,
        ...
    ) {
        ## 1. Validate method argument
        method <- match.arg(
            method,
            c("FOS", "jaccard", "brayCurtis", "redundancyOverlap", "MAAS")
        )

        ## 2. Harmonize metabolite space (binary matrices)
        harmonized <- .harmonizeMetaboliteSpace(x, y)
        xBin <- harmonized$X
        yBin <- harmonized$Y
        union_mets <- harmonized$metabolites

        ## 3. Expand weighted assays if both CMs are weighted
        xWeighted <- NULL
        yWeighted <- NULL
        if (x@Weighted && y@Weighted) {
            xWeighted <- list(
                Consumption = .expandMatrix(
                    assays(x)$Consumption,
                    union_mets
                ),
                Production = .expandMatrix(
                    assays(x)$Production,
                    union_mets
                ),
                nSpecies = .expandMatrix(
                    assays(x)$nSpecies,
                    union_mets
                )
            )
            yWeighted <- list(
                Consumption = .expandMatrix(
                    assays(y)$Consumption,
                    union_mets
                ),
                Production = .expandMatrix(
                    assays(y)$Production,
                    union_mets
                ),
                nSpecies = .expandMatrix(
                    assays(y)$nSpecies,
                    union_mets
                )
            )
        }

        ## 4. Compute all scores
        all_scores <- .computeAllScores(
            xBin,
            yBin,
            xWeighted,
            yWeighted
        )

        ## 5. Determine primary score
        if (method == "MAAS") {
            dots <- list(...)
            weights <- dots$weights
            if (is.null(weights)) {
                primary <- .computeMAAS(all_scores)
            } else {
                primary <- .computeMAAS(
                    all_scores,
                    weights = weights
                )
            }
        } else {
            primary <- all_scores[[method]]
        }

        ## 6. Identify pathway correspondences
        correspondences <- .identifyPathwayCorrespondences(
            xBin,
            yBin,
            x@Pathways,
            y@Pathways
        )

        ## 7. Optional: compute p-value
        pval <- NA_real_
        if (computePvalue) {
            if (method %in% c("brayCurtis", "redundancyOverlap")) {
                cli::cli_warn(
                    paste0(
                        "P-value computation for ",
                        "{.val {method}} not yet ",
                        "supported. Skipping."
                    )
                )
            } else {
                metric_fn <- switch(
                    method,
                    FOS = .functionalOverlap,
                    jaccard = .jaccardIndex,
                    MAAS = {
                        w <- list(...)$weights
                        function(xBinPerm, yBin) {
                            scores <- .computeAllScores(
                                xBinPerm,
                                yBin,
                                xWeighted,
                                yWeighted
                            )
                            if (is.null(w)) {
                                .computeMAAS(scores)
                            } else {
                                .computeMAAS(
                                    scores,
                                    weights = w
                                )
                            }
                        }
                    }
                )
                pval <- .computePvalue(
                    xGraph = x@Graphs[[1L]],
                    yBin = yBin,
                    observed = primary,
                    metricFn = metric_fn,
                    nPerm = nPermutations,
                    metabolites = union_mets
                )
            }
        }

        ## 8. Build and return CMA
        ConsortiumMetabolismAlignment(
            Type = "pairwise",
            Metric = method,
            QueryName = x@Name,
            ReferenceName = y@Name,
            Scores = all_scores,
            PrimaryScore = primary,
            Pvalue = pval,
            SharedPathways = correspondences$shared,
            UniqueQuery = correspondences$uniqueQuery,
            UniqueReference = correspondences$uniqueReference,
            Params = list(
                computePvalue = computePvalue,
                nPermutations = nPermutations
            ),
            Pathways = correspondences$shared,
            Metabolites = data.frame(
                met = union_mets,
                stringsAsFactors = FALSE
            )
        )
    }
)

#' @describeIn align Multiple alignment across all consortia in a
#'   [ConsortiumMetabolismSet]
#'
#' @param x A [ConsortiumMetabolismSet] object.
#' @param y Must be `missing`.
#' @param method Character scalar specifying the metric.
#' @param BPPARAM A [BiocParallel::BiocParallelParam] object.
#'   Default `BiocParallel::SerialParam()`.
#' @param ... Additional arguments (currently unused).
#'
#' @return A [ConsortiumMetabolismAlignment] object of type
#'   `"multiple"`.
#'
#' @importFrom BiocParallel SerialParam bplapply
#' @export
setMethod(
    "align",
    signature(
        x = "ConsortiumMetabolismSet",
        y = "missing"
    ),
    function(x, y, method = "FOS", BPPARAM = BiocParallel::SerialParam(), ...) {
        ## 1. Validate
        method <- match.arg(
            method,
            c("FOS", "jaccard", "brayCurtis", "redundancyOverlap", "MAAS")
        )
        n <- length(x@Consortia)
        if (n < 2L) {
            cli::cli_abort(
                "Multiple alignment requires at least 2 \\
                 consortia."
            )
        }

        cm_names <- names(x@BinaryMatrices)
        cli::cli_inform(
            "Computing multiple alignment for {n} \\
             consortia using {.val {method}}."
        )

        ## 2. Pairwise similarity matrix
        sim_mat <- .computePairwiseSimilarityMatrix(
            x,
            method,
            BPPARAM
        )

        ## 3. Summary scores
        pairwise_vals <- sim_mat[upper.tri(sim_mat)]
        scores <- list(
            mean = mean(pairwise_vals),
            median = stats::median(pairwise_vals),
            min = min(pairwise_vals),
            max = max(pairwise_vals),
            sd = stats::sd(pairwise_vals),
            nPairs = length(pairwise_vals)
        )
        primary <- scores$median

        ## 4. Consensus network and prevalence
        prevalence <- .buildPrevalence(x)

        ## 5. Dendrogram
        dist_mat <- stats::as.dist(1 - sim_mat)
        dend <- stats::hclust(dist_mat) |>
            stats::as.dendrogram()

        ## 6. Build and return CMA
        ConsortiumMetabolismAlignment(
            Type = "multiple",
            Metric = method,
            Scores = scores,
            PrimaryScore = primary,
            SimilarityMatrix = sim_mat,
            ConsensusPathways = prevalence,
            Prevalence = prevalence,
            Dendrogram = list(dend),
            Pathways = prevalence,
            Metabolites = data.frame(
                met = rownames(
                    x@BinaryMatrices[[1L]]
                ),
                stringsAsFactors = FALSE
            ),
            Params = list(
                nConsortia = n,
                consortiaNames = cm_names,
                BPPARAM = as.character(
                    class(BPPARAM)
                )
            )
        )
    }
)

#' @describeIn align Database search -- align a
#'   [ConsortiumMetabolism] against all consortia in a
#'   [ConsortiumMetabolismSet]
#'
#' @param x A [ConsortiumMetabolism] object (query).
#' @param y A [ConsortiumMetabolismSet] object (database).
#' @param method Character scalar specifying the metric.
#' @param BPPARAM A [BiocParallel::BiocParallelParam] object.
#'   Default `BiocParallel::SerialParam()`.
#' @param ... Additional arguments (currently unused).
#'
#' @return A [ConsortiumMetabolismAlignment] object of type
#'   `"search"`.
#'
#' @export
setMethod(
    "align",
    signature(
        x = "ConsortiumMetabolism",
        y = "ConsortiumMetabolismSet"
    ),
    function(x, y, method = "FOS", BPPARAM = BiocParallel::SerialParam(), ...) {
        cli::cli_abort(
            "Database search not yet implemented (future)."
        )
    }
)
