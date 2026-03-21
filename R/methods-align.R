#' @include helpers-align.R

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
    function(x, y, method = "FOS",
             computePvalue = FALSE,
             nPermutations = 999L, ...) {

        ## 1. Validate method argument
        method <- match.arg(
            method,
            c("FOS", "jaccard", "brayCurtis",
              "redundancyOverlap", "MAAS")
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
                    assays(x)$Consumption, union_mets
                ),
                Production = .expandMatrix(
                    assays(x)$Production, union_mets
                ),
                nEdges = .expandMatrix(
                    assays(x)$nEdges, union_mets
                )
            )
            yWeighted <- list(
                Consumption = .expandMatrix(
                    assays(y)$Consumption, union_mets
                ),
                Production = .expandMatrix(
                    assays(y)$Production, union_mets
                ),
                nEdges = .expandMatrix(
                    assays(y)$nEdges, union_mets
                )
            )
        }

        ## 4. Compute all scores
        all_scores <- .computeAllScores(
            xBin, yBin, xWeighted, yWeighted
        )

        ## 5. Determine primary score
        if (method == "MAAS") {
            dots <- list(...)
            weights <- dots$weights
            if (is.null(weights)) {
                primary <- .computeMAAS(all_scores)
            } else {
                primary <- .computeMAAS(
                    all_scores, weights = weights
                )
            }
        } else {
            primary <- all_scores[[method]]
        }

        ## 6. Identify pathway correspondences
        correspondences <- .identifyPathwayCorrespondences(
            xBin, yBin, x@Edges, y@Edges
        )

        ## 7. Optional: compute p-value
        pval <- NA_real_
        if (computePvalue) {
            if (method %in% c("brayCurtis",
                              "redundancyOverlap")) {
                cli::cli_warn(
                    paste0(
                        "P-value computation for ",
                        "{.val {method}} not yet ",
                        "supported. Skipping."
                    )
                )
            } else {
                metric_fn <- switch(method,
                    FOS = .functionalOverlap,
                    jaccard = .jaccardIndex,
                    MAAS = .functionalOverlap
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
            Edges = correspondences$shared,
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
#' @importFrom BiocParallel SerialParam
#' @export
setMethod(
    "align",
    signature(
        x = "ConsortiumMetabolismSet",
        y = "missing"
    ),
    function(x, y, method = "FOS",
             BPPARAM = BiocParallel::SerialParam(),
             ...) {
        cli::cli_abort(
            "Multiple alignment not yet implemented (Phase 2)."
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
    function(x, y, method = "FOS",
             BPPARAM = BiocParallel::SerialParam(),
             ...) {
        cli::cli_abort(
            "Database search not yet implemented (future)."
        )
    }
)
