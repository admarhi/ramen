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

        ## 5b. Coverage ratios (subset detection)
        coverage <- .computeCoverage(xBin, yBin)
        all_scores$coverageQuery <-
            coverage$coverageQuery
        all_scores$coverageReference <-
            coverage$coverageReference

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
        auto_name <- paste0(x@Name, " vs ", y@Name)
        ConsortiumMetabolismAlignment(
            Name = auto_name,
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
#' @param linkage Character scalar specifying the agglomeration
#'   method for hierarchical clustering. Passed to
#'   \code{\link[stats]{hclust}}. One of `"complete"` (default),
#'   `"average"`, `"single"`, or `"ward.D2"`.
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
    function(
        x,
        y,
        method = "FOS",
        linkage = "complete",
        BPPARAM = BiocParallel::SerialParam(),
        ...
    ) {
        ## 1. Validate
        method <- match.arg(
            method,
            c("FOS", "jaccard", "brayCurtis", "redundancyOverlap", "MAAS")
        )
        linkage <- match.arg(
            linkage,
            c("complete", "average", "single", "ward.D2")
        )
        n <- length(x@Consortia)
        if (n < 2L) {
            cli::cli_abort(
                "Multiple alignment requires at \\
                least 2 consortia."
            )
        }

        cm_names <- names(x@BinaryMatrices)
        cli::cli_inform(
            "Computing multiple alignment for \\
            {n} consortia using {.val {method}}."
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
        dend <- stats::hclust(dist_mat, method = linkage) |>
            stats::as.dendrogram()

        ## 6. Build and return CMA
        auto_name <- paste0("Multiple (", n, " consortia)")
        ConsortiumMetabolismAlignment(
            Name = auto_name,
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
#' @param method Character scalar specifying the primary metric
#'   used to rank hits. One of `"FOS"` (default), `"jaccard"`,
#'   `"brayCurtis"`, `"redundancyOverlap"`, or `"MAAS"`.
#' @param metrics Character vector of metrics to compute per
#'   hit. Defaults to all four base metrics
#'   (`"FOS"`, `"jaccard"`, `"brayCurtis"`,
#'   `"redundancyOverlap"`). Restrict to a subset (e.g.
#'   `metrics = "FOS"`) to skip weighted-assay expansion and
#'   speed up large-database searches. `"MAAS"` as the primary
#'   metric requires all four components.
#' @param topK Integer; if non-`NULL`, truncate the ranked
#'   results (and `SimilarityMatrix`) to the top `topK` hits.
#'   The overall top hit's pathway correspondences are always
#'   reported regardless of `topK`. Default `NULL` (all hits).
#' @param computePvalue Logical; compute a permutation p-value
#'   for the top hit? Default `FALSE`. Not supported for
#'   `"brayCurtis"` or `"redundancyOverlap"`.
#' @param nPermutations Integer; number of permutations used
#'   when `computePvalue = TRUE`. Default `999L`.
#' @param BPPARAM A [BiocParallel::BiocParallelParam] object.
#'   Default `BiocParallel::SerialParam()`.
#' @param ... Additional arguments (currently unused).
#'
#' @return A [ConsortiumMetabolismAlignment] object of type
#'   `"search"`. `PrimaryScore`/`ReferenceName` point to the
#'   top hit; the full ranked hit table is stored in
#'   `Scores$ranking` and as a 1-row `SimilarityMatrix`.
#'   Pathway correspondences (`SharedPathways`, `UniqueQuery`,
#'   `UniqueReference`) are for the top hit only.
#'
#' @export
setMethod(
    "align",
    signature(
        x = "ConsortiumMetabolism",
        y = "ConsortiumMetabolismSet"
    ),
    function(
        x,
        y,
        method = "FOS",
        metrics = c(
            "FOS",
            "jaccard",
            "brayCurtis",
            "redundancyOverlap"
        ),
        topK = NULL,
        computePvalue = FALSE,
        nPermutations = 999L,
        BPPARAM = BiocParallel::SerialParam(),
        ...
    ) {
        ## 1. Validate
        valid_metrics <- c(
            "FOS",
            "jaccard",
            "brayCurtis",
            "redundancyOverlap",
            "MAAS"
        )
        method <- match.arg(method, valid_metrics)
        metrics <- match.arg(
            metrics,
            valid_metrics,
            several.ok = TRUE
        )

        if (method == "MAAS") {
            needed <- c(
                "FOS",
                "jaccard",
                "brayCurtis",
                "redundancyOverlap"
            )
            missingM <- setdiff(needed, metrics)
            if (length(missingM) > 0L) {
                cli::cli_abort(
                    "method {.val MAAS} requires {.arg metrics} \\
                    to include all four components. Missing: \\
                    {.val {missingM}}."
                )
            }
        } else if (!method %in% metrics) {
            cli::cli_abort(
                "{.arg method} ({.val {method}}) must be in \\
                {.arg metrics}."
            )
        }

        n <- length(y@Consortia)
        if (n < 1L) {
            cli::cli_abort(
                "Database CMS must contain at least one \\
                consortium."
            )
        }

        ## 2. Harmonize metabolite space: union(CM, CMS universal)
        cms_universal <- rownames(y@BinaryMatrices[[1L]])
        query_mets <- rownames(assays(x)$Binary)
        union_mets <- sort(unique(c(cms_universal, query_mets)))

        ## 3. Expand query
        xBin <- .expandMatrix(assays(x)$Binary, union_mets)

        needsWeighted <- any(
            c("brayCurtis", "redundancyOverlap") %in% metrics
        ) || method == "MAAS"

        xWeighted <- NULL
        if (needsWeighted && x@Weighted) {
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
        }

        ## 4. Expand each database consortium
        cm_names <- names(y@BinaryMatrices)
        yBins <- lapply(
            y@BinaryMatrices,
            function(m) .expandMatrix(m, union_mets)
        )

        yWeightedList <- NULL
        if (needsWeighted) {
            yWeightedList <- lapply(
                y@Consortia,
                function(cm) {
                    if (!cm@Weighted) {
                        return(NULL)
                    }
                    list(
                        Consumption = .expandMatrix(
                            assays(cm)$Consumption,
                            union_mets
                        ),
                        Production = .expandMatrix(
                            assays(cm)$Production,
                            union_mets
                        ),
                        nSpecies = .expandMatrix(
                            assays(cm)$nSpecies,
                            union_mets
                        )
                    )
                }
            )
            names(yWeightedList) <- cm_names
        }

        cli::cli_inform(
            "Searching {.val {n}} consortia using \\
            {.val {method}}."
        )

        ## 5. Score each database consortium
        scores_list <- BiocParallel::bplapply(
            seq_len(n),
            function(i) {
                yB <- yBins[[i]]
                yW <- if (is.null(yWeightedList)) {
                    NULL
                } else {
                    yWeightedList[[i]]
                }
                sc <- .computeAllScores(
                    xBin,
                    yB,
                    xWeighted,
                    yW
                )
                cov <- .computeCoverage(xBin, yB)
                sc$coverageQuery <- cov$coverageQuery
                sc$coverageReference <- cov$coverageReference
                sc
            },
            BPPARAM = BPPARAM
        )

        ## 6. Primary score per hit
        primary_scores <- vapply(
            scores_list,
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

        ## 7. Build ranking tibble
        ranking <- tibble::tibble(
            reference = cm_names,
            score = primary_scores,
            FOS = vapply(
                scores_list,
                function(s) s$FOS,
                numeric(1L)
            ),
            jaccard = vapply(
                scores_list,
                function(s) s$jaccard,
                numeric(1L)
            ),
            brayCurtis = vapply(
                scores_list,
                function(s) s$brayCurtis,
                numeric(1L)
            ),
            redundancyOverlap = vapply(
                scores_list,
                function(s) s$redundancyOverlap,
                numeric(1L)
            ),
            coverageQuery = vapply(
                scores_list,
                function(s) s$coverageQuery,
                numeric(1L)
            ),
            coverageReference = vapply(
                scores_list,
                function(s) s$coverageReference,
                numeric(1L)
            )
        ) |>
            dplyr::arrange(
                dplyr::desc(.data$score)
            )

        if (!is.null(topK)) {
            topK <- as.integer(topK)
            if (
                length(topK) != 1L ||
                    is.na(topK) ||
                    topK < 1L
            ) {
                cli::cli_abort(
                    "{.arg topK} must be a positive integer."
                )
            }
            ranking <- utils::head(ranking, topK)
        }

        ## 8. Top hit (overall best, regardless of topK)
        top_idx <- which.max(primary_scores)
        top_name <- cm_names[top_idx]
        top_score <- primary_scores[top_idx]
        top_cm <- y@Consortia[[top_idx]]

        ## 9. Similarity matrix: 1 x (topK or n)
        query_label <- if (
            length(x@Name) == 1L && !is.na(x@Name)
        ) {
            x@Name
        } else {
            "query"
        }
        sim_mat <- matrix(
            ranking$score,
            nrow = 1L,
            ncol = nrow(ranking),
            dimnames = list(query_label, ranking$reference)
        )

        ## 10. Pathway correspondences for top hit only
        correspondences <- .identifyPathwayCorrespondences(
            xBin,
            yBins[[top_idx]],
            x@Pathways,
            top_cm@Pathways
        )

        ## 11. Optional top-hit p-value
        pval <- NA_real_
        if (computePvalue) {
            if (
                method %in% c("brayCurtis", "redundancyOverlap")
            ) {
                cli::cli_warn(
                    paste0(
                        "P-value computation for ",
                        "{.val {method}} not yet ",
                        "supported. Skipping."
                    )
                )
            } else {
                top_yW <- if (is.null(yWeightedList)) {
                    NULL
                } else {
                    yWeightedList[[top_idx]]
                }
                metric_fn <- switch(
                    method,
                    FOS = .functionalOverlap,
                    jaccard = .jaccardIndex,
                    MAAS = function(xBinPerm, yBinTop) {
                        sc <- .computeAllScores(
                            xBinPerm,
                            yBinTop,
                            xWeighted,
                            top_yW
                        )
                        .computeMAAS(sc)
                    }
                )
                pval <- .computePvalue(
                    xGraph = x@Graphs[[1L]],
                    yBin = yBins[[top_idx]],
                    observed = top_score,
                    metricFn = metric_fn,
                    nPerm = nPermutations,
                    metabolites = union_mets
                )
            }
        }

        ## 12. Build and return CMA
        auto_name <- paste0(x@Name, " vs CMS (", n, " consortia)")
        ConsortiumMetabolismAlignment(
            Name = auto_name,
            Type = "search",
            Metric = method,
            QueryName = x@Name,
            ReferenceName = top_name,
            PrimaryScore = top_score,
            Pvalue = pval,
            Scores = list(
                ranking = ranking,
                top = scores_list[[top_idx]],
                coverageQuery =
                    scores_list[[top_idx]]$coverageQuery,
                coverageReference =
                    scores_list[[top_idx]]$coverageReference
            ),
            SimilarityMatrix = sim_mat,
            SharedPathways = correspondences$shared,
            UniqueQuery = correspondences$uniqueQuery,
            UniqueReference = correspondences$uniqueReference,
            Pathways = correspondences$shared,
            Metabolites = data.frame(
                met = union_mets,
                stringsAsFactors = FALSE
            ),
            Params = list(
                metrics = metrics,
                topK = topK,
                computePvalue = computePvalue,
                nPermutations = nPermutations,
                nDatabase = n,
                BPPARAM = as.character(class(BPPARAM))
            )
        )
    }
)
