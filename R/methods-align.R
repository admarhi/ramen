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
        cli::cli_abort(
            "Pairwise alignment not yet implemented (Phase 1)."
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
