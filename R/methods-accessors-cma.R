#' @include AllClasses.R AllGenerics.R
NULL

## ---- CMA accessor methods ------------------------------------------------

#' @describeIn scores Scores from a
#'   [ConsortiumMetabolismAlignment]
#' @export
setMethod(
    "scores",
    "ConsortiumMetabolismAlignment",
    function(object) object@Scores
)

#' @describeIn sharedPathways Shared pathways from a
#'   pairwise [ConsortiumMetabolismAlignment]
#' @export
setMethod(
    "sharedPathways",
    "ConsortiumMetabolismAlignment",
    function(object) {
        if (object@Type != "pairwise") {
            cli::cli_abort(
                "{.fun sharedPathways} is only available \\
                 for pairwise alignments."
            )
        }
        object@SharedPathways
    }
)

#' @describeIn uniquePathways Unique pathways from a
#'   pairwise [ConsortiumMetabolismAlignment]
#' @export
setMethod(
    "uniquePathways",
    "ConsortiumMetabolismAlignment",
    function(object) {
        if (object@Type != "pairwise") {
            cli::cli_abort(
                "{.fun uniquePathways} is only available \\
                 for pairwise alignments."
            )
        }
        list(
            query = object@UniqueQuery,
            reference = object@UniqueReference
        )
    }
)

#' @describeIn similarityMatrix Similarity matrix from a
#'   multiple [ConsortiumMetabolismAlignment]
#' @export
setMethod(
    "similarityMatrix",
    "ConsortiumMetabolismAlignment",
    function(object) {
        if (object@Type != "multiple") {
            cli::cli_abort(
                "{.fun similarityMatrix} is only available \\
                 for multiple alignments."
            )
        }
        object@SimilarityMatrix
    }
)

#' @describeIn prevalence Pathway prevalence from a
#'   multiple [ConsortiumMetabolismAlignment]
#' @export
setMethod(
    "prevalence",
    "ConsortiumMetabolismAlignment",
    function(object) {
        if (object@Type != "multiple") {
            cli::cli_abort(
                "{.fun prevalence} is only available \\
                 for multiple alignments."
            )
        }
        object@Prevalence
    }
)

#' @describeIn consensusPathways Consensus pathways from
#'   a multiple [ConsortiumMetabolismAlignment]
#' @export
setMethod(
    "consensusPathways",
    "ConsortiumMetabolismAlignment",
    function(object) {
        if (object@Type != "multiple") {
            cli::cli_abort(
                "{.fun consensusPathways} is only \\
                 available for multiple alignments."
            )
        }
        object@ConsensusPathways
    }
)
