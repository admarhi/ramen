#' @include AllClasses.R AllGenerics.R

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

#' @describeIn prevalence Edge prevalence from a
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

#' @describeIn consensusEdges Consensus edges from a
#'   multiple [ConsortiumMetabolismAlignment]
#' @export
setMethod(
    "consensusEdges",
    "ConsortiumMetabolismAlignment",
    function(object) {
        if (object@Type != "multiple") {
            cli::cli_abort(
                "{.fun consensusEdges} is only available \\
                 for multiple alignments."
            )
        }
        object@ConsensusEdges
    }
)
