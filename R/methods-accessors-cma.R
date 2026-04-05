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

