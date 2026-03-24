#' @include AllClasses.R AllGenerics.R
NULL

#' @param object An object of class \code{ConsortiumMetabolism}
#' @describeIn consortia Get the Community
#' @return A list with the community data.
#' @export
setMethod("consortia", "ConsortiumMetabolism", function(object) {
    object@InputData |>
        dplyr::select("met", "species", "flux")
})

#' @describeIn consortia Not applicable for alignments
#' @param object A \code{ConsortiumMetabolismAlignment} object.
#' @export
setMethod(
    "consortia",
    "ConsortiumMetabolismAlignment",
    function(object) {
        cli::cli_abort(
            "{.fun consortia} is not applicable for \\
             {.cls ConsortiumMetabolismAlignment} \\
             objects. Use {.fun scores}, \\
             {.fun sharedPathways}, or \\
             {.fun similarityMatrix} instead."
        )
    }
)
