#' @param object An object of class \code{ConsortiumMetabolism}
#' @describeIn getCo Get the Community
#' @return A list with the community data.
#' @export
setMethod("getCo", "ConsortiumMetabolism", function(object) {
    object@InputData |>
        dplyr::select("met", "species", "flux")
})

#' @describeIn getCo Not applicable for alignments
#' @param object A \code{ConsortiumMetabolismAlignment} object.
#' @export
setMethod(
    "getCo",
    "ConsortiumMetabolismAlignment",
    function(object) {
        cli::cli_abort(
            "{.fun getCo} is not applicable for \\
             {.cls ConsortiumMetabolismAlignment} objects. \\
             Use {.fun scores}, {.fun sharedPathways}, or \\
             {.fun similarityMatrix} instead."
        )
    }
)
