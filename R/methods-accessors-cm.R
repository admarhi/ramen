#' @include AllClasses.R AllGenerics.R
NULL

## ---- CM accessor methods ---------------------------------------------------

#' @describeIn pathways Get Pathways From a
#'   \code{ConsortiumMetabolism} Object
#' @param verbose Logical scalar. If \code{FALSE}
#'   (default), returns a concise summary with columns
#'   \code{consumed}, \code{produced}, and
#'   \code{n_species}. If \code{TRUE}, returns the full
#'   pathway data including flux statistics, effective
#'   values, and per-species detail.
#' @export
setMethod(
    "pathways",
    "ConsortiumMetabolism",
    function(object, verbose = FALSE) {
        if (verbose) {
            return(object@Pathways)
        }
        object@Pathways[,
            c("consumed", "produced", "n_species")
        ]
    }
)

#' @rdname metabolites
setMethod("metabolites", "ConsortiumMetabolism", function(object) {
    object@Metabolites
})

#' @describeIn species Return Species in a Microbiome
#' @param object a \code{ConsortiumMetabolism} object
#' @return A character vector representing the microorganisms.
setMethod("species", "ConsortiumMetabolism", function(object) {
    unique(object@InputData$species)
})

#' @param object An object of class \code{ConsortiumMetabolism}
#' @describeIn consortia Get the Community
#' @return A list with the community data.
#' @export
setMethod("consortia", "ConsortiumMetabolism", function(object) {
    object@InputData |>
        dplyr::select("met", "species", "flux")
})
