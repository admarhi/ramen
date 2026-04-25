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
setMethod(
    "metabolites",
    "ConsortiumMetabolism",
    function(
        object,
        species = NULL,
        direction = c("all", "consumed", "produced")
    ) {
        direction <- match.arg(direction)
        if (is.null(species) && direction == "all") {
            return(object@Metabolites)
        }
        sp_pw <- tidyr::unnest(object@Pathways, "data")
        if (!is.null(species)) {
            if (!is.character(species) || length(species) != 1L) {
                cli::cli_abort(
                    "{.arg species} must be a length-1 character scalar."
                )
            }
            if (!species %in% sp_pw$species) {
                cli::cli_abort(
                    "Species {.val {species}} not found in this consortium."
                )
            }
            sp_pw <- sp_pw[sp_pw$species == species, , drop = FALSE]
        }
        mets <- switch(
            direction,
            all = unique(c(sp_pw$consumed, sp_pw$produced)),
            consumed = unique(sp_pw$consumed),
            produced = unique(sp_pw$produced)
        )
        sort(mets)
    }
)

#' @describeIn species Return Species in a
#'   \code{ConsortiumMetabolism}
#' @param object A \code{ConsortiumMetabolism} object.
setMethod("species", "ConsortiumMetabolism", function(object) {
    unique(object@InputData$species)
})

#' @describeIn speciesSummary Species summary for a
#'   \code{ConsortiumMetabolism}
#' @param object A \code{ConsortiumMetabolism} object.
#' @export
setMethod(
    "speciesSummary",
    "ConsortiumMetabolism",
    function(object, ...) {
        object@Pathways |>
            tidyr::unnest("data") |>
            dplyr::reframe(
                n_pathways = dplyr::n_distinct(
                    paste0(
                        .data$consumed,
                        "-",
                        .data$produced
                    )
                ),
                n_consumed = dplyr::n_distinct(
                    .data$consumed
                ),
                n_produced = dplyr::n_distinct(
                    .data$produced
                ),
                .by = "species"
            ) |>
            dplyr::arrange(dplyr::desc(.data$n_pathways))
    }
)

#' @param object An object of class \code{ConsortiumMetabolism}
#' @describeIn consortia Get the Community
#' @return A list with the community data.
#' @export
setMethod("consortia", "ConsortiumMetabolism", function(object) {
    object@InputData |>
        dplyr::select("met", "species", "flux")
})

#' @rdname growth
#' @importFrom S4Vectors metadata
#' @export
setMethod("growth", "ConsortiumMetabolism", function(object) {
    metadata(object)$growth
})
