#' @include AllClasses.R AllGenerics.R
NULL

## ---- CMS accessor methods --------------------------------------------------

#' @describeIn pathways Get Pathways From a
#'   \code{ConsortiumMetabolismSet} Object
#' @param type Character scalar giving the type of
#'   pathways to output.
#' @param quantileCutoff Numeric scalar between 0 and 1
#'   giving the quantile threshold to use for filtering
#'   pathways. For \code{"pan-cons"} and \code{"core"}
#'   types, pathways above \code{1 - quantileCutoff}
#'   are returned. For \code{"niche"} and \code{"aux"}
#'   types, pathways below \code{quantileCutoff} are
#'   returned. Defaults to 0.1 (i.e., top/bottom 10
#'   percent).
#' @export
setMethod(
    "pathways",
    "ConsortiumMetabolismSet",
    function(
        object,
        type = c(
            "all",
            "pan-cons",
            "niche",
            "core",
            "aux"
        ),
        quantileCutoff = 0.1,
        verbose = FALSE
    ) {
        type <- match.arg(type)

        # Validate quantileCutoff parameter
        if (quantileCutoff <= 0 || quantileCutoff >= 1) {
            cli::cli_abort(
                "{.arg quantileCutoff} must be \\
                between 0 and 1 (exclusive), \\
                not {.val {quantileCutoff}}."
            )
        }

        # Build unified summary: one row per pathway
        # with both n_cons and n_species
        pathway_summary <- object@Pathways |>
            dplyr::reframe(
                n_species = dplyr::n_distinct(
                    .data$species
                ),
                n_cons = dplyr::n_distinct(
                    .data$cm_name
                ),
                .by = c("consumed", "produced")
            ) |>
            dplyr::arrange(dplyr::desc(.data$n_cons))

        total_cons <- length(object@Consortia)

        result <- if (type == "all") {
            pathway_summary
        } else if (type == "pan-cons") {
            quant <- stats::quantile(
                seq_len(total_cons),
                p = 1 - quantileCutoff
            )
            dplyr::filter(
                pathway_summary,
                .data$n_cons > quant
            ) |>
                dplyr::arrange(
                    dplyr::desc(.data$n_cons)
                )
        } else if (type == "niche") {
            quant <- stats::quantile(
                seq_len(total_cons),
                p = quantileCutoff
            )
            dplyr::filter(
                pathway_summary,
                .data$n_cons < quant
            ) |>
                dplyr::arrange(.data$n_cons)
        } else if (type == "core") {
            quant <- stats::quantile(
                pathway_summary$n_species,
                p = 1 - quantileCutoff
            )
            dplyr::filter(
                pathway_summary,
                .data$n_species > quant
            ) |>
                dplyr::arrange(
                    dplyr::desc(.data$n_species)
                )
        } else if (type == "aux") {
            quant <- stats::quantile(
                pathway_summary$n_species,
                p = quantileCutoff
            )
            dplyr::filter(
                pathway_summary,
                .data$n_species < quant
            ) |>
                dplyr::arrange(.data$n_species)
        }

        if (verbose) {
            # Join n_cons back onto full slot data
            n_cons_lookup <- pathway_summary[,
                c("consumed", "produced", "n_cons")
            ]
            dplyr::left_join(
                object@Pathways,
                n_cons_lookup,
                by = c("consumed", "produced")
            ) |>
                dplyr::semi_join(
                    result,
                    by = c("consumed", "produced")
                )
        } else {
            result[,
                c(
                    "consumed",
                    "produced",
                    "n_species",
                    "n_cons"
                )
            ]
        }
    }
)

#' @rdname overlapMatrix
#' @exportMethod overlapMatrix
setMethod(
    "overlapMatrix",
    "ConsortiumMetabolismSet",
    function(object) {
        object@OverlapMatrix
    }
)

#' @rdname metabolites
setMethod("metabolites", "ConsortiumMetabolismSet", function(object) {
    sort(object@Metabolites$met)
})

#' @describeIn consortia Get the list of
#'   \code{ConsortiumMetabolism} objects from a
#'   \code{ConsortiumMetabolismSet}
#' @param object A \code{ConsortiumMetabolismSet} object.
#' @return A named list of \code{ConsortiumMetabolism} objects.
#' @export
setMethod(
    "consortia",
    "ConsortiumMetabolismSet",
    function(object) object@Consortia
)

#' @describeIn growth Get Growth Rates From a
#'   \code{ConsortiumMetabolismSet}
#' @export
setMethod("growth", "ConsortiumMetabolismSet", function(object) {
    stats::setNames(
        lapply(object@Consortia, growth),
        vapply(
            object@Consortia,
            \(x) x@Name,
            character(1L)
        )
    )
})

#' @describeIn species Return Species in a
#'   \code{ConsortiumMetabolismSet}
#' @param object A \code{ConsortiumMetabolismSet} object.
setMethod(
    "species",
    "ConsortiumMetabolismSet",
    function(object, ...) {
        sort(unique(object@Pathways$species))
    }
)

#' @describeIn speciesSummary Species summary for a
#'   \code{ConsortiumMetabolismSet}
#' @param object A \code{ConsortiumMetabolismSet} object.
#' @export
setMethod(
    "speciesSummary",
    "ConsortiumMetabolismSet",
    function(object, ...) {
        object@Pathways |>
            dplyr::mutate(
                pathway_name = paste0(
                    .data$consumed,
                    "-",
                    .data$produced
                )
            ) |>
            dplyr::reframe(
                n_consortia = dplyr::n_distinct(
                    .data$cm_name
                ),
                n_pathways = dplyr::n_distinct(
                    .data$pathway_name
                ),
                .by = "species"
            ) |>
            dplyr::arrange(dplyr::desc(.data$n_pathways))
    }
)

#' @describeIn filterConsortia Filter consortia from a
#'   \code{ConsortiumMetabolismSet}
#' @param object A \code{ConsortiumMetabolismSet} object.
#' @param i Integer vector, character vector of consortium
#'   names, or logical vector.
#' @importFrom S4Vectors metadata
#' @export
setMethod(
    "filterConsortia",
    "ConsortiumMetabolismSet",
    function(object, i) {
        n <- length(object@Consortia)
        cm_names <- vapply(
            object@Consortia,
            \(x) x@Name,
            character(1L)
        )

        idx <- if (is.logical(i)) {
            if (length(i) != n) {
                cli::cli_abort(
                    "Logical {.arg i} must have length \\
                    {n}, not {length(i)}."
                )
            }
            which(i)
        } else if (is.character(i)) {
            not_found <- setdiff(i, cm_names)
            if (length(not_found) > 0L) {
                cli::cli_abort(
                    "Consortium name{?s} not found: \\
                    {.val {not_found}}."
                )
            }
            match(i, cm_names)
        } else if (is.numeric(i)) {
            idx_int <- as.integer(i)
            oob <- idx_int[idx_int < 1L | idx_int > n]
            if (length(oob) > 0L) {
                cli::cli_abort(
                    "Index {.val {oob}} out of bounds \\
                    (1 to {n})."
                )
            }
            idx_int
        } else {
            cli::cli_abort(
                "{.arg i} must be integer, character, \\
                or logical, not {.cls {class(i)}}."
            )
        }

        selected <- object@Consortia[idx]
        lnk <- metadata(object)$linkage
        if (is.null(lnk)) {
            lnk <- "complete"
        }
        ConsortiumMetabolismSet(
            selected,
            name = object@Name,
            desc = object@Description,
            linkage = lnk,
            verbose = FALSE
        )
    }
)
