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
        total_species <- nrow(species(object))

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
                2:total_species,
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
                seq_len(total_species),
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

#' @rdname metabolites
setMethod("metabolites", "ConsortiumMetabolismSet", function(object) {
    Map(
        \(x, y) dplyr::mutate(x, consortium = y),
        lapply(object@Consortia, \(x) {
            tibble::as_tibble(
                SummarizedExperiment::colData(x)
            )
        }),
        vapply(
            object@Consortia,
            \(x) x@Name,
            character(1)
        )
    ) |>
        dplyr::bind_rows() |>
        dplyr::pull("met") |>
        unique() |>
        sort()
})

#' @describeIn species Return Species in a Microbiome
#' @param object a \code{ConsortiumMetabolismSet} Object
#' @param type Character scalar giving the type of species to output.
#' @param quantileCutoff Numeric scalar between 0 and 1 specifying
#'   the fraction of species to return when \code{type} is
#'   "generalists" or "specialists".
#'   For "generalists", the top \code{quantileCutoff} fraction of
#'   species with the most pathways is returned. For
#'   "specialists", the bottom \code{quantileCutoff} fraction
#'   with the fewest pathways is returned.
#'   Defaults to 0.15 (i.e., 15 percent). Ignored when
#'   \code{type = "all"}.
#'
#' @return A character vector representing the microorganisms.
setMethod(
    "species",
    "ConsortiumMetabolismSet",
    function(
        object,
        type = c("all", "generalists", "specialists"),
        quantileCutoff = 0.15
    ) {
        type <- match.arg(type)

        # Validate quantileCutoff parameter
        if (quantileCutoff <= 0 || quantileCutoff >= 1) {
            cli::cli_abort(
                "{.arg quantileCutoff} must be between \\
                0 and 1 (exclusive), not \\
                {.val {quantileCutoff}}."
            )
        }

        tb <- object@Pathways |>
            dplyr::mutate(
                pathway_name = paste0(
                    .data$consumed,
                    "-",
                    .data$produced
                )
            ) |>
            dplyr::reframe(
                n_pathways = dplyr::n_distinct(
                    .data$pathway_name
                ),
                .by = "species"
            ) |>
            dplyr::arrange(dplyr::desc(.data$n_pathways))

        total_species <- length(unique(tb$species))
        if (type == "all") {
            tb
        } else if (type == "generalists") {
            # Get the top quantileCutoff fraction of species
            n_species_to_return <- ceiling(
                total_species * quantileCutoff
            )
            tb |> dplyr::slice_head(n = n_species_to_return)
        } else if (type == "specialists") {
            # Get the bottom quantileCutoff fraction of species
            n_species_to_return <- ceiling(
                total_species * quantileCutoff
            )
            tb |> dplyr::slice_tail(n = n_species_to_return)
        }
    }
)
