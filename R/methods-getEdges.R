### The output of the two methods should be standardised.

#' @describeIn getEdges Get Edges From a \code{ConsortiumMetabolism} Object
#' @export
setMethod("getEdges", "ConsortiumMetabolism", function(object) {
  object@Edges
})

#' @describeIn getEdges Get Edges From a \code{ConsortiumMetabolismSet} Object
#' @param type Character scalar giving the type of edges to output.
#' @param quantileCutoff Numeric scalar between 0 and 1 giving the quantile
#'   threshold to use for filtering edges. For "pan-cons" and "core" types,
#'   edges above \code{1 - quantileCutoff} are returned. For "niche" and "aux"
#'   types, edges below \code{quantileCutoff} are returned. Defaults to 0.1
#'   (i.e., top/bottom 10\%).
#' @export
setMethod(
  "getEdges",
  "ConsortiumMetabolismSet",
  function(
    object,
    type = c("all", "pan-cons", "niche", "core", "aux"),
    quantileCutoff = 0.1
  ) {
    type <- match.arg(type)

    # Validate quantileCutoff parameter
    stopifnot(
      "quantileCutoff must be between 0 and 1" = quantileCutoff > 0 &&
        quantileCutoff < 1
    )

    pathways_cons <- object@Edges |>
      dplyr::reframe(
        n_cons = dplyr::n_distinct(.data$cm_name),
        .by = c("consumed", "produced", "c_ind_alig", "p_ind_alig")
      ) |>
      dplyr::arrange(dplyr::desc(.data$n_cons))

    pathways_species <- object@Edges |>
      dplyr::reframe(
        n_species = dplyr::n_distinct(.data$species),
        .by = c("consumed", "produced", "c_ind_alig", "p_ind_alig")
      )

    # Get the total n of consortia in the set
    total_cons <- length(object@Consortia)

    # Get the total n of species in the set
    total_species <- nrow(getSpecies(object))

    if (type == "all") {
      pathways_cons
    } else {
      if (type == "pan-cons") {
        # Get edges in the top (1 - quantileCutoff) of consortia
        quant <- stats::quantile(1:total_cons, p = 1 - quantileCutoff)
        dplyr::filter(pathways_cons, .data$n_cons > quant) |>
          dplyr::arrange(dplyr::desc(.data$n_cons))
      } else if (type == "niche") {
        # Get edges in the bottom quantileCutoff of consortia
        quant <- stats::quantile(1:total_cons, p = quantileCutoff)
        pathways_cons |>
          dplyr::filter(.data$n_cons < quant) |>
          dplyr::arrange(.data$n_cons)
      } else if (type == "core") {
        # Get edges in the top (1 - quantileCutoff) of species
        quant <- stats::quantile(2:total_species, p = 1 - quantileCutoff)
        pathways_species |>
          dplyr::filter(.data$n_species > quant) |>
          dplyr::arrange(dplyr::desc(.data$n_species))
      } else if (type == "aux") {
        # Get edges in the bottom quantileCutoff of species
        quant <- stats::quantile(1:total_species, p = quantileCutoff)
        pathways_species |>
          dplyr::filter(.data$n_species < quant) |>
          dplyr::arrange(.data$n_species)
      }
    }
  }
)
