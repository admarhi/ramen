### The output of the two methods should be standardised.

#' @describeIn getEdges Get Edges From a \code{ConsortiumMetabolism} Object
#' @export
setMethod("getEdges", "ConsortiumMetabolism", function(object) {
  object@Edges
})

#' @describeIn getEdges Get Edges From a \code{ConsortiumMetabolismSet} Object
#' @param type Character scalar giving the type of edges to output.
#' @param perc Numeric scalar giving the percentage to use for the quantile
#' calculation.
#' @export
setMethod(
  "getEdges",
  "ConsortiumMetabolismSet",
  function(
    object,
    type = c("all", "pan-cons", "niche", "core", "aux"),
    perc = 0.1
  ) {
    type <- match.arg(type)

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
        # Get the upper 10 % based on the number of consortia in the object
        quant <- stats::quantile(1:total_cons, p = 0.9)
        # Returns all pathways that appear in all consortia
        dplyr::filter(pathways_cons, .data$n_cons > quant) |>
          dplyr::arrange(dplyr::desc(.data$n_cons))
      } else if (type == "niche") {
        # Get the lower quantile based on the number of consortia in the object
        quant <- stats::quantile(1:total_cons, p = 0.1)
        pathways_cons |>
          dplyr::filter(.data$n_cons < quant) |>
          dplyr::arrange(.data$n_cons)
      } else if (type == "core") {
        # Get the upper quantile based on number of species in the object
        quant <- stats::quantile(2:total_species, p = 0.8)
        pathways_species |>
          dplyr::filter(.data$n_species > quant) |>
          dplyr::arrange(dplyr::desc(.data$n_species))
      } else if (type == "aux") {
        # Get the lower quantile based on number of species in the object
        quant <- stats::quantile(1:total_species, p = 0.2)
        pathways_species |>
          dplyr::filter(.data$n_species < quant) |>
          dplyr::arrange(.data$n_species)
      }
    }
  }
)
