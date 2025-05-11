### The output of the two methods should be standardised.

#' @describeIn getEdges Get Edges From a \code{ConsortiumMetabolism} Object
#' @export
setMethod("getEdges", "ConsortiumMetabolism", function(object) {
  object@Edges
})

#' @describeIn getEdges Get Edges From a \code{ConsortiumMetabolismSet} Object
#' @export
setMethod(
  "getEdges",
  "ConsortiumMetabolismSet",
  function(
    object,
    type = c("all", "pan-cons", "niche", "core", "aux"),
    perc = 0.25
  ) {
    type <- match.arg(type)

    pathways <- object@Edges |>
      dplyr::reframe(
        n_cons = dplyr::n_distinct(.data$cm_name),
        .by = c("consumed", "produced", "c_ind_alig", "p_ind_alig")
      ) |>
      dplyr::arrange(dplyr::desc(.data$n_cons))

    # Get the quantiles based on length of consortia in object
    quants <- quantile(2:length(object@Consortia))

    if (type == "all") {
      pathways
    } else {
      if (type == "pan-cons") {
        dplyr::filter(pathways, .data$n_cons > quants[4])
      } else if (type == "niche") {
        dplyr::filter(pathways, .data$n_cons < quants[2])
      } else if (type == "core") {
        object@Edges |>
          dplyr::reframe(
            n_species = dplyr::n_distinct(.data$species),
            .by = c("consumed", "produced", "c_ind_alig", "p_ind_alig")
          ) |>
          dplyr::arrange(dplyr::desc(.data$n_species))
      } else if (type == "aux") {
        object@Edges |>
          dplyr::reframe(
            n_species = dplyr::n_distinct(.data$species),
            .by = c("consumed", "produced", "c_ind_alig", "p_ind_alig")
          ) |>
          dplyr::arrange(.data$n_species)
      }
    }
  }
)
