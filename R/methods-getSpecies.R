#' @describeIn getSpecies Return Species in a Microbiome
#' @param object a \code{ConsortiumMetabolism} object
#' @return A character vector representing the microorganisms.
setMethod("getSpecies", "ConsortiumMetabolism", function(object) {
  unique(object@InputData$species)
})


#' @describeIn getSpecies Return Species in a Microbiome
#' @param object a \code{ConsortiumMetabolismSet} Object
#' @param type Character scalar giving the type of species to output.
#' @param quantileCutoff Numeric scalar between 0 and 1 specifying the fraction
#'   of species to return when \code{type} is "generalists" or "specialists".
#'   For "generalists", the top \code{quantileCutoff} fraction of species with
#'   the most edges is returned. For "specialists", the bottom
#'   \code{quantileCutoff} fraction with the fewest edges is returned.
#'   Defaults to 0.15 (i.e., 15\%). Ignored when \code{type = "all"}.
#'
#' @return A character vector representing the microorganisms.
setMethod(
  "getSpecies",
  "ConsortiumMetabolismSet",
  function(
    object,
    type = c("all", "generalists", "specialists"),
    quantileCutoff = 0.15
  ) {
    type <- match.arg(type)

    # Validate quantileCutoff parameter
    stopifnot(
      "quantileCutoff must be between 0 and 1" = quantileCutoff > 0 &&
        quantileCutoff < 1
    )

    tb <- object@Edges |>
      dplyr::mutate(edge_name = paste0(.data$consumed, "-", .data$produced)) |>
      dplyr::reframe(
        n_edges = dplyr::n_distinct(.data$edge_name),
        .by = "species"
      ) |>
      dplyr::arrange(dplyr::desc(.data$n_edges))

    total_species <- length(unique(tb$species))
    if (type == "all") {
      tb
    } else if (type == "generalists") {
      # Get the top quantileCutoff fraction of species
      n_species_to_return <- ceiling(total_species * quantileCutoff)
      tb |> dplyr::slice_head(n = n_species_to_return)
    } else if (type == "specialists") {
      # Get the bottom quantileCutoff fraction of species
      n_species_to_return <- ceiling(total_species * quantileCutoff)
      tb |> dplyr::slice_tail(n = n_species_to_return)
    }
  }
)

#' @describeIn getSpecies Return Species in a Microbiome
#' @param object a \code{ConsortiumMetabolismAlignment} Object
#' @return A character vector representing the microorganisms.
setMethod("getSpecies", "ConsortiumMetabolismAlignment", function(object) {
  ### ToDo
})
