#' @describeIn getSpecies Return Species in a Microbiome
#' @param object a \code{ConsortiumMetabolism} object
#' @return A character vector representing the microorganisms.
setMethod("getSpecies", "ConsortiumMetabolism", function(object) {
  # cat(
  #   length(unique(object@InputData$species)),
  #   " microorganisms in consortia ",
  #   object@Name,
  #   ":\n",
  #   sep = ""
  # )
  # cat(paste0("  - ", unique(object@InputData$species), collapse = "\n"), "\n")
  # invisible(unique(object@InputData$species))
  unique(object@InputData$species)
})


#' @describeIn getSpecies Return Species in a Microbiome
#' @param object a \code{ConsortiumMetabolismSet} Object
#' @param type Character scalar giving the type of species to output.
#'
#' @return A character vector representing the microorganisms.
setMethod(
  "getSpecies",
  "ConsortiumMetabolismSet",
  function(object, type = c("all", "generalists", "specialists")) {
    type <- match.arg(type)

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
      # Get the upper 15 % based on the number of edges in the object
      quant <- stats::quantile(1:total_species, p = 0.15) |> round()
      tb |> dplyr::slice_head(n = quant)
    } else if (type == "specialists") {
      # Get the upper 15 % based on the number of edges in the object
      quant <- stats::quantile(1:total_species, p = 0.15) |> round()
      tb |> dplyr::slice_tail(n = quant)
    }
  }
)

#' @describeIn getSpecies Return Species in a Microbiome
#' @param object a \code{ConsortiumMetabolismAlignment} Object
#' @return A character vector representing the microorganisms.
setMethod("getSpecies", "ConsortiumMetabolismAlignment", function(object) {
  ### ToDo
})
