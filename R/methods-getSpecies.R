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
      )

    if (type == "all" | type == "generalists") {
      tb |>
        dplyr::arrange(dplyr::desc(.data$n_edges))
    } else {
      tb |> dplyr::arrange(.data$n_edges)
    }
  }
)

#' @describeIn getSpecies Return Species in a Microbiome
#' @param object a \code{ConsortiumMetabolismAlignment} Object
#' @return A character vector representing the microorganisms.
setMethod("getSpecies", "ConsortiumMetabolismAlignment", function(object) {
  ### ToDo
})
