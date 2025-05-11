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
#' @return A character vector representing the microorganisms.
setMethod(
  "getSpecies",
  "ConsortiumMetabolismSet",
  function(object, type = c("all", "generalists", "specialists")) {
    type <- match.arg(type)

    object@Edges |>
      dplyr::mutate(edge_name = paste0(.data$consumed, "-", .data$produced)) |>
      dplyr::reframe(
        n_edges = dplyr::n_distinct(.data$edge_name),
        .by = "species"
      ) |>
      dplyr::arrange(dplyr::desc(n_edges))
  }
)

#' @describeIn getSpecies Return Species in a Microbiome
#' @param object a \code{ConsortiumMetabolismAlignment} Object
#' @return A character vector representing the microorganisms.
setMethod("getSpecies", "ConsortiumMetabolismAlignment", function(object) {
  ### ToDo
})
