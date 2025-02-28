#' @describeIn getSpecies Return Species in a Microbiome
#' @param object a \code{ConsortiumMetabolism} object
#' @return A character vector representing the microorganisms.
setMethod("getSpecies", "ConsortiumMetabolism", function(object) {
  cat(
    length(unique(object@InputData$species)),
    " microorganisms in community ",
    object@Name,
    ":\n",
    sep = ""
  )
  cat(paste0("  - ", unique(object@InputData$species), collapse = "\n"))
  invisible(unique(object@InputData$species))
})


#' @describeIn getSpecies Return Species in a Microbiome
#' @param object a \code{ConsortiumMetabolismAlignment} Object
#' @return A character vector representing the microorganisms.
setMethod("getSpecies", "ConsortiumMetabolismAlignment", function(object) {
  ### ToDo
})
