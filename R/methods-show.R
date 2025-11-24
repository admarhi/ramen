#' Show Method for \code{ConsortiumMetabolism} Object
#'
#' @param object An object of class \code{ConsortiumMetabolism}
#' @export
setMethod("show", "ConsortiumMetabolism", function(object) {
  stringr::str_glue(
    "\n{object@Name}: ConsortiumMetabolism Object\n",
    "{ifelse(object@Weighted, 'Weighted', 'Unweighted')} ",
    "metabolic network with {length(object@Metabolites)} metabolites.\n\n"
  ) |>
    cat()
})

#' Show Method for \code{ConsortiumMetabolismSet} Object
#'
#' @param object An object of class \code{ConsortiumMetabolismSet}
#' @export
setMethod("show", "ConsortiumMetabolismSet", function(object) {
  stringr::str_glue(
    "\n{object@Name}: ConsortiumMetabolismSet Object\n",
    "containing {length(object@Consortia)} consortia.\n",
    "Description: {object@Description}\n\n"
  ) |>
    cat()
})

#' Show method for \code{ConsortiumMetabolismAlignment} Objects
#'
#' @param object a \code{ConsortiumMetabolismAlignment} object.
#' @export
setMethod("show", "ConsortiumMetabolismAlignment", function(object) {
  stringr::str_glue(
    "Functional Alignment of {length(object@Graphs)} ConsortiumMetabolism\n\n"
  ) |>
    cat()
})
