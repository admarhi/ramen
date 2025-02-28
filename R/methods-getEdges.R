#' @describeIn getEdges Get Edges From a \code{ConsortiumMetabolism} Object
#' @export
setMethod("getEdges", "ConsortiumMetabolism", function(object) {
  object@Edges
})
