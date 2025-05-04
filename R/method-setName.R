#' @rdname setName
#' @export
setMethod("setName", "ConsortiumMetabolism", function(object, value) {
  object@Name <- value
  object
})

#' @rdname setName
#' @export
setMethod("setName", "ConsortiumMetabolismSet", function(object, value) {
  object@Name <- value
  object
})

#' @rdname setName
#' @export
setMethod("setName", "ConsortiumMetabolismAlignment", function(object, value) {
  object@Name <- value
  object
})
