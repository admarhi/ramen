#' @rdname setDesc
#' @export
setMethod("setDesc", "ConsortiumMetabolism", function(object, value) {
  object@Description <- value
  object
})

#' @rdname setDesc
#' @export
setMethod("setDesc", "ConsortiumMetabolismSet", function(object, value) {
  object@Description <- value
  object
})

#' @rdname setDesc
#' @export
setMethod("setDesc", "ConsortiumMetabolismAlignment", function(object, value) {
  object@Description <- value
  object
})
