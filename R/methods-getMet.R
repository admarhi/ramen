#' @rdname getMet
setMethod("getMet", "ConsortiumMetabolism", function(object) {
  object@Metabolites
})


#' @rdname getMet
setMethod("getMet", "ConsortiumMetabolismAlignment", function(object) {
  ### ToDo
})
