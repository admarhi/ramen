#' @include AllClasses.R AllGenerics.R
NULL

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

#' @rdname description-set
#' @export
setReplaceMethod("description", "ConsortiumMetabolism",
    function(object, value) {
        object@Description <- value
        validObject(object)
        object
    }
)

#' @rdname description-set
#' @export
setReplaceMethod("description", "ConsortiumMetabolismSet",
    function(object, value) {
        object@Description <- value
        validObject(object)
        object
    }
)

#' @rdname description-set
#' @export
setReplaceMethod("description", "ConsortiumMetabolismAlignment",
    function(object, value) {
        object@Description <- value
        validObject(object)
        object
    }
)
