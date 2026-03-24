#' @include AllClasses.R AllGenerics.R
NULL

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

#' @rdname name-set
#' @export
setReplaceMethod("name", "ConsortiumMetabolism",
    function(object, value) {
        object@Name <- value
        validObject(object)
        object
    }
)

#' @rdname name-set
#' @export
setReplaceMethod("name", "ConsortiumMetabolismSet",
    function(object, value) {
        object@Name <- value
        validObject(object)
        object
    }
)

#' @rdname name-set
#' @export
setReplaceMethod("name", "ConsortiumMetabolismAlignment",
    function(object, value) {
        object@Name <- value
        validObject(object)
        object
    }
)
