#' @include AllClasses.R AllGenerics.R
NULL

#' @rdname description
#' @export
setMethod("description", "ConsortiumMetabolism",
    function(object) object@Description
)

#' @rdname description
#' @export
setMethod("description", "ConsortiumMetabolismSet",
    function(object) object@Description
)

#' @rdname description
#' @export
setMethod("description", "ConsortiumMetabolismAlignment",
    function(object) object@Description
)

#' @rdname description
#' @export
setReplaceMethod("description", "ConsortiumMetabolism",
    function(object, value) {
        object@Description <- value
        validObject(object)
        object
    }
)

#' @rdname description
#' @export
setReplaceMethod("description", "ConsortiumMetabolismSet",
    function(object, value) {
        object@Description <- value
        validObject(object)
        object
    }
)

#' @rdname description
#' @export
setReplaceMethod("description", "ConsortiumMetabolismAlignment",
    function(object, value) {
        object@Description <- value
        validObject(object)
        object
    }
)
