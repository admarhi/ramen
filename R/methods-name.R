#' @include AllClasses.R AllGenerics.R
NULL

#' @rdname name
#' @export
setMethod("name", "ConsortiumMetabolism",
    function(object) object@Name
)

#' @rdname name
#' @export
setMethod("name", "ConsortiumMetabolismSet",
    function(object) object@Name
)

#' @rdname name
#' @export
setMethod("name", "ConsortiumMetabolismAlignment",
    function(object) object@Name
)

#' @rdname name
#' @export
setReplaceMethod("name", "ConsortiumMetabolism",
    function(object, value) {
        object@Name <- value
        validObject(object)
        object
    }
)

#' @rdname name
#' @export
setReplaceMethod("name", "ConsortiumMetabolismSet",
    function(object, value) {
        object@Name <- value
        validObject(object)
        object
    }
)

#' @rdname name
#' @export
setReplaceMethod("name", "ConsortiumMetabolismAlignment",
    function(object, value) {
        object@Name <- value
        validObject(object)
        object
    }
)
