#' @include AllClasses.R AllGenerics.R
NULL

## ---- name accessor ---------------------------------------------------------

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

## ---- description accessor --------------------------------------------------

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
