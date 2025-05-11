#' @title Return Species in a Consortium
#'
#' @param object A \code{ConsortiumMetabolism} object
#' @return A character vector containing the names of species in the consortium
#' @export
setGeneric("getSpecies", function(object, ...) standardGeneric("getSpecies"))

#' @title Get Metabolites
#'
#' @description
#' Retrieves the metabolites involved in the metabolic network.
#'
#' @param object A \code{ConsortiumMetabolism} or
#' \code{ConsortiumMetabolismAlignment} object
#' @return A character vector containing the names of metabolites in the network
#' @export
setGeneric("getMet", function(object) standardGeneric("getMet"))

#' @title Get Edges From a \code{ConsortiumMetabolism} Object
#'
#' @description
#' Retrieves the edges representing metabolic interactions between species.
#'
#' @param object A \code{ConsortiumMetabolism} object
#' @return A tibble containing edge information including:
#' \itemize{
#'   \item consumed/produced metabolites
#'   \item number of species involved
#'   \item consumption/production sums
#'   \item effective consumption/production metrics
#' }
#' @export
setGeneric("getEdges", function(object, ...) standardGeneric("getEdges"))

#' @title Get the Consortia
#'
#' @description
#' Returns the consortia data in a tabular format.
#'
#' @param object A \code{ConsortiumMetabolism}, \code{ConsortiumMetabolismSet},
#' or \code{ConsortiumMetabolismAlignment} object
#' @return For \code{ConsortiumMetabolism} objects, returns a tibble with species,
#' metabolite and flux information. For \code{ConsortiumMetabolismSet} and
#' \code{ConsortiumMetabolismAlignment} objects, returns a list of such tibbles.
#' @export
setGeneric("getCo", function(object) standardGeneric("getCo"))

#' @title Align a \code{ConsortiaMetabolismSet} Object
#'
#' @description
#' Creates an alignment of multiple consortium metabolisms to identify common
#' metabolic patterns and interactions across different communities.
#'
#' @param object A \code{ConsortiumMetabolismSet} object containing multiple
#' consortium metabolisms to align
#' @return A \code{ConsortiumMetabolismAlignment} object containing the aligned
#' metabolic networks and associated metrics
#' @export
setGeneric("align", function(object) standardGeneric("align"))

#' @title Modify a \code{ConsortiaMetabolismSet} Object
#'
#' @description
#' Modifies a \code{ConsortiumMetabolismSet} by adding or removing
#' \code{ConsortiumMetabolism} objects. This allows for dynamic updating of
#' consortium sets for comparative analyses.
#'
#' @param object A \code{ConsortiumMetabolismSet} object to modify
#' @return A modified \code{ConsortiumMetabolismSet} object
#' @export
setGeneric("modify", function(object) standardGeneric("modify"))


### Improve the documentation for the methods below

#' @title Cluster a \code{ConsortiaMetabolismSet} Object
#'
#' @description
#' Modifies a \code{ConsortiumMetabolismSet} by clustering the binary matrices
#' and storing the result to the object.
#'
#' @param object A \code{ConsortiumMetabolismSet} object to modify
#' @return A modified \code{ConsortiumMetabolismSet} object
#' @export
setGeneric("cluster", function(object) standardGeneric("cluster"))

#' @title Get a cluster from a \code{ConsortiaMetabolismSet} Object
#'
#' @description
#' Gets a cluster from a \code{ConsortiumMetabolismSet} object.
#'
#' @param object A \code{ConsortiumMetabolismSet} object
#' @param node_id Numeric scalar giving the node to be extracted.
#' @param name Character scalar specifying a name for the selection.
#' @param description Character scalar describing the selection.
#'
#' @return A cluster from a \code{ConsortiumMetabolismSet} object
#' @export
setGeneric(
  "getCluster",
  function(object, node_id, name = NA_character_, description = NA_character_)
    standardGeneric("getCluster")
)


#' @title Set Names or description of ramen Objects
#'
#' @description
#' Sets a new name or description of \code{ConsortiumMetabolism},
#' \code{ConsortiumMetabolismSet}, or \code{ConsortiumMetabolismAlignment}
#' objects
#'
#' @param object A \code{ConsortiumMetabolism},
#' \code{ConsortiumMetabolismSet}, or \code{ConsortiumMetabolismAlignment}
#' object
#' @param value Character scalar specifying a name.
#'
#' @return A cluster from a \code{ConsortiumMetabolismSet} object
#' @export
setGeneric(
  "setName",
  function(object, value = NA_character_) standardGeneric("setName")
)

#' @title Set Names or description of ramen Objects
#'
#' @description
#' Sets a new name or description of \code{ConsortiumMetabolism},
#' \code{ConsortiumMetabolismSet}, or \code{ConsortiumMetabolismAlignment}
#' objects
#'
#' @param object A \code{ConsortiumMetabolism},
#' \code{ConsortiumMetabolismSet}, or \code{ConsortiumMetabolismAlignment}
#' object
#' @param value Character scalar specifying a name.
#'
#' @return A cluster from a \code{ConsortiumMetabolismSet} object
#' @export
setGeneric(
  "setDesc",
  function(object, value = NA_character_) standardGeneric("setDesc")
)
