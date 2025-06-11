#' @title Return Species in a Consortium
#'
#' @param object A \code{ConsortiumMetabolism} object
#' @param ... Object specific arguments.
#'
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
#' The argument `type` can be used to return only specfic types of Edges from
#' a `ConsortiumMetabolismSet` type object.
#' \itemize{
#'   \item `all` will return all edges
#'   \item `pan-cons` will return only edges that exist in all consortia
#'   \item `niche` will return niche consortia. A niche is defined a
#' }
#'
#' @param object A \code{ConsortiumMetabolism} object
#' @param ... Object specific arguments.
#'
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
#' @param name Character scalar giving name of the alignment, if `NULL` inherits
#' from the \code{ConsortiumMetabolismSet} object.
#' @return A \code{ConsortiumMetabolismAlignment} object containing the aligned
#' metabolic networks and associated metrics
#' @export
setGeneric("align", function(object, name) standardGeneric("align"))

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
  function(object, node_id, name = NA_character_, description = NA_character_) {
    standardGeneric("getCluster")
  }
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

#' @title Get Functional Groups
#'
#' @description
#' Calculates and returns functional groups based on metabolic reactions.
#' For \code{ConsortiumMetabolismSet} objects, this involves analyzing shared
#' reactions across species to identify clusters of species with similar
#' metabolic capabilities.
#'
#' @details
#' This method is currently implemented for \code{ConsortiumMetabolismSet}
#' objects. Future versions will extend functionality to
#' \code{ConsortiumMetabolismAlignment} objects to allow for comparative
#' functional group analysis across different alignments.
#'
#' @param object A \code{ConsortiumMetabolismSet} object.
#' @param k An integer scalar specifying the number of clusters to color in the
#'   dendrogram.
#' @param ... Additional arguments to be passed to specific methods.
#'
#' @return A dendrogram object representing the hierarchical clustering of
#' species into functional groups. The plot of the dendrogram is also displayed.
#'
#' @export
setGeneric(
  "getFunctionalGroups",
  function(object, k = 4, ...) standardGeneric("getFunctionalGroups")
)
