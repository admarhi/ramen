#' Return Species in a Microbiome
#'
#' @export
setGeneric("getSpecies", function(object) standardGeneric("getSpecies"))

#' Get Metabolites
#'
#' @param object a \code{ConsortiumMetabolism} or
#' \code{ConsortiumMetabolismAlignment} object
#'
#' @return A character vector representing the metabolites.
#' @export
setGeneric("getMet", function(object) standardGeneric("getMet"))

#' Get Edges From a \code{ConsortiumMetabolism} Object
#'
#' @param object a \code{ConsortiumMetabolism} object
#'
#' @return A list of edges in the community.
#' @export
setGeneric("getEdges", function(object) standardGeneric("getEdges"))

#' Get the Community
#'
#' Returns the community of a single \code{ConsortiumMetabolism} object in a \
#' tibble format or a list of communities in tibble format for
#' \code{ConsortiumMetabolismAlignment} objects.
#'
#' @export
setGeneric("getCo", function(object) standardGeneric("getCo"))
