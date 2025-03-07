#' Return Species in a Consortium
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

#' Get the Consortia
#'
#' Returns the consortia of a single \code{ConsortiumMetabolism} object in a
#' `tibble` format or a list of consortias in `tibble` format for
#' \code{ConsortiumMetabolismSet} and \code{ConsortiumMetabolismAlignment}
#' objects.
#'
#' @export
setGeneric("getCo", function(object) standardGeneric("getCo"))

#' @title Align a \code{ConsortiaMetabolismSet} Object
#'
#' @description
#' A short description...
#'
#' @export
setGeneric("align", function(object) standardGeneric("align"))

#' @title Modify a \code{ConsortiaMetabolismSet} Object
#'
#' @description
#' To be used to remove or add \code{ConsortiumMetabolism} object to or from a
#'
#' @export
setGeneric("modify", function(object) standardGeneric("modify"))
