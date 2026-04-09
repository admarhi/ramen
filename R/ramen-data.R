#' @title Example Data from MiSoSoup 2024
#'
#' @description
#' List of 56 solutions from MiSoSoup for different focal
#' strains in different media.
#'
#' @name misosoup24
#' @aliases misosoup24
#' @docType data
#'
#' @return A list of 56 tibbles.
#'
#' @format A list of 56 tibbles, each with columns:
#' \describe{
#'  \item{species}{Character, species identifier.}
#'  \item{metabolite}{Character, metabolite name (BiGG
#'    identifiers).}
#'  \item{flux}{Numeric, metabolic flux (negative for
#'    consumption, positive for production).}
#' }
#'
#' @examples
#' data("misosoup24")
#' head(misosoup24[[1]])
NULL
