#' @param object An object of class \code{ConsortiumMetabolism}
#' @describeIn getCo Get the Community
#' @return A list with the community data.
#' @export
setMethod("getCo", "ConsortiumMetabolism", function(object) {
  object@InputData |>
    dplyr::select("met", "species", "flux")
})
