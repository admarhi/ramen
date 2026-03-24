#' @include AllClasses.R AllGenerics.R
NULL

#' @rdname getMet
setMethod("getMet", "ConsortiumMetabolism", function(object) {
    object@Metabolites
})

#' @rdname getMet
setMethod("getMet", "ConsortiumMetabolismSet", function(object) {
    Map(
        \(x, y) dplyr::mutate(x, consortium = y),
        lapply(object@Consortia, \(x) {
            tibble::as_tibble(SummarizedExperiment::colData(x))
        }),
        vapply(object@Consortia, \(x) x@Name, character(1))
    ) |>
        dplyr::bind_rows() |>
        dplyr::pull("met") |>
        unique() |>
        sort()

    ### Should go in another method bc the output format differs
    # all_met |>
    #   # Change the value to TRUE to allow for easy filtering
    #   dplyr::mutate(index = TRUE) |>
    #   tidyr::pivot_wider(names_from = "consortium", values_from = "index") |>
    #   dplyr::arrange(.data$met, .locale = "C")
})

#' @rdname getMet
setMethod(
    "getMet",
    "ConsortiumMetabolismAlignment",
    function(object) {
        object@Metabolites$met
    }
)
