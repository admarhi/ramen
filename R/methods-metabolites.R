#' @include AllClasses.R AllGenerics.R
NULL

#' @rdname metabolites
setMethod("metabolites", "ConsortiumMetabolism", function(object) {
    object@Metabolites
})

#' @rdname metabolites
setMethod("metabolites", "ConsortiumMetabolismSet",
    function(object) {
        Map(
            \(x, y) dplyr::mutate(x, consortium = y),
            lapply(object@Consortia, \(x) {
                tibble::as_tibble(
                    SummarizedExperiment::colData(x)
                )
            }),
            vapply(
                object@Consortia,
                \(x) x@Name, character(1)
            )
        ) |>
            dplyr::bind_rows() |>
            dplyr::pull("met") |>
            unique() |>
            sort()
    }
)

#' @rdname metabolites
setMethod(
    "metabolites",
    "ConsortiumMetabolismAlignment",
    function(object) {
        object@Metabolites$met
    }
)
