#' @examples
#' # Create a synthetic consortium
#' cm <- synCM("example", n_species = 3, max_met = 5)
#' cm
#'
#' @keywords internal
"_PACKAGE"

## usethis namespace: start
#' @importFrom dplyr inner_join
#' @importFrom dplyr left_join
#' @importFrom dplyr mutate
#' @importFrom dplyr reframe
#' @importFrom dplyr rename
#' @importFrom dplyr select
#' @importFrom Matrix sparseMatrix
#' @importFrom methods show
#' @importFrom rlang .data
#' @importFrom SummarizedExperiment assay
#' @importFrom SummarizedExperiment assays
#' @importFrom SummarizedExperiment metadata<-
#' @importFrom tibble rowid_to_column
#' @importFrom tibble tibble
#' @importFrom tibble tribble
#' @importFrom tidyr nest
#' @importFrom TreeSummarizedExperiment TreeSummarizedExperiment
## usethis namespace: end
NULL
