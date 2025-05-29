#' Import and Process Misosoup Data
#'
#' Processes raw misosoup data into a structured format containing consortia,
#' media, and growth information. The function handles data cleaning,
#' transformation, and organization of metabolic flux data.
#'
#' @param data A nested list containing misosoup simulation results.
#' The structure should be data[[substrate]][[focal_strain]] where each element
#' contains solution data.
#'
#' @return A list containing three tibbles:
#'   \item{consortia}{Metabolic flux data for each species in the consortium}
#'   \item{media}{Media composition data}
#'   \item{growth}{Growth information for each solution}
#'
#' @seealso \code{\link{overviewMisosoup}} for a summary of the input data
#'
#' @export
importMisosoup <- function(data) {
  tb_import <- overviewMisosoup(data) |>
    dplyr::filter(n_cons != n_zero_growth)

  if (any(tb_import$n_zero_growth > 0)) stop("Zero growth solution detected")

  tb <- purrr::map2(
    tb_import$substrate,
    tb_import$focal_strain,
    \(x, y)
      data[[x]][[y]] |>
        purrr::imap(
          \(z, idx)
            tibble::as_tibble(z$solution) |>
              tidyr::pivot_longer(
                cols = dplyr::everything(),
                names_to = "rxn",
                values_to = "flux"
              ) |>
              dplyr::mutate(
                substrate = x,
                focal_strain = y,
                solution = idx
              )
        ) |>
        dplyr::bind_rows(),
    .progress = TRUE
  ) |>
    dplyr::bind_rows() |>
    dplyr::relocate("rxn", "flux", .after = "solution") |>
    dplyr::mutate(
      cons_id = paste(
        .data$substrate,
        .data$focal_strain,
        .data$solution,
        sep = "_"
      )
    ) |>
    dplyr::relocate("cons_id")

  # Filter growth information
  growth <- tb |>
    dplyr::filter(stringr::str_detect(.data$rxn, "[G,g]rowth")) |>
    dplyr::rename(growth = "rxn") |>
    dplyr::mutate(
      growth = stringr::str_remove(.data$growth, "^Growth_|_growth$")
    )

  # Clean and split the rxn column
  tb <- tb |>
    dplyr::filter(!stringr::str_detect(.data$rxn, "[G,g]rowth")) |>
    tidyr::separate_wider_delim(
      cols = "rxn",
      delim = "_e_",
      names = c("met", "species"),
      too_few = "align_start"
    ) |>
    dplyr::mutate(
      met = stringr::str_remove_all(.data$met, "^R_EX_|_e$"),
      species = stringr::str_remove(.data$species, "_i$")
    )

  # Filter NA values in species as these entries are the media
  consortia <- dplyr::filter(tb, !is.na(.data$species))

  # Compute the media tibble
  media <- tb |>
    dplyr::filter(is.na(.data$species)) |>
    dplyr::select(1:5, "flux")

  list(
    consortia = consortia,
    media = media,
    growth = growth
  )
}

#' Overview of Misosoup Data Structure
#'
#' Provides a summary of the misosoup data structure, including the number of
#' consortia and zero-growth solutions for each substrate and focal strain
#' combination.
#'
#' @inheritParams importMisosoup
#'
#' @return A tibble with columns:
#'   \item{substrate}{Substrate identifier}
#'   \item{focal_strain}{Focal strain identifier}
#'   \item{n_cons}{Number of consortia solutions}
#'   \item{n_zero_growth}{Number of zero-growth solutions}
#'
#' @seealso \code{\link{importMisosoup}} for processing the full data
#'
#' @export
overviewMisosoup <- function(data) {
  purrr::map(data, \(x) names(x)) |>
    tibble::enframe(name = "substrate", value = "focal_strain") |>
    tidyr::unnest_longer(col = "focal_strain") |>
    dplyr::mutate(
      n_cons = purrr::map2_dbl(
        .data$substrate,
        .data$focal_strain,
        \(x, y) length(data[[x]][[y]])
      ),
      n_zero_growth = purrr::map2(
        .data$substrate,
        .data$focal_strain,
        \(x, y)
          data[[x]][[y]] |>
            purrr::map_lgl(
              \(z) length(z) == 1 && z[[1]] == 0
            )
      ) |>
        purrr::map_dbl(\(x) sum(x))
    )
}
