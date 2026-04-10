#' Import and Process Misosoup Data
#'
#' Processes raw misosoup data into a structured format containing consortia,
#' media, and growth information. The function handles data cleaning,
#' transformation, and organization of metabolic flux data.
#'
#' @param data A nested list containing misosoup simulation results.
#' The structure should be data[\[substrate\]][\[focal_strain\]] where each
#' element contains solution data.
#'
#' @return A list containing three tibbles:
#'   \item{consortia}{Metabolic flux data for each species in the consortium}
#'   \item{media}{Media composition data}
#'   \item{growth}{Growth information for each solution}
#'
#' @seealso \code{\link{overviewMisosoup}} for a summary of the input data
#'
#' @export
#'
#' @examples
#' \donttest{
#' # Requires raw MiSoSoup YAML data
#' raw <- yaml::read_yaml("misosoup_output.yaml")
#' result <- importMisosoup(raw)
#' str(result, max.level = 1)
#' }
importMisosoup <- function(data) {
    tb_import <- overviewMisosoup(data) |>
        dplyr::filter(.data$n_cons != .data$n_zero_growth)

    if (any(tb_import$n_zero_growth > 0)) {
        cli::cli_abort("Zero-growth solution detected in input data.")
    }

    tb <- Map(
        \(x, y) {
            entries <- data[[x]][[y]]
            Map(
                \(z, idx) {
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
                },
                entries,
                seq_along(entries)
            ) |>
                dplyr::bind_rows()
        },
        tb_import$substrate,
        tb_import$focal_strain
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
        dplyr::filter(grepl("[Gg]rowth", .data$rxn)) |>
        dplyr::rename(growth = "rxn") |>
        dplyr::mutate(
            growth = sub("^Growth_", "", .data$growth)
        )

    # Clean and split the rxn column
    tb <- tb |>
        dplyr::filter(!grepl("[Gg]rowth", .data$rxn)) |>
        tidyr::separate_wider_delim(
            cols = "rxn",
            delim = "_e_",
            names = c("met", "species"),
            too_few = "align_start"
        ) |>
        dplyr::mutate(
            met = gsub("^R_EX_|_e$", "", .data$met),
            species = sub("_i$", "", .data$species)
        )

    # Filter NA values in species as these entries are the media
    consortia <- dplyr::filter(tb, !is.na(.data$species))

    # Compute the media tibble
    media <- tb |>
        dplyr::filter(is.na(.data$species)) |>
        dplyr::select(
            "cons_id",
            "substrate",
            "focal_strain",
            "solution",
            "met",
            "flux"
        )

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
#'
#' @examples
#' \donttest{
#' # Requires raw MiSoSoup YAML data
#' raw <- yaml::read_yaml("misosoup_output.yaml")
#' overviewMisosoup(raw)
#' }
overviewMisosoup <- function(data) {
    lapply(data, \(x) names(x)) |>
        tibble::enframe(
            name = "substrate",
            value = "focal_strain"
        ) |>
        tidyr::unnest_longer(col = "focal_strain") |>
        dplyr::mutate(
            n_cons = vapply(
                seq_len(dplyr::n()),
                \(i) {
                    x <- .data$substrate[[i]]
                    y <- .data$focal_strain[[i]]
                    length(data[[x]][[y]])
                },
                integer(1)
            ),
            n_zero_growth = vapply(
                seq_len(dplyr::n()),
                \(i) {
                    x <- .data$substrate[[i]]
                    y <- .data$focal_strain[[i]]
                    sum(vapply(
                        data[[x]][[y]],
                        \(z) is.null(z[["solution"]]),
                        logical(1)
                    ))
                },
                numeric(1)
            )
        )
}
