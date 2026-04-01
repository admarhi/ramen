#' Generate Synthetic Consortium Metabolism
#'
#' Creates a synthetic metabolic network with random species-metabolite
#' interactions and flux values. Useful for testing and demonstrating
#' `ramen` package functionality.
#'
#' Each species is assigned a random subset of metabolites (between 2
#' and \code{max_met}) with normally distributed flux values. By
#' default, species are guaranteed to have both positive and negative
#' fluxes (i.e., no dead ends).
#'
#' @param name Character string. Name for the consortium.
#' @param n_species Integer. Number of species in the consortium.
#' @param max_met Integer. Size of the metabolite pool to sample
#'     from.
#' @param scale_fac Integer. Multiplier for the initial species name
#'     pool from which \code{n_species} are sampled. Defaults to 2.
#' @param seed Integer or \code{FALSE}. Random seed for
#'     reproducibility. If \code{FALSE} (default), a random seed is
#'     chosen.
#' @param dead_ends Logical. If \code{FALSE} (default), ensures each
#'     species has both consumed and produced metabolites by flipping
#'     one flux value when all fluxes share the same sign.
#' @param cm Logical. If \code{TRUE} (default), returns a
#'     \code{\linkS4class{ConsortiumMetabolism}} object. If
#'     \code{FALSE}, returns the raw edge list as a tibble.
#'
#' @return A \code{\linkS4class{ConsortiumMetabolism}} object when
#'     \code{cm = TRUE}, or a \code{\link[tibble]{tibble}} with
#'     columns \code{species}, \code{metabolites}, and \code{fluxes}
#'     when \code{cm = FALSE}.
#'
#' @export
#'
#' @importFrom tibble tibble
#'
#' @examples
#' synCM("Ex. Community", n_species = 5, max_met = 10)
#'
#' # Return raw edge list instead
#' synCM("Test", n_species = 3, max_met = 8, cm = FALSE)
synCM <- function(
    name,
    n_species,
    max_met,
    scale_fac = 2,
    seed = FALSE,
    dead_ends = FALSE,
    cm = TRUE
) {
    if (!is.character(name) || length(name) != 1L) {
        cli::cli_abort(
            "{.arg name} must be a single character string."
        )
    }
    if (!is.numeric(n_species) || length(n_species) != 1L ||
        n_species < 1L) {
        cli::cli_abort(
            "{.arg n_species} must be a positive integer,
            not {.val {n_species}}."
        )
    }
    if (!is.numeric(max_met) || length(max_met) != 1L ||
        max_met < 2L) {
        cli::cli_abort(
            "{.arg max_met} must be an integer >= 2,
            not {.val {max_met}}."
        )
    }
    if (!is.numeric(scale_fac) || length(scale_fac) != 1L ||
        scale_fac < 1L) {
        cli::cli_abort(
            "{.arg scale_fac} must be a positive integer,
            not {.val {scale_fac}}."
        )
    }
    .rNames <- function(n = 5000) {
        a <- do.call(paste0, replicate(3, sample(LETTERS, n, TRUE), FALSE))
        paste0(
            a,
            sprintf("%03d", sample(9999, n, TRUE)),
            sample(LETTERS, n, TRUE)
        )
    }

    species_names <- .rNames(n_species * scale_fac)
    met_vec <- paste0("met", 1:max_met)

    if (!seed) {
        seed <- sample(1:1000, 1)
    }
    set.seed(seed)

    # Get sample of species for the community
    species_names <- sample(species_names, size = n_species, replace = FALSE)

    per_species <- lapply(species_names, function(sp) {
        n_mets <- sample(2:max_met, 1)
        species_mets <- sample(met_vec, n_mets)
        flux_vals <- stats::rnorm(n_mets, mean = 0, sd = 3)
        if (!dead_ends && (all(flux_vals < 0) || all(flux_vals > 0))) {
            idx <- sample(seq_along(flux_vals), 1)
            flux_vals[idx] <- flux_vals[idx] * -1
        }
        tibble::tibble(
            species = rep(sp, n_mets),
            metabolites = species_mets,
            fluxes = flux_vals
        )
    })

    community <- do.call(rbind, per_species)

    if (!cm) {
        return(community)
    }

    ConsortiumMetabolism(
        data = tibble::tibble(
            species = community$species,
            met = community$metabolites,
            flux = community$fluxes
        ),
        name = name
    )
}
