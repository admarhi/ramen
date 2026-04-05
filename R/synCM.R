#' @noRd
#' @keywords internal
.rNames <- function(n = 5000L) {
    a <- do.call(
        paste0,
        replicate(3, sample(LETTERS, n, TRUE), FALSE)
    )
    paste0(
        a,
        sprintf("%03d", sample(9999L, n, TRUE)),
        sample(LETTERS, n, TRUE)
    )
}


#' Zipf-distributed metabolite sampling weights
#'
#' Produces a randomly-permuted Zipf weight vector so that a few
#' metabolites act as hubs while the rest are rarer.
#'
#' @param max_met Integer.
#' @param shape Numeric exponent (default 1.5).
#' @return Named numeric vector summing to 1.
#' @noRd
#' @keywords internal
.metWeights <- function(max_met, shape = 1.5) {
    w <- seq_len(max_met)^(-shape)
    w <- sample(w)
    w <- w / sum(w)
    names(w) <- paste0("met", seq_len(max_met))
    w
}


#' Log-normal species degree sequence
#'
#' Draws how many metabolites each species interacts with from a
#' discretised log-normal distribution, clamped to [2, max_met].
#'
#' @param n_species Integer.
#' @param max_met Integer.
#' @param meanlog Numeric (default 1.0).
#' @param sdlog Numeric (default 0.8).
#' @return Integer vector of length \code{n_species}.
#' @noRd
#' @keywords internal
.speciesDegrees <- function(n_species, max_met,
                            meanlog = 1.0, sdlog = 0.8) {
    raw <- stats::rlnorm(n_species, meanlog, sdlog)
    deg <- as.integer(round(raw))
    deg <- pmax(deg, 2L)
    deg <- pmin(deg, max_met)
    deg
}


#' Create cyclic cross-feeding backbone
#'
#' Arranges species in a ring and assigns one shared metabolite
#' per consecutive pair so every species cross-feeds with at least
#' one neighbour.
#'
#' @param species_names Character vector.
#' @param met_vec Character vector of metabolite IDs.
#' @param met_weights Named numeric vector (sampling weights).
#' @return A \code{tibble} with columns \code{species},
#'   \code{metabolites}, \code{fluxes}.
#' @noRd
#' @keywords internal
.wireCrossFeeding <- function(species_names, met_vec,
                              met_weights) {
    n <- length(species_names)
    if (n < 2L) {
        return(tibble::tibble(
            species = character(0L),
            metabolites = character(0L),
            fluxes = numeric(0L)
        ))
    }

    n_chain <- min(n, max(2L, as.integer(
        floor(length(met_vec) * 0.3)
    )))
    chain_mets <- sample(
        met_vec,
        size = n_chain,
        prob = met_weights[met_vec]
    )

    sp_out <- character(0L)
    met_out <- character(0L)
    flux_out <- numeric(0L)

    for (i in seq_len(n)) {
        j <- if (i < n) i + 1L else 1L
        m <- chain_mets[((i - 1L) %% n_chain) + 1L]
        mag <- stats::rlnorm(1L, meanlog = 0.5, sdlog = 0.5)

        # upstream produces
        sp_out <- c(sp_out, species_names[i])
        met_out <- c(met_out, m)
        flux_out <- c(flux_out, mag)

        # downstream consumes
        sp_out <- c(sp_out, species_names[j])
        met_out <- c(met_out, m)
        flux_out <- c(flux_out, -mag)
    }

    tibble::tibble(
        species = sp_out,
        metabolites = met_out,
        fluxes = flux_out
    )
}


#' Fill remaining edges up to target degree
#'
#' For each species, adds random metabolite interactions until its
#' degree target is met.  Consumption/production ratio is drawn
#' per species from a Beta(2,2) distribution.
#'
#' @param species_names Character vector.
#' @param met_vec Character vector of metabolite IDs.
#' @param met_weights Named numeric vector.
#' @param degrees Integer vector (per-species target degree).
#' @param chain_edges Tibble returned by \code{.wireCrossFeeding}.
#' @return A \code{tibble} with columns \code{species},
#'   \code{metabolites}, \code{fluxes}.
#' @noRd
#' @keywords internal
.fillRandomEdges <- function(species_names, met_vec,
                             met_weights, degrees,
                             chain_edges) {
    sp_out <- character(0L)
    met_out <- character(0L)
    flux_out <- numeric(0L)

    for (idx in seq_along(species_names)) {
        sp <- species_names[idx]
        existing <- chain_edges$metabolites[
            chain_edges$species == sp
        ]
        n_have <- length(existing)
        n_need <- max(0L, degrees[idx] - n_have)

        if (n_need == 0L) next

        available <- setdiff(met_vec, existing)
        if (length(available) == 0L) next

        n_add <- min(n_need, length(available))
        w <- met_weights[available]
        w <- w / sum(w)
        new_mets <- sample(
            available,
            size = n_add,
            prob = w
        )

        mags <- stats::rlnorm(
            n_add, meanlog = 0.5, sdlog = 0.8
        )
        r <- stats::rbeta(1L, 2, 2)
        n_cons <- max(
            0L,
            min(n_add, as.integer(round(n_add * r)))
        )
        signs <- c(
            rep(-1, n_cons),
            rep(1, n_add - n_cons)
        )
        signs <- sample(signs)

        sp_out <- c(sp_out, rep(sp, n_add))
        met_out <- c(met_out, new_mets)
        flux_out <- c(flux_out, signs * mags)
    }

    tibble::tibble(
        species = sp_out,
        metabolites = met_out,
        fluxes = flux_out
    )
}


#' Approximate mass balance for exchanged metabolites
#'
#' For metabolites that appear as both produced and consumed,
#' rescales the smaller side so total production and consumption
#' are within a slack tolerance, then adds small noise.
#'
#' @param community Tibble with columns \code{species},
#'   \code{metabolites}, \code{fluxes}.
#' @param slack Numeric tolerance (default 0.2).
#' @return The modified tibble.
#' @noRd
#' @keywords internal
.balanceFluxes <- function(community, slack = 0.2) {
    if (nrow(community) == 0L) return(community)

    mets <- unique(community$metabolites)

    for (m in mets) {
        rows <- which(community$metabolites == m)
        fluxes <- community$fluxes[rows]
        prod_idx <- rows[fluxes > 0]
        cons_idx <- rows[fluxes < 0]

        if (length(prod_idx) == 0L || length(cons_idx) == 0L) {
            next
        }

        total_p <- sum(community$fluxes[prod_idx])
        total_c <- sum(abs(community$fluxes[cons_idx]))

        ratio <- total_p / total_c

        if (ratio < (1 - slack) || ratio > (1 + slack)) {
            # rescale the smaller side toward balance
            if (ratio < (1 - slack)) {
                community$fluxes[prod_idx] <-
                    community$fluxes[prod_idx] *
                    (total_c / total_p)
            } else {
                community$fluxes[cons_idx] <-
                    community$fluxes[cons_idx] *
                    (total_p / total_c)
            }
            # add small noise so values are not exact
            balanced <- c(prod_idx, cons_idx)
            noise <- stats::runif(
                length(balanced),
                1 - slack / 2,
                1 + slack / 2
            )
            community$fluxes[balanced] <-
                community$fluxes[balanced] * noise
        }
    }

    community
}


#' Generate Synthetic Consortium Metabolism
#'
#' Creates a synthetic metabolic network with biologically
#' realistic structure: hub metabolites (Zipf-distributed
#' degree), log-normal species degrees, cyclic cross-feeding
#' backbone for connectivity, approximate mass balance, and
#' no dead-end species by default.
#'
#' @param name Character string. Name for the consortium.
#' @param n_species Integer. Number of species in the
#'     consortium.
#' @param max_met Integer. Size of the metabolite pool to
#'     sample from.
#' @param seed Integer or \code{FALSE}. Random seed for
#'     reproducibility. If \code{FALSE} (default), a random
#'     seed is chosen.
#' @param dead_ends Logical. If \code{FALSE} (default), ensures
#'     each species has both consumed and produced metabolites
#'     by flipping one flux value when all fluxes share the
#'     same sign.
#' @param cm Logical. If \code{TRUE} (default), returns a
#'     \code{\linkS4class{ConsortiumMetabolism}} object. If
#'     \code{FALSE}, returns the raw edge list as a tibble.
#'
#' @return A \code{\linkS4class{ConsortiumMetabolism}} object
#'     when \code{cm = TRUE}, or a \code{\link[tibble]{tibble}}
#'     with columns \code{species}, \code{metabolites}, and
#'     \code{fluxes} when \code{cm = FALSE}.
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
    seed = FALSE,
    dead_ends = FALSE,
    cm = TRUE
) {
    ## ---- input validation ----
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

    ## ---- seed handling ----
    if (!isFALSE(seed)) {
        set.seed(seed)
    } else {
        seed <- sample(seq_len(1000L), 1L)
        set.seed(seed)
    }

    ## ---- generate names ----
    species_names <- .rNames(n_species)
    met_vec <- paste0("met", seq_len(max_met))

    ## ---- stage 1: metabolite hub weights ----
    met_weights <- .metWeights(max_met)

    ## ---- stage 2: species degree sequence ----
    degrees <- .speciesDegrees(n_species, max_met)

    ## ---- stage 3: cross-feeding backbone ----
    chain_edges <- .wireCrossFeeding(
        species_names, met_vec, met_weights
    )

    ## ---- stage 4: fill random edges ----
    random_edges <- .fillRandomEdges(
        species_names, met_vec, met_weights,
        degrees, chain_edges
    )

    ## ---- stage 5: approximate mass balance ----
    community <- .balanceFluxes(
        rbind(chain_edges, random_edges)
    )

    ## ---- dead-end fix ----
    if (!dead_ends) {
        for (sp in species_names) {
            rows <- which(community$species == sp)
            fluxes <- community$fluxes[rows]
            if (all(fluxes > 0) || all(fluxes < 0)) {
                idx <- sample(rows, 1L)
                community$fluxes[idx] <-
                    community$fluxes[idx] * -1
            }
        }
    }

    ## ---- return ----
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
