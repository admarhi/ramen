#' @include AllClasses.R AllGenerics.R
NULL

## ---- Internal helpers ------------------------------------------------------

#' Compute functional groups from a species-pathway table
#'
#' Shared workhorse for both
#' \code{functionalGroups,ConsortiumMetabolism-method} and
#' \code{functionalGroups,ConsortiumMetabolismSet-method}.
#'
#' @param sp_pw A data.frame with columns \code{species}
#'   and \code{pathway} (one row per species-pathway
#'   participation). Duplicates are tolerated and removed
#'   internally.
#' @param linkage Character scalar; agglomeration method
#'   passed to \code{\link[stats]{hclust}}.
#'
#' @return A named list with elements \code{dendrogram},
#'   \code{similarity_matrix}, \code{incidence_matrix},
#'   and \code{reactions_per_species}. When fewer than
#'   two species are present, \code{dendrogram} is
#'   \code{NULL}.
#'
#' @noRd
#' @keywords internal
.functionalGroupsFromTable <- function(sp_pw, linkage) {
    sp_pw <- unique(sp_pw[, c("species", "pathway")])

    sp_levels <- sort(unique(sp_pw$species))
    pw_levels <- unique(sp_pw$pathway)
    n_sp <- length(sp_levels)
    n_pw <- length(pw_levels)

    sp_idx <- match(sp_pw$species, sp_levels)
    pw_idx <- match(sp_pw$pathway, pw_levels)

    incidence <- Matrix::sparseMatrix(
        i = sp_idx,
        j = pw_idx,
        x = 1,
        dims = c(n_sp, n_pw),
        dimnames = list(sp_levels, pw_levels)
    )

    if (n_sp < 2L) {
        ## Single-species (or zero) edge case: similarity
        ## is the trivial 1x1 (or 0x0) matrix and there
        ## is no clustering to do.
        similarity_matrix <- if (n_sp == 1L) {
            matrix(
                1,
                nrow = 1L,
                ncol = 1L,
                dimnames = list(sp_levels, sp_levels)
            )
        } else {
            matrix(
                numeric(0L),
                nrow = 0L,
                ncol = 0L
            )
        }
        return(list(
            dendrogram = NULL,
            similarity_matrix = similarity_matrix,
            incidence_matrix = incidence,
            reactions_per_species = sp_pw
        ))
    }

    ## Jaccard via crossprod:
    ##   intersection = A %*% t(A)
    ##   union = |A_i| + |A_j| - intersection
    intersection <- Matrix::tcrossprod(incidence)
    pathway_counts <- Matrix::rowSums(incidence)
    union_mat <- outer(
        pathway_counts,
        pathway_counts,
        "+"
    ) -
        as.matrix(intersection)
    similarity_matrix <-
        as.matrix(intersection) / union_mat
    ## Fix 0/0 = NaN on diagonal or empty species
    similarity_matrix[is.nan(similarity_matrix)] <- 0
    diag(similarity_matrix) <- 1

    dend <- similarity_matrix |>
        stats::dist() |>
        stats::hclust(method = linkage) |>
        stats::as.dendrogram()

    list(
        dendrogram = dend,
        similarity_matrix = similarity_matrix,
        incidence_matrix = incidence,
        reactions_per_species = sp_pw
    )
}

#' Resolve common functionalGroups() arguments
#'
#' Centralises the deprecation warning for plotting
#' arguments and the \code{linkage} default so both
#' methods stay in lockstep.
#'
#' @param dots A list (typically captured via
#'   \code{list(...)} inside an S4 method).
#'
#' @return A character scalar: the validated linkage
#'   method.
#'
#' @noRd
#' @keywords internal
.resolveFunctionalGroupsArgs <- function(dots) {
    linkage <- if (is.null(dots$linkage)) {
        "complete"
    } else {
        dots$linkage
    }

    viz_args <- intersect(
        names(dots),
        c("k", "label_size", "label_colours")
    )
    if (length(viz_args) > 0L) {
        cli::cli_warn(
            c(
                "!" = paste0(
                    "Argument{?s} {.arg {viz_args}}",
                    " moved to",
                    " {.fun plotFunctionalGroups}."
                ),
                "i" = paste0(
                    "Use ",
                    "{.code plotFunctionalGroups",
                    "(fg, ...)} to visualize."
                )
            )
        )
    }

    match.arg(
        linkage,
        c("complete", "average", "single", "ward.D2")
    )
}

## ---- ConsortiumMetabolism method -------------------------------------------

#' @describeIn functionalGroups Functional groups within a
#'   single \code{ConsortiumMetabolism}. Builds a species
#'   x pathway incidence from the consortium's
#'   \code{Pathways} slot (one pathway per unique
#'   \code{(consumed, produced)} pair) and clusters
#'   species by Jaccard similarity over their pathway
#'   sets. If the consortium contains fewer than two
#'   species, a \code{cli::cli_warn} is emitted and the
#'   returned list has \code{dendrogram = NULL}.
#' @export
setMethod(
    "functionalGroups",
    "ConsortiumMetabolism",
    function(object, ...) {
        dots <- list(...)
        linkage <- .resolveFunctionalGroupsArgs(dots)

        ## Unnest per-species detail from the CM Pathways
        ## table, mirroring the species-pathway pairing
        ## the CMS variant uses (one pair per species per
        ## (consumed, produced) row).
        rxns_per_species <- object@Pathways |>
            dplyr::select(
                "consumed",
                "produced",
                "data"
            ) |>
            tidyr::unnest("data") |>
            dplyr::mutate(
                pathway = paste0(
                    .data$consumed,
                    "-",
                    .data$produced
                )
            ) |>
            dplyr::select("species", "pathway") |>
            unique()

        n_sp <- length(unique(rxns_per_species$species))
        if (n_sp < 2L) {
            cli::cli_warn(
                c(
                    "!" = paste0(
                        "Consortium has {n_sp} species; ",
                        "need at least 2 to compute ",
                        "functional groups."
                    ),
                    "i" = paste0(
                        "Returning incidence matrix only;",
                        " {.field dendrogram} will be ",
                        "{.code NULL}."
                    )
                )
            )
        }

        invisible(
            .functionalGroupsFromTable(
                rxns_per_species,
                linkage
            )
        )
    }
)

## ---- ConsortiumMetabolismSet method ----------------------------------------

#' @describeIn functionalGroups Functional groups across a
#'   \code{ConsortiumMetabolismSet}. Pools species-pathway
#'   pairs from every consortium's \code{Pathways} table
#'   and clusters species by Jaccard similarity.
#' @export
setMethod(
    "functionalGroups",
    "ConsortiumMetabolismSet",
    function(object, ...) {
        dots <- list(...)
        linkage <- .resolveFunctionalGroupsArgs(dots)

        rxns_per_species <- object@Pathways |>
            dplyr::select(
                "consumed",
                "produced",
                "species"
            ) |>
            unique() |>
            dplyr::mutate(
                pathway = paste0(
                    .data$consumed,
                    "-",
                    .data$produced
                )
            ) |>
            dplyr::select("species", "pathway")

        invisible(
            .functionalGroupsFromTable(
                rxns_per_species,
                linkage
            )
        )
    }
)
