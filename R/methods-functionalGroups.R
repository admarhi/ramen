#' @include AllClasses.R AllGenerics.R
NULL

#' @rdname functionalGroups
#' @aliases functionalGroups,ConsortiumMetabolismSet-method
#' @export
setMethod(
    "functionalGroups",
    "ConsortiumMetabolismSet",
    function(object, ...) {
        dots <- list(...)
        linkage <- if (is.null(dots$linkage)) {
            "complete"
        } else {
            dots$linkage
        }

        # Warn about deprecated visualization arguments
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

        linkage <- match.arg(
            linkage,
            c("complete", "average", "single",
              "ward.D2")
        )
        # Extract pathways from the CMS
        tb <- object@Pathways

        # Build unique species-pathway pairs
        rxns_per_species <- tb |>
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

        # Build binary species-by-pathway incidence
        # matrix
        sp_levels <- sort(
            unique(rxns_per_species$species)
        )
        pw_levels <- unique(
            rxns_per_species$pathway
        )
        sp_idx <- match(
            rxns_per_species$species,
            sp_levels
        )
        pw_idx <- match(
            rxns_per_species$pathway,
            pw_levels
        )
        incidence <- Matrix::sparseMatrix(
            i = sp_idx,
            j = pw_idx,
            x = 1,
            dims = c(
                length(sp_levels),
                length(pw_levels)
            ),
            dimnames = list(sp_levels, pw_levels)
        )

        # Jaccard via crossprod:
        #   intersection = A %*% t(A)
        #   union = |A_i| + |A_j| - intersection
        intersection <- Matrix::tcrossprod(incidence)
        pathway_counts <- Matrix::rowSums(incidence)
        union_mat <- outer(
            pathway_counts,
            pathway_counts,
            "+"
        ) - as.matrix(intersection)
        similarity_matrix <-
            as.matrix(intersection) / union_mat
        # Fix 0/0 = NaN on diagonal or empty species
        similarity_matrix[
            is.nan(similarity_matrix)
        ] <- 0
        diag(similarity_matrix) <- 1

        # Hierarchical clustering on distance matrix
        dend <- similarity_matrix |>
            stats::dist() |>
            stats::hclust(method = linkage) |>
            stats::as.dendrogram()

        # Return list with dendrogram and data
        result <- list(
            dendrogram = dend,
            similarity_matrix = similarity_matrix,
            incidence_matrix = incidence,
            reactions_per_species = rxns_per_species
        )

        invisible(result)
    }
)
