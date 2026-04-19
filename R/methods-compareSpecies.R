#' @include AllClasses.R AllGenerics.R
NULL

## ---- compareSpecies helpers ------------------------------------------------

#' Extract pathway set for one species from a CM
#' @noRd
.speciesPathwaySet <- function(cm, sp) {
    if (!sp %in% unique(cm@InputData$species)) {
        cli::cli_abort(
            "Species {.val {sp}} not found in \\
            {.val {cm@Name}}."
        )
    }
    cm@Pathways |>
        tidyr::unnest("data") |>
        dplyr::filter(.data$species == sp) |>
        dplyr::select("consumed", "produced") |>
        dplyr::distinct()
}

#' Compute similarity metrics between two pathway sets
#' @noRd
.comparePathwaySets <- function(set1, set2) {
    pairs1 <- paste0(
        set1$consumed,
        "-",
        set1$produced
    )
    pairs2 <- paste0(
        set2$consumed,
        "-",
        set2$produced
    )

    n1 <- length(pairs1)
    n2 <- length(pairs2)
    n_shared <- length(intersect(pairs1, pairs2))
    min_n <- min(n1, n2)
    union_n <- n1 + n2 - n_shared

    list(
        fos = if (min_n == 0L) 0 else n_shared / min_n,
        jaccard = if (union_n == 0L) 0 else n_shared / union_n,
        n_shared = n_shared,
        n_unique_sp1 = n1 - n_shared,
        n_unique_sp2 = n2 - n_shared
    )
}

## ---- compareSpecies methods ------------------------------------------------

#' @describeIn compareSpecies Compare two species within
#'   the same [ConsortiumMetabolism]
#' @param x A [ConsortiumMetabolism] object.
#' @param y Character scalar; name of species 1.
#' @param sp2 Character scalar; name of species 2.
#' @export
setMethod(
    "compareSpecies",
    signature(
        x = "ConsortiumMetabolism",
        y = "character"
    ),
    function(x, y, sp2, ...) {
        set1 <- .speciesPathwaySet(x, y)
        set2 <- .speciesPathwaySet(x, sp2)
        .comparePathwaySets(set1, set2)
    }
)

#' @describeIn compareSpecies Compare one species from
#'   each of two [ConsortiumMetabolism] objects
#' @param x A [ConsortiumMetabolism] object (first).
#' @param y A [ConsortiumMetabolism] object (second).
#' @param sp1 Character scalar; species from \code{x}.
#' @param sp2 Character scalar; species from \code{y}.
#' @export
setMethod(
    "compareSpecies",
    signature(
        x = "ConsortiumMetabolism",
        y = "ConsortiumMetabolism"
    ),
    function(x, y, sp1, sp2, ...) {
        set1 <- .speciesPathwaySet(x, sp1)
        set2 <- .speciesPathwaySet(y, sp2)
        .comparePathwaySets(set1, set2)
    }
)
