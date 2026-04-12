#' Normalize Bigg exchange reaction IDs to bare metabolite names
#' @param ids Character vector of reaction/metabolite identifiers.
#' @return Character vector of cleaned metabolite names.
#' @keywords internal
.normalizeBiggIds <- function(ids) {
    ids |>
        # Misosoup
        sub("^R_EX_", "", x = _) |>
        # MICOM, COMETS, BacArena
        sub("^EX_", "", x = _) |>
        # SMETANA SBML Prefix
        sub("^M_", "", x = _) |>
        # BacArena Compartment
        sub("\\(e\\)$", "", x = _) |>
        # COMETS bracket compartment
        sub("\\[e\\]$", "", x = _) |>
        # Standard BiGG compartment
        sub("_e$", "", x = _)
}
