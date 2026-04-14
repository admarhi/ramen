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

#' Decode COBRApy/SBML `__NN__` escape sequences
#'
#' COBRApy encodes non-alphanumeric characters in reaction and metabolite
#' IDs as `__<decimal-ASCII-code>__` sequences when writing SBML. For
#' example, `(` becomes `__40__` and `)` becomes `__41__`. This helper
#' reverses that encoding so downstream parsing sees the literal
#' characters.
#'
#' @param ids Character vector of identifiers possibly containing
#'   `__NN__` escape sequences.
#' @return Character vector with escape sequences replaced by their
#'   literal characters. IDs with no escape sequences pass through
#'   unchanged.
#' @keywords internal
.decodeBiggEscapes <- function(ids) {
    matches <- regmatches(
        ids,
        gregexpr("__\\d+__", ids, perl = TRUE)
    )
    codes <- unique(unlist(matches))
    for (code in codes) {
        num <- as.integer(sub("__(\\d+)__", "\\1", code))
        ids <- gsub(code, intToUtf8(num), ids, fixed = TRUE)
    }
    ids
}
