#' @include AllClasses.R AllGenerics.R
NULL

## ---- CMA accessor methods ------------------------------------------------

#' @describeIn scores Scores from a
#'   [ConsortiumMetabolismAlignment]
#' @export
setMethod(
    "scores",
    "ConsortiumMetabolismAlignment",
    function(object) {
        s <- object@Scores
        if (
            !is.na(object@Metric) &&
                object@Metric == "MAAS" &&
                !is.na(object@PrimaryScore)
        ) {
            s$MAAS <- object@PrimaryScore
        }
        if (!is.na(object@Pvalue)) {
            s$pvalue <- object@Pvalue
        }
        s
    }
)


#' @describeIn similarityMatrix Similarity matrix from a
#'   [ConsortiumMetabolismAlignment]. For multiple alignments
#'   this is the symmetric n x n pairwise matrix; for database
#'   search alignments it is a 1 x n row vector of the query's
#'   scores against each database member.
#' @export
setMethod(
    "similarityMatrix",
    "ConsortiumMetabolismAlignment",
    function(object) {
        if (
            is.na(object@Type) ||
                !object@Type %in% c("multiple", "search")
        ) {
            cli::cli_abort(
                "{.fun similarityMatrix} is only \\
                available for multiple or search \\
                alignments."
            )
        }
        object@SimilarityMatrix
    }
)

#' @describeIn prevalence Pathway prevalence from a
#'   multiple [ConsortiumMetabolismAlignment]
#' @export
setMethod(
    "prevalence",
    "ConsortiumMetabolismAlignment",
    function(object) {
        if (is.na(object@Type) || object@Type != "multiple") {
            cli::cli_abort(
                "{.fun prevalence} is only \\
                available for multiple alignments."
            )
        }
        object@Prevalence
    }
)

## ---- Additional CMA accessors ----------------------------------------------

#' @describeIn pathways Get Pathways From a
#'   \code{ConsortiumMetabolismAlignment} Object
#' @export
setMethod(
    "pathways",
    "ConsortiumMetabolismAlignment",
    function(
        object,
        type = c(
            "all",
            "shared",
            "unique",
            "consensus"
        ),
        verbose = FALSE
    ) {
        type <- match.arg(type)
        alnType <- object@Type

        if (type == "shared") {
            if (!alnType %in% c("pairwise", "search")) {
                cli::cli_abort(
                    "{.arg type} = {.val shared} is \\
                    only available for pairwise or \\
                    search alignments, not \\
                    {.val {alnType}}."
                )
            }
            pw <- object@SharedPathways
            if (verbose) {
                return(pw)
            }
            pw[, c("consumed", "produced")]
        } else if (type == "unique") {
            if (!alnType %in% c("pairwise", "search")) {
                cli::cli_abort(
                    "{.arg type} = {.val unique} is \\
                    only available for pairwise or \\
                    search alignments, not \\
                    {.val {alnType}}."
                )
            }
            if (verbose) {
                list(
                    query = object@UniqueQuery,
                    reference = object@UniqueReference
                )
            } else {
                list(
                    query = object@UniqueQuery[,
                        c("consumed", "produced")
                    ],
                    reference = object@UniqueReference[,
                        c(
                            "consumed",
                            "produced"
                        )
                    ]
                )
            }
        } else if (type == "consensus") {
            if (alnType != "multiple") {
                cli::cli_abort(
                    "{.arg type} = {.val consensus} \\
                    is only available for multiple \\
                    alignments, not \\
                    {.val {alnType}}."
                )
            }
            pw <- object@ConsensusPathways
            if (verbose) {
                return(pw)
            }
            pw[,
                c(
                    "consumed",
                    "produced",
                    "nConsortia",
                    "proportion"
                )
            ]
        } else {
            # type == "all"
            if (verbose) {
                return(object@Pathways)
            }
            object@Pathways[,
                c("consumed", "produced")
            ]
        }
    }
)

#' @rdname metabolites
setMethod(
    "metabolites",
    "ConsortiumMetabolismAlignment",
    function(object) {
        object@Metabolites$met
    }
)

#' @describeIn species Return species from a
#'   \code{ConsortiumMetabolismAlignment}
#' @param object A \code{ConsortiumMetabolismAlignment} object.
setMethod(
    "species",
    "ConsortiumMetabolismAlignment",
    function(object) {
        if (object@Type == "pairwise") {
            sp <- unique(c(
                unlist(lapply(
                    object@SharedPathways$querySpecies,
                    `[[`,
                    "species"
                )),
                unlist(lapply(
                    object@SharedPathways$referenceSpecies,
                    `[[`,
                    "species"
                ))
            ))
            return(sort(sp[!is.na(sp)]))
        }
        cli::cli_abort(
            "{.fun species} is not available for \\
            multiple alignments. Use the original \\
            {.cls ConsortiumMetabolismSet} instead."
        )
    }
)

#' @describeIn consortia Not applicable for alignments
#' @param object A \code{ConsortiumMetabolismAlignment} object.
#' @export
setMethod(
    "consortia",
    "ConsortiumMetabolismAlignment",
    function(object) {
        cli::cli_abort(
            "{.fun consortia} is not applicable for \\
            {.cls ConsortiumMetabolismAlignment} \\
            objects. Use {.fun scores}, \\
            {.fun pathways}, or \\
            {.fun similarityMatrix} instead."
        )
    }
)
