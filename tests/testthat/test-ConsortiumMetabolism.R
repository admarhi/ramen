test_that("ConsortiumMetabolism constructor works with basic data", {
    # Create simple test data
    test_data <- tibble::tibble(
        species = c("s1", "s1", "s2", "s2"),
        metabolite = c("m1", "m2", "m1", "m3"),
        flux = c(-1, 1, -1, 1)
    )

    # Test constructor
    cm <- ConsortiumMetabolism(test_data, name = "test_cm")

    # Basic assertions
    expect_s4_class(cm, "ConsortiumMetabolism")
    expect_equal(cm@Name, "test_cm")
    expect_true(length(cm@Metabolites) > 0)
    expect_s3_class(cm@Pathways, "data.frame")
})

test_that("ConsortiumMetabolism handles weighted vs unweighted networks", {
    # Unweighted data (all fluxes = 1)
    unweighted_data <- tibble::tibble(
        species = c("s1", "s1", "s2", "s2"),
        metabolite = c("m1", "m2", "m1", "m3"),
        flux = c(-1, 1, -1, 1)
    )

    cm_unweighted <- ConsortiumMetabolism(unweighted_data, name = "unweighted")
    expect_false(cm_unweighted@Weighted)

    # Weighted data (varying fluxes)
    weighted_data <- tibble::tibble(
        species = c("s1", "s1", "s2", "s2"),
        metabolite = c("m1", "m2", "m1", "m3"),
        flux = c(-2, 3, -1.5, 2.5)
    )

    cm_weighted <- ConsortiumMetabolism(weighted_data, name = "weighted")
    expect_true(cm_weighted@Weighted)
})

test_that("species returns correct species", {
    test_data <- tibble::tibble(
        species = c("s1", "s1", "s2", "s2", "s3", "s3"),
        metabolite = c("m1", "m2", "m1", "m3", "m2", "m4"),
        flux = c(-1, 1, -1, 1, -1, 1)
    )

    cm <- ConsortiumMetabolism(test_data, name = "test")
    sp <- species(cm)

    expect_type(sp, "character")
    expect_true(length(sp) >= 2)
    expect_true(all(c("s1", "s2", "s3") %in% sp))
})

test_that("metabolites returns metabolites", {
    test_data <- tibble::tibble(
        species = c("s1", "s1", "s2", "s2"),
        metabolite = c("m1", "m2", "m1", "m3"),
        flux = c(-1, 1, -1, 1)
    )

    cm <- ConsortiumMetabolism(test_data, name = "test")
    mets <- metabolites(cm)

    expect_type(mets, "character")
    expect_true(length(mets) > 0)
})

test_that("metabolites filters by species and direction", {
    test_data <- tibble::tibble(
        species = c("s1", "s1", "s2", "s2"),
        metabolite = c("m1", "m2", "m1", "m3"),
        flux = c(-1, 1, -1, 1)
    )
    cm <- ConsortiumMetabolism(test_data, name = "test")

    ## All metabolites involved with s1
    s1_all <- metabolites(cm, species = "s1")
    expect_true(all(c("m1", "m2") %in% s1_all))

    ## s1 consumes m1, produces m2
    expect_equal(
        metabolites(cm, species = "s1", direction = "consumed"),
        "m1"
    )
    expect_equal(
        metabolites(cm, species = "s1", direction = "produced"),
        "m2"
    )

    ## Direction without species applies globally
    consumed_all <- metabolites(cm, direction = "consumed")
    produced_all <- metabolites(cm, direction = "produced")
    expect_true("m1" %in% consumed_all) ## consumed by s1 and s2
    expect_true(all(c("m2", "m3") %in% produced_all))

    ## Errors on unknown species
    expect_error(
        metabolites(cm, species = "no_such_species"),
        "not found"
    )
    ## Errors on non-character species
    expect_error(
        metabolites(cm, species = 1L),
        "character scalar"
    )
})

test_that("pathways returns concise output by default", {
    test_data <- tibble::tibble(
        species = c("s1", "s1", "s2", "s2"),
        metabolite = c("m1", "m2", "m1", "m3"),
        flux = c(-1, 1, -1, 1)
    )

    cm <- ConsortiumMetabolism(test_data, name = "test")
    pw <- pathways(cm)

    expect_s3_class(pw, "data.frame")
    expect_true(nrow(pw) > 0)
    expect_named(pw, c("consumed", "produced", "n_species"))
})

test_that("pathways verbose returns full detail", {
    test_data <- tibble::tibble(
        species = c("s1", "s1", "s2", "s2"),
        metabolite = c("m1", "m2", "m1", "m3"),
        flux = c(-1, 1, -1, 1)
    )

    cm <- ConsortiumMetabolism(test_data, name = "test")
    pw <- pathways(cm, verbose = TRUE)

    expect_s3_class(pw, "data.frame")
    expect_true(ncol(pw) > 3)
    expect_true("c_sum" %in% names(pw))
    expect_true("c_eff" %in% names(pw))
    expect_true("data" %in% names(pw))
})

test_that("CM constructor errors on non-data.frame input", {
    expect_error(
        ConsortiumMetabolism("not a df", name = "x"),
        "must be a data.frame"
    )
    expect_error(
        ConsortiumMetabolism(list(a = 1), name = "x"),
        "must be a data.frame"
    )
})

## ---- input validation: NA / type / shape ------------------------------------

test_that("CM constructor errors on NA in species column", {
    bad <- tibble::tibble(
        species = c("s1", NA, "s2"),
        metabolite = c("m1", "m2", "m1"),
        flux = c(-1, 1, -1)
    )
    expect_error(
        ConsortiumMetabolism(bad, name = "x"),
        "NA"
    )
})

test_that("CM constructor errors on NA in metabolite column", {
    bad <- tibble::tibble(
        species = c("s1", "s1", "s2"),
        metabolite = c("m1", NA, "m1"),
        flux = c(-1, 1, -1)
    )
    expect_error(
        ConsortiumMetabolism(bad, name = "x"),
        "NA"
    )
})

test_that("CM constructor errors on NA in flux column", {
    bad <- tibble::tibble(
        species = c("s1", "s1", "s2"),
        metabolite = c("m1", "m2", "m1"),
        flux = c(-1, NA, -1)
    )
    err <- expect_error(
        ConsortiumMetabolism(bad, name = "x")
    )
    ## Should be the friendly NA-row message, not the validity-leak
    ## message about 'Weighted'.
    expect_match(conditionMessage(err), "NA")
    expect_false(grepl("Weighted", conditionMessage(err)))
})

test_that("Weighted validity message is user-facing, not slot-named", {
    ## Build a valid CM, then break the @Weighted slot to trigger
    ## validity. The message must guide users to the constructor
    ## rather than leak slot internals.
    cm <- synCM("a", n_species = 2, max_met = 3, seed = 1)
    cm@Weighted <- logical(0L)
    err <- expect_error(validObject(cm))
    msg <- conditionMessage(err)
    expect_match(msg, "ConsortiumMetabolism\\(")
    expect_match(msg, "TRUE or FALSE")
})

test_that("CM constructor errors on empty data.frame", {
    empty <- tibble::tibble(
        species = character(),
        metabolite = character(),
        flux = numeric()
    )
    expect_error(
        ConsortiumMetabolism(empty, name = "x"),
        "0 rows"
    )
})

test_that("CM constructor errors on numeric species column", {
    bad <- data.frame(
        species = c(1, 2),
        metabolite = c("a", "b"),
        flux = c(-1, 1)
    )
    expect_error(
        ConsortiumMetabolism(bad, name = "x"),
        "character"
    )
})

test_that("CM constructor errors on numeric metabolite column", {
    bad <- data.frame(
        species = c("s1", "s2"),
        metabolite = c(1, 2),
        flux = c(-1, 1)
    )
    expect_error(
        ConsortiumMetabolism(bad, name = "x"),
        "character"
    )
})

test_that("CM constructor rejects misspelled named arguments", {
    test_data <- tibble::tibble(
        species = c("s1", "s1", "s2", "s2"),
        metabolite = c("m1", "m2", "m1", "m3"),
        flux = c(-1, 1, -1, 1)
    )
    ## Capital N — should NOT be silently swallowed; with `...` removed
    ## from the constructor, R's own "unused argument" error fires.
    expect_error(
        ConsortiumMetabolism(test_data, Name = "test_cm"),
        "unused argument"
    )
})

test_that("CM constructor errors when all flux values are zero", {
    bad <- tibble::tibble(
        species = c("s1", "s1", "s2"),
        metabolite = c("m1", "m2", "m1"),
        flux = c(0, 0, 0)
    )
    expect_error(
        ConsortiumMetabolism(bad, name = "x"),
        "zero"
    )
})

test_that("CM constructor still works with valid input (positive control)", {
    good <- tibble::tibble(
        species = c("s1", "s1", "s2", "s2"),
        metabolite = c("m1", "m2", "m1", "m3"),
        flux = c(-1, 1, -1, 1)
    )
    expect_s4_class(
        ConsortiumMetabolism(good, name = "ok"),
        "ConsortiumMetabolism"
    )
})

## ---- growth parameter -------------------------------------------------------

test_that("CM constructor accepts growth parameter", {
    test_data <- tibble::tibble(
        species = c("s1", "s1", "s2", "s2"),
        metabolite = c("m1", "m2", "m1", "m3"),
        flux = c(-1, 1, -1, 1)
    )
    gr <- c(s1 = 0.5, s2 = 0.3)
    cm <- ConsortiumMetabolism(
        test_data,
        name = "test",
        growth = gr
    )
    expect_equal(growth(cm), gr)
})

test_that("growth() returns NULL when no growth supplied", {
    test_data <- tibble::tibble(
        species = c("s1", "s1", "s2", "s2"),
        metabolite = c("m1", "m2", "m1", "m3"),
        flux = c(-1, 1, -1, 1)
    )
    cm <- ConsortiumMetabolism(test_data, name = "test")
    expect_null(growth(cm))
})

test_that("growth validation: must be numeric", {
    test_data <- tibble::tibble(
        species = c("s1", "s1", "s2", "s2"),
        metabolite = c("m1", "m2", "m1", "m3"),
        flux = c(-1, 1, -1, 1)
    )
    expect_error(
        ConsortiumMetabolism(
            test_data,
            name = "test",
            growth = c(s1 = "high", s2 = "low")
        ),
        "must be a numeric"
    )
})

test_that("growth validation: must be named", {
    test_data <- tibble::tibble(
        species = c("s1", "s1", "s2", "s2"),
        metabolite = c("m1", "m2", "m1", "m3"),
        flux = c(-1, 1, -1, 1)
    )
    expect_error(
        ConsortiumMetabolism(
            test_data,
            name = "test",
            growth = c(0.5, 0.3)
        ),
        "named numeric"
    )
})

test_that("growth validation: names must match species", {
    test_data <- tibble::tibble(
        species = c("s1", "s1", "s2", "s2"),
        metabolite = c("m1", "m2", "m1", "m3"),
        flux = c(-1, 1, -1, 1)
    )
    expect_error(
        ConsortiumMetabolism(
            test_data,
            name = "test",
            growth = c(s1 = 0.5, s99 = 0.3)
        ),
        "not found"
    )
})

test_that("growth validation: no duplicate names", {
    test_data <- tibble::tibble(
        species = c("s1", "s1", "s2", "s2"),
        metabolite = c("m1", "m2", "m1", "m3"),
        flux = c(-1, 1, -1, 1)
    )
    gr <- c(s1 = 0.5, s1 = 0.3)
    expect_error(
        ConsortiumMetabolism(
            test_data,
            name = "test",
            growth = gr
        ),
        "duplicate"
    )
})

test_that("growth validation: non-negative values", {
    test_data <- tibble::tibble(
        species = c("s1", "s1", "s2", "s2"),
        metabolite = c("m1", "m2", "m1", "m3"),
        flux = c(-1, 1, -1, 1)
    )
    expect_error(
        ConsortiumMetabolism(
            test_data,
            name = "test",
            growth = c(s1 = 0.5, s2 = -0.1)
        ),
        "non-negative"
    )
})

test_that("growth accepts subset of species", {
    test_data <- tibble::tibble(
        species = c("s1", "s1", "s2", "s2"),
        metabolite = c("m1", "m2", "m1", "m3"),
        flux = c(-1, 1, -1, 1)
    )
    gr <- c(s1 = 0.5)
    cm <- ConsortiumMetabolism(
        test_data,
        name = "test",
        growth = gr
    )
    expect_equal(growth(cm), gr)
})

test_that("synCM generates synthetic communities", {
    # Test synthetic community generation
    syn <- synCM(name = "synthetic", n_species = 5, max_met = 10)

    expect_s4_class(syn, "ConsortiumMetabolism")
    expect_equal(syn@Name, "synthetic")
    expect_true(nrow(syn@Pathways) > 0)
})

## ---- speciesSummary(CM) -----------------------------------------------------

test_that("speciesSummary(CM) returns tibble with correct columns", {
    test_data <- tibble::tibble(
        species = c("s1", "s1", "s2", "s2"),
        metabolite = c("m1", "m2", "m1", "m3"),
        flux = c(-1, 1, -1, 1)
    )
    cm <- ConsortiumMetabolism(test_data, name = "test")
    ss <- speciesSummary(cm)
    expect_s3_class(ss, "tbl_df")
    expect_named(
        ss,
        c("species", "n_pathways", "n_consumed", "n_produced")
    )
})

test_that("speciesSummary(CM) counts are non-negative integers", {
    cm <- synCM("test", n_species = 4, max_met = 6, seed = 1)
    ss <- speciesSummary(cm)
    expect_true(all(ss$n_pathways >= 1L))
    expect_true(all(ss$n_consumed >= 1L))
    expect_true(all(ss$n_produced >= 1L))
})

test_that("speciesSummary(CM) species match species(cm)", {
    cm <- synCM("test", n_species = 3, max_met = 5, seed = 2)
    ss <- speciesSummary(cm)
    expect_setequal(ss$species, species(cm))
})

## ---- effective-species assays ------------------------------------------

test_that("CM exposes nEffectiveSpecies* assays alongside Effective*", {
    cm <- synCM("test", n_species = 3, max_met = 5, seed = 1)
    expect_true(
        all(
            c(
                "EffectiveConsumption",
                "EffectiveProduction",
                "nEffectiveSpeciesConsumption",
                "nEffectiveSpeciesProduction"
            ) %in%
                SummarizedExperiment::assayNames(cm)
        )
    )
})

test_that("EffectiveConsumption equals Consumption * nEffectiveSpecies (toy)", {
    ## Two species, both consume m1 with equal flux (so perplexity = 2),
    ## both produce m2 with unequal flux (perplexity < 2).
    test_data <- tibble::tibble(
        species = c("s1", "s1", "s2", "s2"),
        metabolite = c("m1", "m2", "m1", "m2"),
        flux = c(-1, 3, -1, 1)
    )
    cm <- ConsortiumMetabolism(test_data, name = "toy")
    a <- SummarizedExperiment::assays(cm)

    cons <- as.matrix(a$Consumption)
    eff_c <- as.matrix(a$EffectiveConsumption)
    nef_c <- as.matrix(a$nEffectiveSpeciesConsumption)

    nz <- cons != 0
    expect_equal(eff_c[nz], round(cons[nz] * nef_c[nz], 2))

    ## Hand-check: at the (m1, m2) cell consumption flux sum is 1+1=2,
    ## even split -> perplexity 2, so effective consumption = 4.
    expect_equal(unname(cons["m1", "m2"]), 2)
    expect_equal(unname(nef_c["m1", "m2"]), 2)
    expect_equal(unname(eff_c["m1", "m2"]), 4)
})

test_that("EffectiveProduction equals Production * nEffectiveSpecies", {
    test_data <- tibble::tibble(
        species = c("s1", "s1", "s2", "s2"),
        metabolite = c("m1", "m2", "m1", "m2"),
        flux = c(-1, 3, -1, 1)
    )
    cm <- ConsortiumMetabolism(test_data, name = "toy")
    a <- SummarizedExperiment::assays(cm)

    prod <- as.matrix(a$Production)
    eff_p <- as.matrix(a$EffectiveProduction)
    nef_p <- as.matrix(a$nEffectiveSpeciesProduction)

    nz <- prod != 0
    expect_equal(eff_p[nz], round(prod[nz] * nef_p[nz], 2))

    ## Production at (m1, m2): fluxes 3 and 1, sum 4, p = (0.75, 0.25),
    ## perplexity 2^(-(0.75 log2 0.75 + 0.25 log2 0.25)) ~= 1.7549.
    expect_equal(unname(prod["m1", "m2"]), 4)
    expect_equal(
        unname(nef_p["m1", "m2"]),
        round(2^(-(0.75 * log2(0.75) + 0.25 * log2(0.25))), 2)
    )
    expect_equal(
        unname(eff_p["m1", "m2"]),
        round(4 * unname(nef_p["m1", "m2"]), 2)
    )
})

test_that("nEffectiveSpecies* values lie in [1, nSpecies] cell-wise", {
    cm <- synCM("test", n_species = 4, max_met = 6, seed = 7)
    a <- SummarizedExperiment::assays(cm)
    n_sp <- as.matrix(a$nSpecies)
    nef_c <- as.matrix(a$nEffectiveSpeciesConsumption)
    nef_p <- as.matrix(a$nEffectiveSpeciesProduction)
    nz <- n_sp != 0
    expect_true(all(nef_c[nz] >= 1 - 1e-8))
    expect_true(all(nef_p[nz] >= 1 - 1e-8))
    expect_true(all(nef_c[nz] <= n_sp[nz] + 1e-8))
    expect_true(all(nef_p[nz] <= n_sp[nz] + 1e-8))
})
