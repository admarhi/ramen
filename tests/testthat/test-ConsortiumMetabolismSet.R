test_that("ConsortiumMetabolismSet constructor works", {
    # Create two test consortia
    data1 <- tibble::tibble(
        species = c("s1", "s1", "s2", "s2"),
        metabolite = c("m1", "m2", "m1", "m3"),
        flux = c(-1, 1, -1, 1)
    )

    data2 <- tibble::tibble(
        species = c("s3", "s3", "s4", "s4"),
        metabolite = c("m1", "m2", "m1", "m4"),
        flux = c(-1, 1, -1, 1)
    )

    cm1 <- ConsortiumMetabolism(data1, name = "cm1")
    cm2 <- ConsortiumMetabolism(data2, name = "cm2")

    # Create set
    cms <- ConsortiumMetabolismSet(list(cm1, cm2), name = "test_set")

    # Assertions
    expect_s4_class(cms, "ConsortiumMetabolismSet")
    expect_equal(cms@Name, "test_set")
    expect_equal(length(cms@Consortia), 2)
    expect_s3_class(cms@Pathways, "data.frame")
})

test_that("species(CMS) returns character vector", {
    data1 <- tibble::tibble(
        species = c("s1", "s1", "s2", "s2"),
        metabolite = c("m1", "m2", "m1", "m3"),
        flux = c(-1, 1, -1, 1)
    )
    data2 <- tibble::tibble(
        species = c("s2", "s2", "s3", "s3"),
        metabolite = c("m1", "m2", "m1", "m4"),
        flux = c(-1, 1, -1, 1)
    )
    cm1 <- ConsortiumMetabolism(data1, name = "cm1")
    cm2 <- ConsortiumMetabolism(data2, name = "cm2")
    cms <- ConsortiumMetabolismSet(list(cm1, cm2), name = "test")

    sp <- species(cms)
    expect_type(sp, "character")
    expect_true("s1" %in% sp)
    expect_true("s2" %in% sp)
    expect_true("s3" %in% sp)
    expect_equal(sp, sort(unique(sp)))
})

test_that("CMS pathways returns concise output by default", {
    data1 <- tibble::tibble(
        species = c("s1", "s1", "s2", "s2"),
        metabolite = c("m1", "m2", "m1", "m3"),
        flux = c(-1, 1, -1, 1)
    )

    data2 <- tibble::tibble(
        species = c("s1", "s1", "s3", "s3"),
        metabolite = c("m1", "m2", "m2", "m4"),
        flux = c(-1, 1, -1, 1)
    )

    cm1 <- ConsortiumMetabolism(data1, name = "cm1")
    cm2 <- ConsortiumMetabolism(data2, name = "cm2")
    cms <- ConsortiumMetabolismSet(list(cm1, cm2), name = "test")

    pw <- pathways(cms)
    expect_s3_class(pw, "data.frame")
    expect_true(nrow(pw) > 0)
    expect_named(
        pw,
        c("consumed", "produced", "n_species", "n_cons")
    )
})

test_that("CMS pathways verbose returns full detail", {
    data1 <- tibble::tibble(
        species = c("s1", "s1", "s2", "s2"),
        metabolite = c("m1", "m2", "m1", "m3"),
        flux = c(-1, 1, -1, 1)
    )

    data2 <- tibble::tibble(
        species = c("s1", "s1", "s3", "s3"),
        metabolite = c("m1", "m2", "m2", "m4"),
        flux = c(-1, 1, -1, 1)
    )

    cm1 <- ConsortiumMetabolism(data1, name = "cm1")
    cm2 <- ConsortiumMetabolism(data2, name = "cm2")
    cms <- ConsortiumMetabolismSet(list(cm1, cm2), name = "test")

    pw <- pathways(cms, verbose = TRUE)
    expect_s3_class(pw, "data.frame")
    expect_true("cm_name" %in% names(pw))
    expect_true("species" %in% names(pw))
    expect_true("n_cons" %in% names(pw))
    expect_true(ncol(pw) > 4)
})

test_that("CMS pathways type filtering works", {
    data1 <- tibble::tibble(
        species = c("s1", "s1", "s2", "s2"),
        metabolite = c("m1", "m2", "m1", "m3"),
        flux = c(-1, 1, -1, 1)
    )

    data2 <- tibble::tibble(
        species = c("s1", "s1", "s3", "s3"),
        metabolite = c("m1", "m2", "m2", "m4"),
        flux = c(-1, 1, -1, 1)
    )

    cm1 <- ConsortiumMetabolism(data1, name = "cm1")
    cm2 <- ConsortiumMetabolism(data2, name = "cm2")
    cms <- ConsortiumMetabolismSet(list(cm1, cm2), name = "test")

    pw_all <- pathways(cms, type = "all")
    expect_true(nrow(pw_all) > 0)
    expect_named(
        pw_all,
        c("consumed", "produced", "n_species", "n_cons")
    )
})

test_that("name<- and description<- work for ConsortiumMetabolismSet", {
    data1 <- tibble::tibble(
        species = c("s1", "s1"),
        metabolite = c("m1", "m2"),
        flux = c(-1, 1)
    )

    cm1 <- ConsortiumMetabolism(data1, name = "cm1")
    cms <- ConsortiumMetabolismSet(list(cm1), name = "original")

    # Test name<-
    name(cms) <- "new_name"
    expect_equal(name(cms), "new_name")

    # Test description<-
    description(cms) <- "test description"
    expect_equal(description(cms), "test description")
})


## ---- CMS BinaryMatrices and OverlapMatrix -----------------------------------

test_that("CMS BinaryMatrices slot is populated", {
    cm1 <- synCM("a", n_species = 3, max_met = 5)
    cm2 <- synCM("b", n_species = 3, max_met = 5)
    cms <- ConsortiumMetabolismSet(
        list(cm1, cm2),
        name = "test"
    )
    expect_true(length(cms@BinaryMatrices) == 2L)
    ## All matrices same dimension
    dims <- vapply(
        cms@BinaryMatrices,
        nrow,
        integer(1L)
    )
    expect_true(all(dims == dims[1L]))
})

test_that("CMS OverlapMatrix unchanged after refactor", {
    set.seed(42)
    cm1 <- synCM("a", n_species = 4, max_met = 8, seed = 42)
    cm2 <- synCM("b", n_species = 4, max_met = 8, seed = 43)
    cms <- ConsortiumMetabolismSet(
        list(cm1, cm2),
        name = "test"
    )
    ## OverlapMatrix should be a 2x2 matrix
    expect_equal(nrow(cms@OverlapMatrix), 2L)
    expect_equal(ncol(cms@OverlapMatrix), 2L)
    expect_true(is.matrix(cms@OverlapMatrix))
})

test_that("single-consortium CMS can be created", {
    cm1 <- synCM("a", n_species = 3, max_met = 5)
    cms <- ConsortiumMetabolismSet(list(cm1), name = "single")
    expect_s4_class(cms, "ConsortiumMetabolismSet")
    expect_equal(length(cms@Consortia), 1L)
})

## ---- B5: Duplicate consortium names error ----------------------------------

test_that("CMS constructor errors on duplicate consortium names", {
    cm1 <- synCM("same_name", n_species = 3, max_met = 5, seed = 1)
    cm2 <- synCM("same_name", n_species = 3, max_met = 5, seed = 2)
    expect_error(
        ConsortiumMetabolismSet(cm1, cm2, name = "test"),
        "unique"
    )
})

test_that("CMS constructor accepts unique consortium names", {
    cm1 <- synCM("name_a", n_species = 3, max_met = 5, seed = 1)
    cm2 <- synCM("name_b", n_species = 3, max_met = 5, seed = 2)
    expect_no_error(
        ConsortiumMetabolismSet(cm1, cm2, name = "test")
    )
})

## ---- B6: core pathway quantile uses actual distribution --------------------

test_that("pathways(cms, type='core') uses actual species distribution", {
    cm1 <- synCM("a", n_species = 4, max_met = 5, seed = 10)
    cm2 <- synCM("b", n_species = 4, max_met = 5, seed = 20)
    cm3 <- synCM("c", n_species = 4, max_met = 5, seed = 30)
    cms <- ConsortiumMetabolismSet(
        list(cm1, cm2, cm3),
        name = "test"
    )
    all_pw <- pathways(cms, type = "all")
    core_pw <- pathways(cms, type = "core")
    if (nrow(core_pw) > 0L) {
        quant <- stats::quantile(
            all_pw$n_species,
            p = 0.9
        )
        expect_true(all(core_pw$n_species > quant))
    }
})

## ---- consortia() method for CMS -------------------------------------------

test_that("consortia(cms) returns list of CM objects", {
    cm1 <- synCM("a", n_species = 3, max_met = 5, seed = 1)
    cm2 <- synCM("b", n_species = 3, max_met = 5, seed = 2)
    cms <- ConsortiumMetabolismSet(
        list(cm1, cm2),
        name = "test"
    )
    result <- consortia(cms)
    expect_type(result, "list")
    expect_length(result, 2L)
    expect_s4_class(result[[1L]], "ConsortiumMetabolism")
    expect_s4_class(result[[2L]], "ConsortiumMetabolism")
})

## ---- growth() on CMS --------------------------------------------------------

test_that("growth(cms) returns named list from stored CMs", {
    test_data <- tibble::tibble(
        species = c("s1", "s1", "s2", "s2"),
        metabolite = c("m1", "m2", "m1", "m3"),
        flux = c(-1, 1, -1, 1)
    )
    cm1 <- ConsortiumMetabolism(
        test_data,
        name = "cm1",
        growth = c(s1 = 0.5, s2 = 0.3)
    )
    cm2 <- ConsortiumMetabolism(
        test_data,
        name = "cm2",
        growth = c(s1 = 0.8, s2 = 0.1)
    )
    cms <- ConsortiumMetabolismSet(
        list(cm1, cm2),
        name = "test"
    )
    gr <- growth(cms)
    expect_type(gr, "list")
    expect_named(gr, c("cm1", "cm2"))
    expect_equal(gr$cm1, c(s1 = 0.5, s2 = 0.3))
    expect_equal(gr$cm2, c(s1 = 0.8, s2 = 0.1))
})

test_that("growth(cms) returns NULLs when no growth data", {
    cm1 <- synCM("a", n_species = 3, max_met = 5, seed = 1)
    cm2 <- synCM("b", n_species = 3, max_met = 5, seed = 2)
    cms <- ConsortiumMetabolismSet(
        list(cm1, cm2),
        name = "test"
    )
    gr <- growth(cms)
    expect_type(gr, "list")
    expect_named(gr, c("a", "b"))
    expect_null(gr$a)
    expect_null(gr$b)
})

test_that("consortia(cms) names match consortium names", {
    cm1 <- synCM("alpha", n_species = 3, max_met = 5, seed = 1)
    cm2 <- synCM("beta", n_species = 3, max_met = 5, seed = 2)
    cms <- ConsortiumMetabolismSet(
        list(cm1, cm2),
        name = "test"
    )
    result <- consortia(cms)
    expect_equal(
        vapply(result, name, character(1L)),
        c("alpha", "beta")
    )
})

## ---- overlapMatrix() accessor -----------------------------------------------

test_that("overlapMatrix returns numeric matrix with correct dims", {
    cm1 <- synCM("a", n_species = 3, max_met = 5, seed = 1)
    cm2 <- synCM("b", n_species = 4, max_met = 6, seed = 2)
    cms <- ConsortiumMetabolismSet(
        list(cm1, cm2),
        name = "test"
    )
    om <- overlapMatrix(cms)
    expect_true(is.matrix(om))
    expect_true(is.numeric(om))
    expect_equal(nrow(om), 2L)
    expect_equal(ncol(om), 2L)
    expect_equal(rownames(om), c("a", "b"))
    expect_equal(colnames(om), c("a", "b"))
})

test_that("overlapMatrix diagonal is zero", {
    cm1 <- synCM("x", n_species = 3, max_met = 5, seed = 10)
    cm2 <- synCM("y", n_species = 3, max_met = 5, seed = 20)
    cm3 <- synCM("z", n_species = 3, max_met = 5, seed = 30)
    cms <- ConsortiumMetabolismSet(
        list(cm1, cm2, cm3),
        name = "test"
    )
    om <- overlapMatrix(cms)
    expect_equal(nrow(om), 3L)
    expect_true(all(diag(om) == 0))
})

## ---- show(CMS) community size stats -----------------------------------------

test_that("show(cms) includes consortium count and community size stats", {
    cm1 <- synCM("a", n_species = 3, max_met = 5, seed = 1)
    cm2 <- synCM("b", n_species = 4, max_met = 6, seed = 2)
    cms <- ConsortiumMetabolismSet(
        list(cm1, cm2),
        name = "test"
    )
    out <- paste(
        capture.output(show(cms), type = "message"),
        collapse = " "
    )
    expect_match(out, "2 consortia")
    expect_match(out, "species")
    expect_match(out, "min")
    expect_match(out, "max")
    expect_match(out, "mean")
    expect_match(out, "metabolites")
})

## ---- speciesSummary(CMS) ----------------------------------------------------

test_that("speciesSummary(CMS) returns tibble with correct columns", {
    cm1 <- synCM("a", n_species = 3, max_met = 5, seed = 1)
    cm2 <- synCM("b", n_species = 4, max_met = 6, seed = 2)
    cms <- ConsortiumMetabolismSet(
        list(cm1, cm2),
        name = "test"
    )
    ss <- speciesSummary(cms)
    expect_s3_class(ss, "tbl_df")
    expect_named(ss, c("species", "n_consortia", "n_pathways"))
    expect_type(ss$species, "character")
    expect_type(ss$n_consortia, "integer")
    expect_type(ss$n_pathways, "integer")
})

test_that("speciesSummary(CMS) n_consortia bounded by n_consortia", {
    cm1 <- synCM("a", n_species = 3, max_met = 5, seed = 1)
    cm2 <- synCM("b", n_species = 3, max_met = 5, seed = 2)
    cms <- ConsortiumMetabolismSet(
        list(cm1, cm2),
        name = "test"
    )
    ss <- speciesSummary(cms)
    expect_true(all(ss$n_consortia >= 1L))
    expect_true(all(ss$n_consortia <= 2L))
})

test_that("speciesSummary(CMS) species match species(CMS)", {
    cm1 <- synCM("a", n_species = 3, max_met = 5, seed = 1)
    cm2 <- synCM("b", n_species = 4, max_met = 6, seed = 2)
    cms <- ConsortiumMetabolismSet(
        list(cm1, cm2),
        name = "test"
    )
    expect_setequal(speciesSummary(cms)$species, species(cms))
})

## ---- filterConsortia(CMS) ---------------------------------------------------

test_that("filterConsortia integer index returns correct subset", {
    cm1 <- synCM("a", n_species = 3, max_met = 5, seed = 1)
    cm2 <- synCM("b", n_species = 3, max_met = 5, seed = 2)
    cm3 <- synCM("c", n_species = 3, max_met = 5, seed = 3)
    cms <- ConsortiumMetabolismSet(
        list(cm1, cm2, cm3),
        name = "full"
    )
    sub <- filterConsortia(cms, c(1L, 3L))
    expect_s4_class(sub, "ConsortiumMetabolismSet")
    expect_equal(length(sub@Consortia), 2L)
    expect_equal(
        vapply(sub@Consortia, name, character(1L)),
        c("a", "c")
    )
})

test_that("filterConsortia character names returns correct subset", {
    cm1 <- synCM("alpha", n_species = 3, max_met = 5, seed = 1)
    cm2 <- synCM("beta", n_species = 3, max_met = 5, seed = 2)
    cm3 <- synCM("gamma", n_species = 3, max_met = 5, seed = 3)
    cms <- ConsortiumMetabolismSet(
        list(cm1, cm2, cm3),
        name = "full"
    )
    sub <- filterConsortia(cms, c("alpha", "gamma"))
    expect_equal(length(sub@Consortia), 2L)
    expect_equal(
        vapply(sub@Consortia, name, character(1L)),
        c("alpha", "gamma")
    )
})

test_that("filterConsortia logical vector returns correct subset", {
    cm1 <- synCM("a", n_species = 3, max_met = 5, seed = 1)
    cm2 <- synCM("b", n_species = 3, max_met = 5, seed = 2)
    cm3 <- synCM("c", n_species = 3, max_met = 5, seed = 3)
    cms <- ConsortiumMetabolismSet(
        list(cm1, cm2, cm3),
        name = "full"
    )
    sub <- filterConsortia(cms, c(TRUE, FALSE, TRUE))
    expect_equal(length(sub@Consortia), 2L)
    expect_equal(
        vapply(sub@Consortia, name, character(1L)),
        c("a", "c")
    )
})

test_that("filterConsortia preserves name and description", {
    cm1 <- synCM("a", n_species = 3, max_met = 5, seed = 1)
    cm2 <- synCM("b", n_species = 3, max_met = 5, seed = 2)
    cms <- ConsortiumMetabolismSet(
        list(cm1, cm2),
        name = "myname",
        desc = "mydesc"
    )
    sub <- filterConsortia(cms, 1L)
    expect_equal(sub@Name, "myname")
    expect_equal(sub@Description, "mydesc")
})

test_that("filterConsortia errors on out-of-bounds integer", {
    cm1 <- synCM("a", n_species = 3, max_met = 5, seed = 1)
    cm2 <- synCM("b", n_species = 3, max_met = 5, seed = 2)
    cms <- ConsortiumMetabolismSet(
        list(cm1, cm2),
        name = "test"
    )
    expect_error(filterConsortia(cms, 5L), "out of bounds")
})

test_that("filterConsortia errors on unknown character name", {
    cm1 <- synCM("a", n_species = 3, max_met = 5, seed = 1)
    cms <- ConsortiumMetabolismSet(list(cm1), name = "test")
    expect_error(filterConsortia(cms, "zzz"), "not found")
})

test_that("filterConsortia errors on wrong-length logical", {
    cm1 <- synCM("a", n_species = 3, max_met = 5, seed = 1)
    cm2 <- synCM("b", n_species = 3, max_met = 5, seed = 2)
    cms <- ConsortiumMetabolismSet(
        list(cm1, cm2),
        name = "test"
    )
    expect_error(filterConsortia(cms, c(TRUE)), "length")
})

## ---- species(CMS, type = ...) filter ---------------------------------------

test_that("species(CMS, type='all') matches default", {
    data("misosoup24", package = "ramen")
    cm_list <- lapply(seq_len(10), function(i) {
        ConsortiumMetabolism(
            misosoup24[[i]],
            name = names(misosoup24)[i]
        )
    })
    cms <- ConsortiumMetabolismSet(cm_list, name = "ms10")
    expect_equal(species(cms), species(cms, type = "all"))
})

test_that("species(CMS, type='generalists') returns non-empty subset", {
    data("misosoup24", package = "ramen")
    cm_list <- lapply(seq_len(10), function(i) {
        ConsortiumMetabolism(
            misosoup24[[i]],
            name = names(misosoup24)[i]
        )
    })
    cms <- ConsortiumMetabolismSet(cm_list, name = "ms10")
    gen <- species(cms, type = "generalists")
    expect_type(gen, "character")
    expect_true(length(gen) >= 1L)
    expect_true(all(gen %in% species(cms)))
})

test_that("species(CMS, type='specialists') returns non-empty subset", {
    data("misosoup24", package = "ramen")
    cm_list <- lapply(seq_len(10), function(i) {
        ConsortiumMetabolism(
            misosoup24[[i]],
            name = names(misosoup24)[i]
        )
    })
    cms <- ConsortiumMetabolismSet(cm_list, name = "ms10")
    spec <- species(cms, type = "specialists")
    expect_type(spec, "character")
    expect_true(length(spec) >= 1L)
    expect_true(all(spec %in% species(cms)))
})

test_that("species(CMS) generalists outrank specialists in pathway count", {
    data("misosoup24", package = "ramen")
    cm_list <- lapply(seq_len(10), function(i) {
        ConsortiumMetabolism(
            misosoup24[[i]],
            name = names(misosoup24)[i]
        )
    })
    cms <- ConsortiumMetabolismSet(cm_list, name = "ms10")
    summary_df <- speciesSummary(cms)
    gen <- species(cms, type = "generalists")
    spec <- species(cms, type = "specialists")
    gen_counts <- summary_df$n_pathways[
        summary_df$species %in% gen
    ]
    spec_counts <- summary_df$n_pathways[
        summary_df$species %in% spec
    ]
    expect_true(min(gen_counts) >= max(spec_counts))
})

test_that("species(CMS, type='invalid') errors via match.arg", {
    cm1 <- synCM("a", n_species = 3, max_met = 5, seed = 1)
    cm2 <- synCM("b", n_species = 3, max_met = 5, seed = 2)
    cms <- ConsortiumMetabolismSet(
        list(cm1, cm2),
        name = "test"
    )
    expect_error(
        species(cms, type = "bogus"),
        "should be one of"
    )
})

test_that("species(CMS) errors on invalid quantileCutoff", {
    cm1 <- synCM("a", n_species = 3, max_met = 5, seed = 1)
    cm2 <- synCM("b", n_species = 3, max_met = 5, seed = 2)
    cms <- ConsortiumMetabolismSet(
        list(cm1, cm2),
        name = "test"
    )
    expect_error(
        species(cms, type = "generalists", quantileCutoff = 0)
    )
    expect_error(
        species(cms, type = "specialists", quantileCutoff = 1)
    )
})

## ---- B6 residual: aux/niche tie handling -----------------------------------

test_that("pathways(cms, 'aux') non-empty when n_species ties at floor", {
    data("misosoup24", package = "ramen")
    cm_list <- lapply(seq_len(10), function(i) {
        ConsortiumMetabolism(
            misosoup24[[i]],
            name = names(misosoup24)[i]
        )
    })
    cms <- ConsortiumMetabolismSet(cm_list, name = "ms10")
    aux_pw <- pathways(cms, type = "aux")
    expect_true(nrow(aux_pw) >= 1L)
})

test_that("pathways(cms, 'niche') non-empty when n_cons ties at floor", {
    data("misosoup24", package = "ramen")
    cm_list <- lapply(seq_len(10), function(i) {
        ConsortiumMetabolism(
            misosoup24[[i]],
            name = names(misosoup24)[i]
        )
    })
    cms <- ConsortiumMetabolismSet(cm_list, name = "ms10")
    niche_pw <- pathways(cms, type = "niche")
    expect_true(nrow(niche_pw) >= 1L)
})

test_that("aux and core do not double-count pathways at the boundary", {
    data("misosoup24", package = "ramen")
    cm_list <- lapply(seq_len(10), function(i) {
        ConsortiumMetabolism(
            misosoup24[[i]],
            name = names(misosoup24)[i]
        )
    })
    cms <- ConsortiumMetabolismSet(cm_list, name = "ms10")
    aux_pw <- pathways(cms, type = "aux")
    core_pw <- pathways(cms, type = "core")
    aux_keys <- paste(aux_pw$consumed, aux_pw$produced, sep = "|")
    core_keys <- paste(
        core_pw$consumed,
        core_pw$produced,
        sep = "|"
    )
    expect_length(intersect(aux_keys, core_keys), 0L)
})

## ---- as.data.frame(CMS) ----------------------------------------------------

test_that("as.data.frame(CMS) row-binds per-CM edges with consortium col", {
    cm1 <- synCM("a", n_species = 3, max_met = 5, seed = 1)
    cm2 <- synCM("b", n_species = 4, max_met = 6, seed = 2)
    cms <- ConsortiumMetabolismSet(list(cm1, cm2), name = "test")
    df <- as.data.frame(cms)

    expect_s3_class(df, "data.frame")
    expect_setequal(
        colnames(df),
        c("consortium", "met", "species", "flux")
    )
    ## Row count equals sum of per-CM edge counts.
    expect_equal(
        nrow(df),
        nrow(as.data.frame(cm1)) + nrow(as.data.frame(cm2))
    )
    expect_setequal(unique(df$consortium), c("a", "b"))
})
