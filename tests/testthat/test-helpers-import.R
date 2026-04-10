test_that("import helper extracts BiGG metabolite IDs correctly", {
    expect_equal(.normalizeBiggIds("R_EX_ac_e"), "ac")
    expect_equal(.normalizeBiggIds("EX_glc__D_e"), "glc__D")
    expect_equal(.normalizeBiggIds("M_ac_e"), "ac")
    expect_equal(.normalizeBiggIds("EX_lac_D(e)"), "lac_D")
    expect_equal(.normalizeBiggIds("ac_e"), "ac")
    expect_equal(.normalizeBiggIds("ac"), "ac")
    expect_equal(
        .normalizeBiggIds(c(
            "R_EX_ac_e",
            "EX_glc__D_e",
            "M_ac_e",
            "EX_lac_D(e)",
            "ac_e",
            "ac"
        )),
        c("ac", "glc__D", "ac", "lac_D", "ac", "ac")
    )
})
