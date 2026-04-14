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

test_that(".decodeBiggEscapes decodes __NN__ sequences", {
    # (40) = '(', (41) = ')'
    expect_equal(.decodeBiggEscapes("acisnzd__40__e__41__"), "acisnzd(e)")
    expect_equal(.decodeBiggEscapes("foo__40__bar__41__"), "foo(bar)")
    expect_equal(.decodeBiggEscapes("no_escapes_here"), "no_escapes_here")
    expect_equal(
        .decodeBiggEscapes(c("a__40__b", "plain", "x__41__y")),
        c("a(b", "plain", "x)y")
    )
})
