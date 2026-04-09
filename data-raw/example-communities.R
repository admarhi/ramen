## Build misosoup24 package data from inst/extdata/misosoup/ CSVs
## Each CSV has columns: metabolite, species, flux

misosoup24 <- list()
for (f in sort(list.files("inst/extdata/misosoup", full.names = TRUE))) {
    name <- tools::file_path_sans_ext(basename(f))
    tb <- tibble::as_tibble(read.csv(f))
    misosoup24[[name]] <- tb
}
usethis::use_data(misosoup24, overwrite = TRUE)
