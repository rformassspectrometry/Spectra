library("testthat")
library("Spectra")

sciex_file <- dir(system.file("sciex", package = "msdata"), full.names = TRUE)
cdf_file <- dir(system.file("cdf", package = "msdata"), full.names = TRUE)

test_check("Spectra")
