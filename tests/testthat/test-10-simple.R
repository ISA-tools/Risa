# vi: fdm=marker

context('Basic tests.')

# Load faahKO ISA {{{1
################################################################

load_faahko_isa <- function() {
	faahkoISA <- readISAtab(find.package("faahKO"))
	testthat::expect_is(faahkoISA, "ISATab")
}

# Build XCMS set from faahKO {{{1
################################################################

build_xcms_set_from_faahko <- function() {
	faahkoISA <- readISAtab(find.package("faahKO"))
	testthat::expect_is(faahkoISA, "ISATab")
	assay.filename <- faahkoISA["assay.filenames"][1]
	testthat::expect_is(assay.filename, "character")
	faahkoXset <- processAssayXcmsSet(faahkoISA, assay.filename)
	testthat::expect_is(faahkoXset, "xcmsSet")
}

# Main {{{1
################################################################

test_that("Load ISA.", load_faahko_isa())
test_that("We can build an XCMS set.", build_xcms_set_from_faahko())
