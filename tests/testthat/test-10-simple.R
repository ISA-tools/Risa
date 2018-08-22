# vi: fdm=marker

context('Basic tests.')

# Constants {{{1
################################################################

TEST.DIR <- file.path(getwd(), '..')
RES.DIR  <- file.path(TEST.DIR, 'res')

# Load faahKO ISA {{{1
################################################################

test_load_faahko_isa <- function() {
	faahkoISA <- readISAtab(find.package("faahKO"))
	testthat::expect_is(faahkoISA, "ISATab")
}

# Build XCMS set from faahKO {{{1
################################################################

test_build_xcms_set_from_faahko <- function() {
	faahkoISA <- readISAtab(find.package("faahKO"))
	testthat::expect_is(faahkoISA, "ISATab")
	assay.filename <- faahkoISA["assay.filenames"][1]
	testthat::expect_is(assay.filename, "character")
	faahkoXset <- processAssayXcmsSet(faahkoISA, assay.filename)
	testthat::expect_is(faahkoXset, "xcmsSet")
}

# Test isa2w4m on faahKO {{{1
################################################################

test_isa2w4m_faahko <- function() {
	faahkoISA <- readISAtab(find.package("faahKO"))
	testthat::expect_is(faahkoISA, "ISATab")
	w4m <- isa2w4m(faahkoISA)
	testthat::expect_is(w4m, 'NULL') # faahKO contains no measurement file (m_*.txt), so isa2w4m() fails.
}

# Test isa2w4m on MTBLS404 {{{1
################################################################

test_isa2w4m_mtbls404 <- function() {

	# Load ISA
	isa <- readISAtab(file.path(RES.DIR, 'MTBLS404'), na.strings = c('', 'NA'))
	testthat::expect_is(isa, "ISATab")

	# Convert to W4M
	w4m <- isa2w4m(isa)
	testthat::expect_is(w4m, 'list')
	testthat::expect_length(w4m, 3)
	testthat::expect_false(is.null(names(w4m)))
	testthat::expect_true(all(c('samp', 'var', 'mat') %in% names(w4m)))

	# Load expected outputs
	samp <- read.table(file = file.path(RES.DIR, 'MTBLS404', 'MTBLS404-w4m-sample-metadata.tsv'), header = TRUE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE, comment.char = '')
	var <- read.table(file = file.path(RES.DIR, 'MTBLS404', 'MTBLS404-w4m-variable-metadata.tsv'), header = TRUE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE, comment.char = '')
	mat <- read.table(file = file.path(RES.DIR, 'MTBLS404', 'MTBLS404-w4m-sample-variable-matrix.tsv'), header = TRUE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE, comment.char = '')

	# Convert factors to strings
	out.var <- data.frame(lapply(w4m$var, function(v) if (is.factor(v)) as.character(v) else v), stringsAsFactors = FALSE)

	# Compare outputs
	testthat::expect_identical(w4m$samp, samp)
	testthat::expect_identical(out.var, var)
	testthat::expect_identical(w4m$mat, mat)
}

# Main {{{1
################################################################

test_that("Load ISA.", test_load_faahko_isa())
test_that("We can build an XCMS set.", test_build_xcms_set_from_faahko())
test_that("Conversion from ISA to W4M format for faahKO fails.", test_isa2w4m_faahko())
test_that("Conversion from MTBLS404 ISA to W4M format works.", test_isa2w4m_mtbls404())
