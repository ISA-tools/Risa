# vi: fdm=marker

context('Basic tests.')

# Constants {{{1
################################################################

TEST.DIR <- file.path(getwd(), '..')
RES.DIR  <- file.path(TEST.DIR, 'res')
OUT.DIR  <- file.path(TEST.DIR, 'output')

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

# Test isa2w4m on MTBLS404, with normalization {{{1
################################################################

test_isa2w4m_mtbls404_normalize <- function() {

	# Load ISA
	isa <- readISAtab(file.path(RES.DIR, 'MTBLS404'), na.strings = c('', 'NA'))
	testthat::expect_is(isa, "ISATab")

	# Convert to W4M
	w4m <- isa2w4m(isa, normalize = TRUE)
	testthat::expect_is(w4m, 'list')
	testthat::expect_length(w4m, 3)
	testthat::expect_false(is.null(names(w4m)))
	testthat::expect_true(all(c('samp', 'var', 'mat') %in% names(w4m)))

	# Save W4M data frames
	for (x in c('samp', 'var', 'mat'))
		write.table(w4m[[x]], file = file.path(TEST.DIR, paste0('test_isa2w4m_mtbls404-', x, '_normalized.tsv')), row.names = FALSE, sep = "\t")

	# Load expected outputs
	samp <- read.table(file = file.path(RES.DIR, 'MTBLS404', 'MTBLS404-w4m-sample-metadata_normalized.tsv'), header = TRUE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE, comment.char = '')
	var <- read.table(file = file.path(RES.DIR, 'MTBLS404', 'MTBLS404-w4m-variable-metadata_normalized.tsv'), header = TRUE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE, comment.char = '')
	mat <- read.table(file = file.path(RES.DIR, 'MTBLS404', 'MTBLS404-w4m-sample-variable-matrix_normalized.tsv'), header = TRUE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE, comment.char = '')

	# Convert factors to strings
	out.var <- data.frame(lapply(w4m$var, function(v) if (is.factor(v)) as.character(v) else v), stringsAsFactors = FALSE)

	# Compare outputs
	testthat::expect_identical(w4m$samp, samp)
	testthat::expect_identical(out.var, var)
	testthat::expect_identical(w4m$mat, mat)
}

# Test isa2w4m on MTBLS404, without normalization {{{1
################################################################

test_isa2w4m_mtbls404_dont_normalize <- function() {

	# Load ISA
	isa <- readISAtab(file.path(RES.DIR, 'MTBLS404'), na.strings = c('', 'NA'))
	testthat::expect_is(isa, "ISATab")

	# Convert to W4M
	w4m <- isa2w4m(isa, normalize = FALSE)
	testthat::expect_is(w4m, 'list')
	testthat::expect_length(w4m, 3)
	testthat::expect_false(is.null(names(w4m)))
	testthat::expect_true(all(c('samp', 'var', 'mat') %in% names(w4m)))

	# Save W4M data frames
	for (x in c('samp', 'var', 'mat'))
		write.table(w4m[[x]], file = file.path(TEST.DIR, paste0('test_isa2w4m_mtbls404-', x, '_not_normalized.tsv')), row.names = FALSE, sep = "\t")

	# Check column names, no names must end with .x, .y, .[0-9]
	for (x in c('samp', 'var', 'mat'))
		testthat::expect_length(grep('^.*\\.([0-9xy])$', colnames(w4m[[x]])), 0)

	# Load expected outputs
	samp <- read.table(file = file.path(RES.DIR, 'MTBLS404', 'MTBLS404-w4m-sample-metadata_not_normalized.tsv'), header = TRUE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE, comment.char = '')
	var <- read.table(file = file.path(RES.DIR, 'MTBLS404', 'MTBLS404-w4m-variable-metadata_not_normalized.tsv'), header = TRUE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE, comment.char = '')
	mat <- read.table(file = file.path(RES.DIR, 'MTBLS404', 'MTBLS404-w4m-sample-variable-matrix_not_normalized.tsv'), header = TRUE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE, comment.char = '')

	# Convert factors to strings
	out.var <- data.frame(lapply(w4m$var, function(v) if (is.factor(v)) as.character(v) else v), stringsAsFactors = FALSE)
	colnames(out.var) <- colnames(w4m$var)

	# Compare outputs
	testthat::expect_identical(w4m$samp, samp)
	testthat::expect_identical(out.var, var)
	testthat::expect_identical(w4m$mat, mat)
}

# Test ISA writing {{{1
################################################################

test_isa_writing <- function() {

	# Load ISA
	isa <- readISAtab(file.path(RES.DIR, 'MTBLS404'), na.strings = c('', 'NA'))
	testthat::expect_is(isa, "ISATab")

	# Write ISA
	output.dir <- file.path(OUT.DIR, 'MTBLS404')
	dir.create(output.dir, recursive = TRUE, showWarnings = FALSE)
	write.ISAtab(isa, output.dir)
	# TODO test that content of output files is the same as content of reference files (even if they are not written exactly the same).
}

# Test w4m2isa {{{1
################################################################

test_w4m2isa <- function() {

	# Load ISA
	isa.ref <- readISAtab(file.path(RES.DIR, 'MTBLS404'), na.strings = c('', 'NA'))
	testthat::expect_is(isa.ref, "ISATab")

	# Convert to W4M
	w4m <- isa2w4m(isa.ref)
	# TODO test also with drop = FALSE

	# Convert back to ISA
	isa <- w4m2isa(isa.ref, w4m)

	# Test that ISA has not changed
	study.name <- isa.ref@study.identifiers[[1]]
	testthat::expect_identical(isa.ref@study.files[[study.name]], isa@study.files[[study.name]])
	testthat::expect_identical(isa.ref@assay.files.per.study[[study.name]][[1]], isa@assay.files.per.study[[study.name]][[1]])
	assay.filename <- isa@assay.filenames.per.study[[study.name]][[1]]
	maf.files <- isa@maf.filenames.per.assay.filename[[assay.filename]]
	testthat::expect_identical(isa.ref@maf.dataframes[[maf.files[[1]]]], isa@maf.dataframes[[maf.files[[1]]]])
}

# Test w4m2isa with added columns {{{1
################################################################

test_w4m2isa_added_cols <- function() {

	# Load ISA
	isa.ref <- readISAtab(file.path(RES.DIR, 'MTBLS404'), na.strings = c('', 'NA'))
	testthat::expect_is(isa.ref, "ISATab")

	# Convert to W4M
	w4m <- isa2w4m(isa.ref)

	# Add columns
	new.samp.col <- 'new.col'
	w4m$samp[[new.samp.col]] <- 1

	# Convert back to ISA
	isa <- w4m2isa(isa.ref, w4m)
	study.name <- isa@study.identifiers[[1]]
	testthat::expect_identical(isa.ref@study.files[[study.name]], isa@study.files[[study.name]])
	testthat::expect_true(new.samp.col %in% colnames(isa@assay.files.per.study[[study.name]][[1]]))
	assay.filename <- isa@assay.filenames.per.study[[study.name]][[1]]
	maf.files <- isa@maf.filenames.per.assay.filename[[assay.filename]]
	testthat::expect_identical(isa.ref@maf.dataframes[[maf.files[[1]]]], isa@maf.dataframes[[maf.files[[1]]]])
}

# Test w4m2isa with removed samples {{{1
################################################################

test_w4m2isa_with_removed_samples <- function() {

	# Load ISA
	isa.ref <- readISAtab(file.path(RES.DIR, 'MTBLS404'), na.strings = c('', 'NA'))
	testthat::expect_is(isa.ref, "ISATab")

	# Convert to W4M
	w4m <- isa2w4m(isa.ref)

	# Remove samples
	samp.name.col <- 'Sample Name'
	testthat::expect_true(samp.name.col %in% colnames(w4m$samp))
	samp.to.remove <- w4m$samp[[samp.name.col]][[5]]
	w4m$samp <- w4m$samp[w4m$samp[[samp.name.col]] == samp.to.remove, ]
	w4m$mat <- w4m$mat[ ! colnames(w4m$mat) %in% samp.to.remove]

	# Convert back to ISA
	isa <- w4m2isa(isa.ref, w4m)
	study.name <- isa@study.identifiers[[1]]
	testthat::expect_identical(isa.ref@study.files[[study.name]], isa@study.files[[study.name]])
	testthat::expect_true(nrow(isa@assay.files.per.study[[study.name]][[1]]) == nrow(isa.ref@assay.files.per.study[[study.name]][[1]]) - 1)
	assay.filename <- isa@assay.filenames.per.study[[study.name]][[1]]
	maf.files <- isa@maf.filenames.per.assay.filename[[assay.filename]]
	testthat::expect_false(samp.to.remove %in% colnames(isa@maf.dataframes[[maf.files[[1]]]]))
}

# Main {{{1
################################################################

test_that("Load ISA.", test_load_faahko_isa())
test_that("We can build an XCMS set.", test_build_xcms_set_from_faahko())
test_that("Conversion from ISA to W4M format for faahKO fails.", test_isa2w4m_faahko())
test_that("Conversion from MTBLS404 ISA to W4M format works, with normalization.", test_isa2w4m_mtbls404_normalize())
test_that("Conversion from MTBLS404 ISA to W4M format works, without normalization.", test_isa2w4m_mtbls404_dont_normalize())
test_that("We can write ISA ?_*.txt files on disk.", test_isa_writing())
test_that("w4m2isa runs on unmodified W4M data frames gives back the same ISA data frames.", test_w4m2isa())
test_that("w4m2isa runs correctly when columns have been added.", test_w4m2isa_added_cols())
test_that("w4m2isa runs correctly when samples have been removed.", test_w4m2isa_with_removed_samples())
# TODO test w4m2isa when samples have been removed.
# TODO test w4m2isa when variables have been removed.
# TODO test w4m2isa when both samples and variables have been removed.
