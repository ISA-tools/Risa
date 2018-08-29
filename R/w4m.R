# vi: fdm=marker

# Constants {{{1
################################################################

.NA.STRINGS <- c('NA', '')

# Get study {{{1
################################################################

.get.study <- function(isa, study.filename = NULL) {

	if ( ! is.null(study.filename)) {
		if ( ! study.filename %in% isa@study.filenames)
			stop(paste("Cannot find study \"", study.filename, "\".", sep = ''))
		study.name <- isa@study.identifiers[[which(study.filename == isa@study.filenames)]]
	} else { # Take the first study
		cat("No study name set, choose first study.\n")
		study.name <- isa@study.identifiers[[1]]
	}
	if (length(study.name) != 1)
		stop("No study found.")
	cat(paste("Study \"", study.name, "\" has been selected.\n", sep = ""))

	study.df <- isa@study.files[[study.name]]

	return(list(name = study.name, df = study.df))
}

# Get number of assays {{{1
################################################################

.get.nb.assays <- function(isa, study.name) {
	study.assays <- isa@assay.files.per.study[[study.name]]
	return(length(study.assays))
}

# Get chosen assay index {{{1
################################################################

.get.chosen.assay.index <- function(isa, study.name, assay.filename = NULL) {

	study.assays <- isa@assay.files.per.study[[study.name]]
	if (length(study.assays) == 0)
		stop(paste("No assay in study ", study.name, ".", sep = ''))
	if ( ! is.null(assay.filename)) {
		assay.index <- which(assay.filename == isa@assay.filenames)
		if (length(assay.index) == 0)
			stop(paste("Found no assay file \"", assay.filename, "\" in study \"", study.name, "\".", sep = ''))
		if (length(assay.index) > 1)
			stop(paste("Found more than one assay file named \"", assay.filename, "\" in study \"", study.name, "\".", sep = ''))
	} else
		assay.index <- 1

	return(assay.index)
}

# Get assay {{{1
################################################################

.get.assay <- function(isa, study.name, assay.index) {

	study.assay.df <- isa@assay.files.per.study[[study.name]][[assay.index]]
	study.assay.filename <- isa@assay.filenames.per.study[[study.name]][[assay.index]]

	return(list(filename = study.assay.filename, df = study.assay.df))
}

# Make sample metadata {{{1
################################################################

.make.sample.metadata <- function(study, assay, sample.names, normalize = TRUE) {

	# Create sample metadata by merging assay and study metadata
	colnames(study$df) <- make.names(colnames(study$df), uniq = TRUE)
	colnames(assay$df) <- make.names(colnames(assay$df), uniq = TRUE)
	sample.metadata <- merge(assay$df, study$df, by = "Sample.Name", sort = FALSE)

	# Normalize
	if (normalize) {
		sample.metadata <- cbind(data.frame(sample.name = make.names(sample.names, uniq = TRUE), stringsAsFactors = FALSE), sample.metadata)
		colnames(sample.metadata) <- make.names(colnames(sample.metadata), uniq = TRUE)
	}

	return(sample.metadata)
}

# Get measures {{{1
################################################################

.get.measures <- function(isa, assay) {

	measures <- NULL

	if (assay$filename %in% names(isa@maf.filenames.per.assay.filename)) {

		maf_files <- isa@maf.filenames.per.assay.filename[[assay$filename]]
		if ( ! is.null(maf_files)) {

			if (length(maf_files) != 1)
				stop(paste("More than one metabolite assignement file found in assay \"", assay$filename, "\": ", paste(maf_files, collapse = ", "), ".", sep = ''))

			measures <- list(df = isa@maf.dataframes[[maf_files[[1]]]], file = maf_files[[1]])
		}
	}

	return(measures)
}

# Get sample names {{{1
################################################################

.get.sample.names <- function (assay, measures) {

	sample.names <- NULL
	sample.names.field <- NULL
	for (sample.field in colnames(assay$df)) {
		if (all(assay$df[[sample.field]] %in% colnames(measures$df)) && all(! duplicated(assay$df[[sample.field]]))) {
			sample.names.field <- sample.field
			break
		}
	}
	if (is.null(sample.names.field))
		stop(paste("Impossible to find a column for sample names. Either such a column does not exist, or it contains duplicates.", sep = ''))
	sample.names <- assay$df[[sample.names.field]]

	return(sample.names)
}

# Make variable metadata {{{1
################################################################

.make.variable.metadata <- function(measures, sample.names, variable.names, normalize = TRUE) {

	variable.metadata <- measures$df[ ! colnames(measures$df) %in% sample.names]

	# Add variable names as columns
	variable.metadata <- cbind(data.frame(variable.name = variable.names, stringsAsFactors = FALSE), variable.metadata)

	# Normalize
	if (normalize)
		colnames(variable.metadata) <- make.names(colnames(variable.metadata), uniq = TRUE)

	return(variable.metadata)
}

# Make matrix {{{1
################################################################

.make.matrix <- function(measures, sample.names, variable.names, normalize = TRUE) {

	if (any( ! sample.names %in% colnames(measures$df)))
		stop(paste("Cannot find sample names in column names of \"", measures$file, "\".", sep = ''))
	sample.variable.matrix <- measures$df[sample.names]

	# Add variable names as columns
	sample.variable.matrix <- cbind(data.frame(variable.name = variable.names, stringsAsFactors = FALSE), sample.variable.matrix)

	# Normalize sample names
	if (normalize)
		colnames(sample.variable.matrix) <- c('variable.name', make.names(sample.names, uniq = TRUE))

	return(sample.variable.matrix)
}

# Make variable names {{{1
################################################################

.make.variable.names <- function(measures) {

	variable.names <- character(nrow(measures$df))
	for (f in c('mass_to_charge', 'retention_time'))
		if (f %in% colnames(measures$df)) {
			col.vals <- as.character(measures$df[[f]])
			col.vals <- ifelse(is.na(col.vals), '', col.vals)
			variable.names <- paste(variable.names, ifelse(nchar(variable.names) > 0 & nchar(col.vals) > 0, '_', ''), col.vals, sep = '')
		}
	variable.names <- make.names(variable.names, uniq = TRUE)

	return(variable.names)
}

# Convert W4M 3 files format to ISA {{{1
################################################################

w4m2isa <- function(isa, w4m, study.filename = NULL, assays = NULL, assay.filename = NULL) {

	# An existing ISA object will be used for updating.

	# Get study & assays
	study <- .get.study(isa, study.filename = study.filename)
	assay.indices <- if (is.null(assays)) seq(.get.nb.assays(isa, study$name)) else .get.chosen.assay.index(isa, study.name = study$name, assay.filename = assay.filename)
	if (length(assay.indices) != length(w4m))
		stop(paste0(length(assay.indices), ' assays selected, but only ', length(w4m), ' W4M inputs provided.'))

	# Loop on assay indices
	for (assay.index in (assay.indices)) {
	}

	return(isa)
}

# Convert ISA to W4M 3 files format {{{1
################################################################

isa2w4m <- function(isa, study.filename = NULL, assays = NULL, assay.filename = NULL, drop = TRUE) {

	output <- NULL

	# Get study & assays
	study <- .get.study(isa, study.filename = study.filename)
	assay.indices <- if (is.null(assays)) seq(.get.nb.assays(isa, study$name)) else .get.chosen.assay.index(isa, study.name = study$name, assay.filename = assay.filename)

	# Loop on assay indices to extract
	for (assay.index in (assay.indices)) {

		# Get assay
		assay <- .get.assay(isa, study$name, assay.index)

		# Get measures
		measures <- .get.measures(isa, assay)
		if (is.null(measures))
			next

		# Create variable names
		variable.names <- .make.variable.names(measures)

		# Get sample names
		sample.names <- .get.sample.names(assay, measures)

		# Extract sample metadata
		sample.metadata <- .make.sample.metadata(study, assay, sample.names = sample.names, normalize = TRUE)

		# Extract variable metadata
		variable.metadata <- .make.variable.metadata(measures = measures, sample.names = sample.names, variable.names = variable.names, normalize = TRUE)

		# Extract matrix
		sample.variable.matrix <- .make.matrix(measures = measures, sample.names = sample.names, variable.names = variable.names, normalize = TRUE)

		# Build output
		x <- list(samp = sample.metadata, var = variable.metadata, mat = sample.variable.matrix)
		output <- c(output, list(x))
	}

	# Drop
	if (drop && ! is.null(output) && is.null(names(output)) && length(output))
		output <- output[[1]]

	return(output)
}
