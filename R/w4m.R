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
	colnames(study$df) <- make.names(colnames(study$df), uniq = TRUE) # TODO make.names should be optional
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

.make.variable.names <- function(maf.df) {

	variable.names <- character(nrow(maf.df))
	for (f in c('mass_to_charge', 'retention_time'))
		if (f %in% colnames(maf.df)) {
			col.vals <- as.character(maf.df[[f]])
			col.vals <- ifelse(is.na(col.vals), '', col.vals)
			variable.names <- paste(variable.names, ifelse(nchar(variable.names) > 0 & nchar(col.vals) > 0, '_', ''), col.vals, sep = '')
		}
	variable.names <- make.names(variable.names, uniq = TRUE)

	return(variable.names)
}

# Update assay {{{1
################################################################

.update.assay <- function(isa, study.name, assay.index, w4m.samp) {

	# Get study and assay data frame
	study.df <- isa@study.files[[study.name]]
	study.assay.df <- isa@assay.files.per.study[[study.name]][[assay.index]]

	# TODO Remove rows when some samples have been deleted

	# Modify assay data frame
	cols <- colnames(w4m.samp)
	cols <- cols[ ! cols %in% colnames(study.df)]
	cols <- cols[ ! cols %in% colnames(study.assay.df)]
	if ( ! "Sample Name" %in% colnames(study.assay.df))
		stop("Cannot find column \"Sample Name\" into ISA assay data frame.")
 	if ( ! "Sample.Name" %in% colnames(w4m.samp))
		stop("Cannot find column \"Sample.Name\" into W4M sample data frame.")
	if ( ! identical(study.assay.df[["Sample Name"]], w4m.samp[["Sample.Name"]]))
		stop("\"Sample Name\" column of ISA assay data frame and \"Sample.Name\" column of W4M sample data frame aren't identical.")
	study.assay.df <- cbind(study.assay.df, w4m.samp[cols])
#	study.assay.df <- merge(study.assay.df, w4m.samp[cols], by.x = "Sample Name", by.y = "Sample.Name", sort = FALSE)

	# Update assay data frame
	isa@assay.files.per.study[[study.name]][[assay.index]] <- study.assay.df 
}

# Update MAF {{{1
################################################################

.update.maf <- function(isa, study.name, assay.index, w4m.var, w4m.mat) {

	assay.filename <- isa@assay.filenames.per.study[[study.name]][[assay.index]]
	# Get MAF data frame
	if (assay.filename %in% names(isa@maf.filenames.per.assay.filename)) {
		maf.files <- isa@maf.filenames.per.assay.filename[[assay.filename]]
		if ( ! is.null(maf.files)) {

			if (length(maf.files) != 1)
				stop(paste("More than one metabolite assignement file found in assay \"", assay.filename, "\": ", paste(maf.files, collapse = ", "), ".", sep = ''))

			# Get MAF data frame
			maf.df <- isa@maf.dataframes[[maf.files[[1]]]]

			# Add variable names column
			maf.df["Variable Name"] <- .make.variable.names(maf.df)

			# Remove sample columns
			# TODO

			# Remove rows
			# TODO

			# Add new columns from W4M variable data frame
			cols <- colnames(w4m.var)
			cols <- cols[ ! cols %in% colnames(maf.df)]
			if ( ! "Variable Name" %in% colnames(maf.df))
				stop("Cannot find column \"Variable Name\" into ISA MAF data frame.")
 			if ( ! "variable.name" %in% colnames(w4m.var))
				stop("Cannot find column \"variable.name\" into W4M variable data frame.")
			if ( ! identical(maf.df[["Variable Name"]], w4m.var[["variable.name"]]))
				stop("\"Variable Name\" column of ISA MAF data frame and \"variable.name\" column of W4M variable data frame aren't identical.")
	# TODO check variable names in columns (they must be of same length, and ordered the same way)
	# TODO use cbind instead of merge
			maf.df <- cbind(maf.df, w4m.var[cols])
#			maf.df <- merge(maf.df, w4m.var[cols], by.x = "Variable Name", by.y = "variable.name", sort = FALSE)

			# Update MAF data frame
			isa@maf.dataframes[[maf.files[[1]]]] <- maf.df 
		}
	}
}

# Convert W4M 3 files format to ISA {{{1
################################################################

# TODO write documentation
w4m2isa <- function(isa, w4m, study.filename = NULL, assays = NULL, assay.filename = NULL) {

	# The existing ISA object is updated with W4M values.

	# Get study & assays
	study <- .get.study(isa, study.filename = study.filename)
	assay.indices <- if (is.null(assays)) seq(.get.nb.assays(isa, study$name)) else .get.chosen.assay.index(isa, study.name = study$name, assay.filename = assay.filename)
	if ( ! (length(assay.indices) == length(w4m) || (length(assay.indices) == 1 && all(c('samp', 'var', 'mat') %in% names(w4m)))))
		stop(paste0(length(assay.indices), ' assays selected, but only ', length(w4m), ' W4M inputs provided.'))

	# Loop on assay indices
	for (assay.index in (assay.indices)) {

		if (all(c('samp', 'var', 'mat') %in% names(w4m))) {
			w4m.samp <- w4m$samp
			w4m.var <- w4m$var
			w4m.mat <- w4m$mat
		}
		else {
			w4m.samp <- w4m[[assay.index]]$samp
			w4m.var <- w4m[[assay.index]]$var
			w4m.mat <- w4m[[assay.index]]$mat
		}

		# No need to update ISA study data frame. W4M tools only add new columns, and never modify values on existing columns. All new columns of W4M sample data frame will be added to ISA assay data frame.
#		.update.study(isa, study.name = study$name, w4m.samp = w4m.samp)
		# XXX Except when lines have been removed in W4M sample data frame.

		.update.assay(isa, study.name = study$name, assay.index = assay.index, w4m.samp = w4m.samp)
		.update.maf(isa, study.name = study$name, assay.index = assay.index, w4m.var = w4m.var, w4m.mat = w4m.mat)
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
		variable.names <- .make.variable.names(measures$df)

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
