# vi: fdm=marker
# Convert W4M 3 files format to ISA {{{1
################################################################

w4m2isa <- function(isa, w4m, study.filename = NULL, assay.filename = NULL) {

	# The existing ISA object is updated with W4M values.

	# Transform w4m into a list of W4M inputs
	if (all(c('samp', 'var', 'mat') %in% names(w4m)))
		w4m <- list(w4m)

	# Get study & assays
	study <- .get.study(isa, study.filename = study.filename)
	assay.indices <- if (is.null(assay.filename)) seq(.get.nb.assays(isa, study$name)) else .get.chosen.assay.index(isa, study.name = study$name, assay.filename = assay.filename)
	if (length(assay.indices) != length(w4m))
		stop(paste0(length(assay.indices), ' assays selected, but only ', length(w4m), ' W4M inputs provided.'))

	# Loop on assay indices
	for (i in seq_along(assay.indices)) {

		assay.index <- assay.indices[[i]]

		# No need to update ISA study data frame. W4M tools only add new columns, and never modify values on existing columns. All new columns of W4M sample data frame will be added to ISA assay data frame.
#		.update.study(isa, study.name = study$name, w4m.samp = w4m.samp)
		# XXX Except when lines have been removed in W4M sample data frame.

		isa <- .update.study(isa, study.name = study$name, w4m = w4m[[i]])
		isa <- .update.maf(isa, study.name = study$name, assay.index = assay.index, w4m = w4m[[i]]) # Update MAF data frame before updating asssay data frame, since it uses the assay data frame.
		isa <- .update.assay(isa, study.name = study$name, assay.index = assay.index, w4m = w4m[[i]])
	}

	return(isa)
}

# Convert ISA to W4M 3 files format {{{1
################################################################

isa2w4m <- function(isa, study.filename = NULL, assay.filename = NULL, drop = TRUE, normalize = FALSE) {

	output <- NULL

	# Get study & assays
	study <- .get.study(isa, study.filename = study.filename)
	assay.indices <- if (is.null(assay.filename)) seq(.get.nb.assays(isa, study$name)) else .get.chosen.assay.index(isa, study.name = study$name, assay.filename = assay.filename)

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
		sample.names <- .get.sample.names(assay.df = assay$df, maf.df = measures$df)

		# Extract sample metadata
		sample.metadata <- .make.sample.metadata(study, assay, sample.names = sample.names, normalize = normalize)

		# Extract variable metadata
		variable.metadata <- .make.variable.metadata(measures = measures, sample.names = sample.names, variable.names = variable.names, normalize = normalize)

		# Extract matrix
		sample.variable.matrix <- .make.matrix(measures = measures, sample.names = sample.names, variable.names = variable.names, normalize = normalize)

		# Build output
		x <- list(samp = sample.metadata, var = variable.metadata, mat = sample.variable.matrix)
		output <- c(output, list(x))
	}

	# Drop
	if (drop && ! is.null(output) && is.null(names(output)) && length(output))
		output <- output[[1]]

	return(output)
}

# Get extract name column {{{1
################################################################

getExtractNameColumn <- function(w4m) {
	"The extract name is like the sample name, but is unique for each extract analysed, even if a sample has been analysed twice."

	samp.name.col <- NULL

	# Get sample names
	sample.names <- colnames(w4m$mat)

	# Loop on all columns of sample data frame
	for (c in colnames(w4m$samp))
		if (all(w4m$samp[[c]] %in% sample.names) && all( ! duplicated(w4m$samp[[c]]))) {
			samp.name.col <- c
			break
		}

	return(samp.name.col)
}

# PRIVATE METHODS {{{1
################################################################

# Get study {{{2
################################################################

.get.study <- function(isa, study.filename = NULL) {

	if ( ! is.null(study.filename)) {
		if ( ! study.filename %in% isa@study.filenames)
			stop(paste("Cannot find study \"", study.filename, "\".", sep = ''))
		study.name <- isa@study.identifiers[[which(study.filename == isa@study.filenames)]]
	} else { # Take the first study
		study.name <- isa@study.identifiers[[1]]
	}
	if (length(study.name) != 1)
		stop("No study found.")

	study.df <- isa@study.files[[study.name]]

	return(list(name = study.name, df = study.df))
}

# Get number of assays {{{2
################################################################

.get.nb.assays <- function(isa, study.name) {
	study.assays <- isa@assay.files.per.study[[study.name]]
	return(length(study.assays))
}

# Get chosen assay index {{{2
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

# Get assay {{{2
################################################################

.get.assay <- function(isa, study.name, assay.index) {

	study.assay.df <- isa@assay.files.per.study[[study.name]][[assay.index]]
	study.assay.filename <- isa@assay.filenames.per.study[[study.name]][[assay.index]]

	return(list(filename = study.assay.filename, df = study.assay.df))
}

# Make sample metadata {{{2
################################################################

.make.sample.metadata <- function(study, assay, sample.names, normalize = TRUE) {

	# Create sample metadata by merging assay and study metadata
	sample.name.col <- 'Sample Name'
	if (normalize) {
		colnames(study$df) <- make.names(colnames(study$df), unique = TRUE) # TODO make.names should be optional
		colnames(assay$df) <- make.names(colnames(assay$df), unique = TRUE)
		sample.name.col <- 'Sample.Name'
	}
	study.df <- merge(assay$df[sample.name.col], study$df, by = sample.name.col, sort = FALSE, no.dups = FALSE)
	colnames(study.df) <- c(sample.name.col, colnames(study$df)[colnames(study$df) != sample.name.col])
	cols <- colnames(assay$df)
	cols.select <- ! cols %in% sample.name.col
	assay.sub.df <- assay$df[cols.select]
	colnames(assay.sub.df) <- cols[cols.select]
	sample.metadata <- cbind(study.df, assay.sub.df)

	# Normalize
	if (normalize) {
		sample.metadata <- cbind(data.frame(sample.name = make.names(sample.names, unique = TRUE), stringsAsFactors = FALSE), sample.metadata)
		colnames(sample.metadata) <- make.names(colnames(sample.metadata), unique = TRUE)
	}

	return(sample.metadata)
}

# Get measures {{{2
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


# Get sample names {{{2
################################################################

.get.sample.names <- function (assay.df, samp.name.col = NULL, maf.df = NULL) {

	sample.names <- NULL

	if ((is.null(samp.name.col) || ! samp.name.col %in% colnames(assay.df)) && ! is.null(maf.df))
		for (sample.field in colnames(assay.df)) {
			if (all(assay.df[[sample.field]] %in% colnames(maf.df)) && all(! duplicated(assay.df[[sample.field]]))) {
				samp.name.col <- sample.field
				break
			}
		}

	if (is.null(samp.name.col))
		stop(paste("Impossible to find a column for sample names. Either such a column does not exist, or it contains duplicates.", sep = ''))

	sample.names <- assay.df[[samp.name.col]]

	return(sample.names)
}

# Make variable metadata {{{2
################################################################

.make.variable.metadata <- function(measures, sample.names, variable.names, normalize = TRUE) {

	variable.metadata <- measures$df[ ! colnames(measures$df) %in% sample.names]

	# Add variable names as columns
	variable.metadata <- cbind(data.frame(variable.name = variable.names, stringsAsFactors = FALSE), variable.metadata)
	colnames(variable.metadata)[[1]] <- 'Variable Name'

	# Normalize
	if (normalize)
		colnames(variable.metadata) <- make.names(colnames(variable.metadata), unique = TRUE)

	return(variable.metadata)
}

# Make matrix {{{2
################################################################

.make.matrix <- function(measures, sample.names, variable.names, normalize = TRUE) {

	if (any( ! sample.names %in% colnames(measures$df)))
		stop(paste("Cannot find sample names in column names of \"", measures$file, "\".", sep = ''))
	sample.variable.matrix <- measures$df[sample.names]

	# Add variable names as columns
	sample.variable.matrix <- cbind(data.frame(variable.name = variable.names, stringsAsFactors = FALSE), sample.variable.matrix)
	colnames(sample.variable.matrix)[[1]] <- 'Variable Name'

	# Normalize sample names
	if (normalize)
		colnames(sample.variable.matrix) <- make.names(colnames(sample.variable.matrix), unique = TRUE)

	return(sample.variable.matrix)
}

# Make variable names {{{2
################################################################

.make.variable.names <- function(maf.df) {

	variable.names <- character(nrow(maf.df))
	for (f in c('mass_to_charge', 'retention_time'))
		if (f %in% colnames(maf.df)) {
			col.vals <- as.character(maf.df[[f]])
			col.vals <- ifelse(is.na(col.vals), '', col.vals)
			variable.names <- paste(variable.names, ifelse(nchar(variable.names) > 0 & nchar(col.vals) > 0, '_', ''), col.vals, sep = '')
		}
	variable.names <- make.names(variable.names, unique = TRUE)

	return(variable.names)
}

# Update assay {{{2
################################################################

.update.assay <- function(isa, study.name, assay.index, w4m) {

	# Get study and assay data frame
	study.df <- isa@study.files[[study.name]]
	study.assay.df <- isa@assay.files.per.study[[study.name]][[assay.index]]

	# Get real sample name column
	real.samp.name.col <- getExtractNameColumn(w4m)
	if ( ! real.samp.name.col %in% colnames(study.assay.df))
		stop(paste0("Cannot find column \"", real.samp.name.col, "\" into ISA assay data frame."))
 	if ( ! real.samp.name.col %in% colnames(w4m$samp))
		stop(paste0("Cannot find column \"", real.samp.name.col, "\" into W4M data frame."))

	# Remove rows when some samples have been deleted
	removed.samples = ! study.assay.df[[real.samp.name.col]] %in% w4m$samp[[real.samp.name.col]]
	if (any(removed.samples))
		study.assay.df <- study.assay.df[ ! removed.samples, ]

	# Add new columns to assay data frame
	cols <- colnames(w4m$samp)
	cols <- cols[ ! cols %in% colnames(study.df)]
	cols <- cols[ ! cols %in% colnames(study.assay.df)]
	if ( ! identical(study.assay.df[[real.samp.name.col]], w4m$samp[[real.samp.name.col]]))
		stop(paste0("\"", real.samp.name.col, "\" column of ISA assay data frame and \"", real.samp.name.col,"\" column of W4M sample data frame aren't identical."))
	study.assay.df <- cbind(study.assay.df, w4m$samp[cols])

	# Update assay data frame
	isa@assay.files.per.study[[study.name]][[assay.index]] <- study.assay.df 

	return(isa)
}

# Update MAF {{{2
################################################################

.update.maf <- function(isa, study.name, assay.index, w4m) {

	assay.filename <- isa@assay.filenames.per.study[[study.name]][[assay.index]]
	# Get MAF data frame
	if (assay.filename %in% names(isa@maf.filenames.per.assay.filename)) {
		maf.files <- isa@maf.filenames.per.assay.filename[[assay.filename]]
		if ( ! is.null(maf.files)) {

			if (length(maf.files) != 1)
				stop(paste("More than one metabolite assignement file found in assay \"", assay.filename, "\": ", paste(maf.files, collapse = ", "), ".", sep = ''))

			# Get MAF data frame
			maf.df <- isa@maf.dataframes[[maf.files[[1]]]]

			# Get assay data frame
			assay.df <- isa@assay.files.per.study[[study.name]][[assay.index]]

			# Get variable names
			maf.df.var.names <- if ('Variable Name' %in% colnames(maf.df)) maf.df[['Variable Name']] else .make.variable.names(maf.df)
 			if ( ! "Variable Name" %in% colnames(w4m$var))
				stop("Cannot find column \"Variable Name\" into W4M variable data frame.")

			# Get real sample name column
			real.samp.name.col <- getExtractNameColumn(w4m)
			if (is.null(real.samp.name.col))
				stop("Cannot find a column that contains all sample names used in MAF data frame.")
			if ( ! real.samp.name.col %in% colnames(assay.df))
				stop(paste0("Cannot find column \"", real.samp.name.col, "\" into ISA assay data frame."))

			# Remove sample columns
			samp.name.to.remove <- assay.df[[real.samp.name.col]][ ! assay.df[[real.samp.name.col]] %in% w4m$samp[[real.samp.name.col ]]]
			if (length(samp.name.to.remove) > 0)
				maf.df <- maf.df[, ! colnames(maf.df) %in% samp.name.to.remove]

			# Remove variable rows
			var.name.to.remove <- maf.df.var.names[ ! maf.df.var.names %in% w4m$var[['Variable Name']]]
			if (length(var.name.to.remove) > 0)
				maf.df <- maf.df[ ! maf.df.var.names %in% var.name.to.remove, ]

			# Add new columns from W4M variable data frame
			cols <- colnames(w4m$var)
			cols <- cols[ ! cols %in% 'Variable Name']
			cols <- cols[ ! cols %in% colnames(maf.df)]

			if (length(cols) > 0) {
				if ( ! identical(maf.df.var.names, w4m$var[["Variable Name"]]))
					stop("Variable names found for ISA MAF data frame and \"Variable Name\" column of W4M variable data frame aren't identical.")

				# Append new columns found in variable data frame
				maf.df <- cbind(maf.df, w4m$var[cols])
			}

			# Update MAF data frame
			isa@maf.dataframes[[maf.files[[1]]]] <- maf.df 
		}
	}

	return(isa)
}

# Update study {{{2
################################################################

.update.study <- function(isa, study.name, w4m) {

	# Get study data frame
	study.df <- isa@study.files[[study.name]]

	# Check that "Sample Name" column exists
	if ( ! "Sample Name" %in% colnames(study.df))
		stop("Cannot find column \"Sample Name\" into ISA assay data frame.")
 	if ( ! "Sample Name" %in% colnames(w4m$samp))
		stop("Cannot find column \"Sample Name\" into W4M sample data frame.")

	# Remove rows when some samples have been deleted
	removed.samples = ! study.df[["Sample Name"]] %in% w4m$samp[["Sample Name"]]
	if (any(removed.samples))
		study.df <- study.df[ ! removed.samples, ]

	# Update study data frame
	isa@study.files[[study.name]] <- study.df

	return(isa)
}

