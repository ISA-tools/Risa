\name{getExtractNameColumn}
\alias{getExtractNameColumn}
\docType{data}
\title{
Get the name of the column containing the extract names used as column names inside the MAF data frame and inside the W4M data matrix.
}
\description{
Using a W4M object, finds the column of the sample data frame that contains all extract/sample names used inside the data matrix.
This is the same names that are used inside the ISA MAF data frame.
}
\usage{getExtractNameColumn(w4m)}
\arguments{
	\item{w4m}{A W4M object: a list with three fields (samp, var and mat), each containing the corresponding data frame.}
}
\value{
The name of the column containing the extract names, in the sample data frame.
}
\author{
	Pierrick Roger <pk.roger@icloud.com>
}
\examples{
\dontrun{
## Use after a conversion
myw4m <- isa2w4m(myisa)
extract.name.col <- getExtractNameColumn(myw4m)
}
}
