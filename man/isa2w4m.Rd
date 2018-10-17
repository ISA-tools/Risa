\name{isa2w4m}
\alias{isa2w4m}
\title{
Convert an ISA object into a W4M object.
}
\description{
Search inside an ISA object for studies and assays to extract and transform into a W4M object (a list with samp, var and mat fields).
}
\usage{isa2w4m(isa, study.filename = NULL, assay.filename = NULL, drop = TRUE, normalize = FALSE)}
\arguments{
	\item{isa}{A ISA object.}
	\item{study.filename}{The filename of the study to extract. If unset, the first study will be extracted.}
	\item{assay.filename}{The filename of the assay to extract. If unset, all assays will be extracted.}
	\item{drop}{If TRUE and only one W4M object is extracted, then simplifies the output and returns only this W4M object instead of returning a list containing this W4M object.}
	\item{normalize}{If TRUE, normalize the variable names and the sample names by using the method make.names(). Be careful that, if you use this option, you will not be able to use the w4m2isa() method, since the names will not be recognized.}
}
\value{
A list of W4M objects, or a single W4M object in case the drop parameter has been set to TRUE and only one W4M object is extracted.
}
\author{
	Pierrick Roger <pk.roger@icloud.com>
}
\examples{
\dontrun{
## Convert an ISA
myw4m <- isa2w4m(myisa)
sample.df <- myw4m$samp # access the sample data frame
variable.df <- myw4m$var # access the variable data frame
matrix.df <- myw4m$mat # access the matrix data frame
}
}
