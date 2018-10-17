\name{w4m2isa}
\alias{w4m2isa}
\docType{data}
\title{
Update an ISA object, by integrating changes made on a W4M object obtained from isa2w4m().
}
\description{
In this method, a W4M object, obtained previously from isa2w4m(), is used to update the ISA object that was given as input to isa2w4m().
The sample, variable and matrix data frames of the W4M object are analysed in order to find added and removed columns, and update accordingly the study, assay and MAF data frames of the ISA object. 
}
\usage{w4m2isa(isa, w4m, study.filename = NULL, assay.filename = NULL)}
\arguments{
	\item{isa}{The ISA object used in isa2w4m().}
	\item{w4m}{The W4M object obtained from isa2w4m().}
	\item{study.filename}{The study filename used in isa2w4m(), if any.}
	\item{assay.filename}{The assay filename used in isa2w4m(), if any.}
}
\value{
An updated ISA object.
}
\author{
	Pierrick Roger <pk.roger@icloud.com>
}
\examples{
\dontrun{
## Convert an ISA
## We suppose the ISA object contains only one study and one assay.
myw4m <- isa2w4m(myisa)

## Add a column
w4m$samp[['My new column']] <- 1234

## Update the ISA
newisa <- w4m2isa(myisa, myw4m)
}
}
