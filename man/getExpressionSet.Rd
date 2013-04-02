\name{getExpressionSet}
\alias{getExpressionSet}
\title{
Retrieves an ExpressionSet for a DNA microarray assay.
}
\description{
Retrieves an ExpressionSet for a DNA microarray assay.
}
\usage{
getExpressionSet(isa, assay.filename)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{isa}{
An object from the \code{ISATab-class}.
}
  \item{assay.filename}{
A character vector with the assay filename of the microarray assay.
}
}
%\details{
%}
\value{
It returns an \code{ExpressionSet}.
}
\references{
\code{affy}
}
\author{
Alejandra Gonzalez-Beltran
}
%\note{
%}
\seealso{
\code{\link{ISATab-class}}
}
\examples{
data.dir <- system.file("extdata", package="Risa")
isazip <- system.file("extdata","ARMSTRONG-S-3-metadata.zip", package="Risa")
path <- file.path(data.dir,"ARMSTRONG_ISA")
isa <- readISAtab(zipfile = isazip, path = path, verbose =TRUE)
path <- file.path(data.dir,"ARMSTRONG_ISA/")

cel.files <- c("GSM510377_SA2006121413.CEL", "GSM510378_SA2006121414.CEL", "GSM510379_SA2006121415.CEL", "GSM510380_SA2006121416.CEL",
                "GSM510381_SA2006121417.CEL", "GSM510382_SA2006121418.CEL","GSM510383_SA2006121419.CEL","GSM510384_SA2006121420.CEL",
                "GSM510385_SA2006121421.CEL","GSM510386_SA2006121427.CEL", "GSM510387_SA2006121428.CEL", "GSM510388_SA2006121429.CEL",
                "GSM510389_SA2006121440.CEL","GSM510390_SA2006121441.CEL", "GSM510391_SA2006121442.CEL")

base.url <- "http://perdera.sph.harvard.edu/biidata/microarray/study_ARMSTRONG-S-3_3Zd4KeJNJ7/raw_data/"

###########download all the cel files - this is commented out, uncomment the following two lines if you want to build the expression set 
#lapply(cel.files, function(x) download.file(paste(base.url, x, sep=""), destfile=paste(path, x, sep="")))
#getExpressionSet(isa, isa@assay.filenames[[1]])
                
}
\keyword{microarray}
\keyword{affy}
