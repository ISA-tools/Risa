\name{isatab2bioc}
\alias{isatab2bioc}
\docType{data}
\title{ Build R objects from ISA-tab format files }
\description{
  \code{isatab2bioc} reads a ISA-tab file.
}
\usage{
isatab2bioc(path = getwd())
}
\arguments{
  \item{path}{ the name of the directory containing ISAtab files. The
    default is the working directory.}
  }
\value{
  The output is a list composed, from a data.frame containing the metadata information originally contained in the s_ file of ISAtab, followed by an object per assay of the study which will be an \code{\link[affy:AffyBatch-class]{AffyBatch}} if the assay is a raw Affymetrix microarray dataset, an \code{\link[Biobase:class.ExpressionSet]{ExpressionSet}} for a raw one colour microarray dataset or for a processed microarray dataset, a \code{\link[Biobase:class.NChannelSet]{NChannelSet}} from a raw two colour dataset, a \code{\link[flowCore:flowSet-class]{flowSet}} for a flow cytometry dataset and a \code{\link[xcms:xcmsSet-class]{xcmsSet}} for a mass spectrometry dataset.
}
\author{
  Audrey Kauffmann / 
  Maintainer: Audrey Kauffmann, Gabriella Rustici, Philippe Rocca-Serra, <isatools@googlegroups.com>  
}
\seealso{ \code{isatab2bioczip} }

\examples{
     ## Example for mass spectrometry
     temp = tempdir()
     datafiles = c(file.path(system.file("cdf/KO", package="faahKO"), grep("CDF",dir(system.file("cdf/KO", package="faahKO")),ignore.case=TRUE,value=TRUE)), file.path(system.file("cdf/WT", package="faahKO"), grep("CDF",dir(system.file("cdf/WT", package="faahKO")),value=TRUE)))
     file.copy(datafiles, temp, recursive=TRUE)
     isafiles = file.path(system.file(package="faahKO"), grep("txt",dir(system.file(package="faahKO")),value=TRUE))
     file.copy(isafiles, temp, recursive=TRUE)
     isaset = isatab2bioc(path = temp)
}
