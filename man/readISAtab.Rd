\name{readISAtab}
\alias{readISAtab}
\docType{data}
\title{
Reads an ISAtab dataset into an R object.
}
\description{
Reads an ISAtab dataset given as a zip file or as a set of files in a specific folder, and builds an object from the ISAtab class.
}
\usage{readISAtab(path = getwd(), zipfile = NULL, verbose = FALSE)}
\arguments{
  \item{path}{the name of the directory in which the ISAtab files are located (if the parameter zipfile is not provided or if it is equal to NULL), or the name of the directory where the zip file will be extracted (if the parameter zipfile is not NULL.) The default is the working directory.}
  \item{zipfile}{ a zip archive containing ISAtab files. The default value is NULL.}
  \item{verbose}{a boolean indicating to show messages for the different steps, if TRUE, or not to show them, if FALSE}
}
\value{
The result is an object from the \code{ISAtab} class.
}
\author{
  Alejandra Gonzalez-Beltran, Audrey Kauffmann / 
  Maintainer: Alejandra Gonzalez-Beltran, ISA Team <isatools@googlegroups.com> 
}
\examples{
  ## Example for mass spectrometry dataset
     temp = tempdir()
     datafiles = c(file.path(system.file("cdf/KO", package="faahKO"), grep("CDF",dir(system.file("cdf/KO", package="faahKO")),ignore.case=TRUE,value=TRUE)), file.path(system.file("cdf/WT", package="faahKO"), grep("CDF",dir(system.file("cdf/WT", package="faahKO")),value=TRUE)))
     file.copy(datafiles, temp, recursive=TRUE)
     isafiles = file.path(system.file(package="faahKO"), grep("txt",dir(system.file(package="faahKO")),value=TRUE))
     file.copy(isafiles, temp, recursive=TRUE)
     isaObject1 = readISAtab(path = temp)
     
  ## Example of readISAtab for a mass spectrometry experiment
  isazip <- system.file("extdata","faahKO-metadata.zip", package="Risa")
  isaObject2 <- readISAtab(zipfile = isazip, path = file.path(system.file("extdata", package="Risa"),"ISAexample"), verbose =TRUE)
}
\keyword{ISAtab datasets}