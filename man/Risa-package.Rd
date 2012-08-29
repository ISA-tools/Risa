\name{Risa}
\alias{Risa}
\alias{Risa-package}
\docType{package}
\title{ Risa package overview }
\description{
  The Risa package builds R objects from ISA-tab format files and it is a component of the ISA software infrastructure \url{http://isa-tools.org/}. The source code and latest version can be found in the GitHub repository \url{https://github.com/ISA-tools/Risa}. Please, submit all 'bugs' and feature requests through \url{https://github.com/ISA-tools/Risa/issues}.
  
  The Risa package contains two functions:\code{isatab2bioczip} reads ISA-tab files from a local zip archive; and \code{isatab2bioc} reads ISA-tab files from a local folder.
}
\author{
  Alejandra Gonzalez-Beltran, Steffen Neumann, Audrey Kauffmann, Gabriella Rustici, Philippe Rocca-Serra, Eamonn Maguire, Susanna-Asunta Sansone.
  Maintainer: Alejandra Gonzalez-Beltran \email{isatools@googlegroups.com}  
}
\seealso{ \code{isatab2bioczip}, \code{isatab2bioc}, \code{processAssayXcmsSet}, \code{addAssayMetadata}, \code{write.isatab} }
\examples{
     ## Example of isatab2bioczip for a mass spectrometry experiment
     isazip = system.file("extdata","faahKO-metadata.zip", package="Risa")
     isaset = isatab2bioczip(isazip, path = file.path(system.file("extdata", package="Risa"),"ISAexample"), verbose =TRUE)

     require("faahKO")
     ## Example of isatab2bioc for mass spectrometry data
     temp = tempdir()
     datafiles = c(file.path(system.file("cdf/KO", package="faahKO"), grep("CDF",dir(system.file("cdf/KO", package="faahKO")),value=TRUE)), file.path(system.file("cdf/WT", package="faahKO"), grep("CDF",dir(system.file("cdf/WT", package="faahKO")),value=TRUE)))
     file.copy(datafiles, temp, recursive=TRUE)
     isafiles = file.path(system.file(package="faahKO"), grep("txt",dir(system.file(package="faahKO")),value=TRUE))
     file.copy(isafiles, temp, recursive=TRUE)
     isatab.object = isatab2bioc(path = temp)
}
