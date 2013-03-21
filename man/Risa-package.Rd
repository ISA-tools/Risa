\name{Risa}
\alias{Risa}
\alias{Risa-package}
\docType{package}
\title{ Risa package overview }
\description{
  The Risa package builds R objects from ISA-tab format files and it is a component of the ISA software infrastructure \url{http://isa-tools.org/}. The source code and latest version can be found in the GitHub repository \url{https://github.com/ISA-tools/Risa}. Please, submit all 'bugs' and feature requests through \url{https://github.com/ISA-tools/Risa/issues}.
  
  The function \code{readISAtab} reads ISA-tab files from a local zip archive or from a local folder.
  
  It also contains methods to save the ISAtab files and to update the metadata after analysis. For mass spectrometry assays, the package includes several methods to build necessary object for processing using the \code{xcms} package.  
  
  In addition, there is a function called \code{suggestBiocPackage} that, given an object from the ISAtab-class - i.e. an ISAtab dataset - it will retrieve a list of potentially useful packages from Bioconductor, according to the assays measurement and technology types. This functionality is based on the annotations provided by the BiocViews.
  
}
\author{
  Alejandra Gonzalez-Beltran, Steffen Neumann, Audrey Kauffmann, Gabriella Rustici, Philippe Rocca-Serra, Eamonn Maguire, Susanna-Asunta Sansone.
  Maintainer: Alejandra Gonzalez-Beltran \email{isatools@googlegroups.com}  
}
\seealso{ \code{readISAtab}, \code{processAssayXcmsSet}, \code{updateAssayMetadata}, \code{write.ISAtab} }
\examples{
     ## Example of readISAtab for a mass spectrometry experiment from a zip file
     isazip = system.file("extdata","faahKO-metadata.zip", package="Risa")
     isaset = readISAtab(zipfile = isazip, path = file.path(system.file("extdata", package="Risa"),"ISAexample"), verbose =TRUE)

     require("faahKO")
     ## Example of readISAtab for mass spectrometry data
     temp = tempdir()
     datafiles = c(file.path(system.file("cdf/KO", package="faahKO"), grep("CDF",dir(system.file("cdf/KO", package="faahKO")),value=TRUE)), file.path(system.file("cdf/WT", package="faahKO"), grep("CDF",dir(system.file("cdf/WT", package="faahKO")),value=TRUE)))
     file.copy(datafiles, temp, recursive=TRUE)
     isafiles = file.path(system.file(package="faahKO"), grep("txt",dir(system.file(package="faahKO")),value=TRUE))
     file.copy(isafiles, temp, recursive=TRUE)
     isatab.object = readISAtab(path = temp)
}
