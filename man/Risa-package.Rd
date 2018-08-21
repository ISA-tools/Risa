\name{Risa}
\alias{Risa}
\alias{Risa-package}
\docType{package}
\title{ Risa package overview }
\description{
  The Risa package builds R objects from ISA-tab format files and it is a component of the ISA software infrastructure (see \url{http://isa-tools.org/}). ISA stands for Investigation, Study, Assay and the infrastructure revolves around the ISA-Tab format, which supports descriptions of biological experiments, with focus but not limited to high-throughput experiments. The format is supported by a series of open-source tools (all available in There is a growing international user community of the ISA infrastructure, featured in the ISA commons website (\url{http://isacommons.org}).
  
  The source code and latest version of the Risa package can be found in the GitHub repository \url{https://github.com/ISA-tools/Risa}. Please, submit all 'bugs' and feature requests through \url{https://github.com/ISA-tools/Risa/issues}.
  
  The function \code{readISAtab} reads ISA-tab files from a local zip archive or from a local folder.
  
  The package also contains methods to save the ISA-tab files and to update the metadata after analysis. For mass spectrometry assays, the package includes several methods to build necessary object for processing using the \code{xcms} package (available in Bioconductor).  
  
  There are also functions to deal with DNA microarray-based assays, such as building an \code{ExpressionSet} object, relying on the \code{affy} package (available in Bioconductor)
  
  In addition, there is a function called \code{suggestBiocPackage} that, given an object from the ISAtab-class - i.e. an ISAtab dataset - it will retrieve a list of potentially useful packages from Bioconductor, according to the assays measurement and technology types. This functionality is based on the annotations provided by the BiocViews.
  
  The Risa package bridges the ISA metadata with analysis pipelines in R, and this is the reason of its name: R + ISA = Risa. Anecdotically, we note that Risa is the Pleasure Planet in Star Trek (see \url{http://en.memory-alpha.org/wiki/Risa}) and that 'Risa' means laughter in Spanish.
  
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
