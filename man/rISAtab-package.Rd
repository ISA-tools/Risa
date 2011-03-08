\name{rISAtab}
\alias{rISAtab}
\alias{rISAtab-package}
\docType{package}
\title{ rISAtab package overview }
\description{
  The rISAtab package builds R objects from ISA-tab format files. The package contains two functions:  \code{risazip} reads ISA-tab files from a local zip archive; and \code{risatab} reads ISA-tab files from a local folder.
}
\author{
  Audrey Kauffmann / 
  Maintainer: Audrey Kauffmann, Gabriella Rustici, Philippe Rocca-Serra, <isatools@googlegroups.com>  
}
\seealso{ \code{risazip}, \code{risatab} }
\examples{
     ## Example of risazip for microarray and flow cytometry experiment
     isazip = system.file("extdata","Yox1.zip", package="rISAtab")
     isaset = risazip(isazip, path = "ISAexample")

     ## Example of risatab for mass spectrometry data
     temp = tempdir()
     datafiles = c(file.path(system.file("cdf/KO", package="faahKO"), grep("CDF",dir(system.file("cdf/KO", package="faahKO")),value=TRUE)), file.path(system.file("cdf/WT", package="faahKO"), grep("CDF",dir(system.file("cdf/WT", package="faahKO")),value=TRUE)))
     file.copy(datafiles, temp, recursive=TRUE)
     isafiles = file.path(system.file(package="faahKO"), grep("txt",dir(system.file(package="faahKO")),value=TRUE))
     file.copy(isafiles, temp, recursive=TRUE)
     isaset = risatab(path = temp)
}
