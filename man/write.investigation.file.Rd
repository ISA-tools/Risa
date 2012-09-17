\name{write.investigation.file}
\alias{write.investigation.file}
\docType{data}
\title{ Writes the Investigation File of an ISA-Tab dataset }
\description{
  \code{write.investigation.file} reads a ISA-tab file.
}
\usage{
write.investigation.file(isa, path = getwd())
}
\arguments{
  \item{isa}{ the ISA-Tab object (as retrieved by the \link{readISAtab} function)}
  \item{path}{the output path in which the study file is going to be written, by default is the working directory}
  }
\author{
  Alejandra Gonzalez-Beltran / 
  Maintainer: Alejandra Gonzalez-Beltran, ISA Team <isatools@googlegroups.com>  
}
\seealso{ \code{write.ISAtab}, \code{write.study.file}, \code{write.assay.file}, \code{readISAtab} }
\examples{
    ### Example of writing the study file of faahKO ISA-Tab dataset into the temp director
    faahkoISA <- readISAtab(find.package("faahKO"))
    temp = tempdir()
    write.investigation.file(faahkoISA, temp)
}