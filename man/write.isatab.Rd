\name{write.isatab}
\alias{write.isatab}
\docType{data}
\title{ Write an ISA-Tab dataset }
\description{
  \code{write.isatab} writes a ISA-Tab dataset.
}
\usage{
write.isatab(isa, path = getwd())
}
\arguments{
  \item{isa}{ the ISA-Tab object as built with the functions \code{isatab2bioc} or \code{isatab2bioczip}.}
  \item{path}{ the directory where the files are going to be written, by default is the working directory}
  }
\author{
  Alejandra Gonzalez-Beltran / 
  Maintainer: Alejandra Gonzalez-Beltran, ISA Team <isatools@googlegroups.com> 
}
\seealso{ \code{write.investigation.file}, \code{write.study.file}, \code{write.assay.file} }
\examples{
    ### Example of writing the study file of faahKO ISA-Tab dataset into the temp director
    faahkoISA <- isatab2bioc(find.package("faahKO"))
    temp = tempdir()
    write.isatab(faahkoISA, temp)
}