\name{addAssayMetadata}
\alias{addAssayMetadata}
\docType{data}
\title{Adds metadata into the assay file}
\description{
  \code{addAssayMetadata} Adds metadata into the assay file.
}
\usage{
addAssayMetadata(isa, assay.filename, col.name, values)
}
\arguments{
  \item{isa}{ An isatab object, as retrieved by the isatab2bioc function.}
  \item{assay.filename}{ the filename of the assay file to be augmented/modified}
  \item{col.name}{ the name of the column of the assay file to be modified}
  \item{values}{ the values to be added to the column of the assay file: it could be a single value, and in this case the value is repeated across the column, or it could be a list of values (whose length must match the number of rows of the assay file) }
  }
\value{
  The updated ISA-Tab dataset
}
\author{
  Alejandra Gonzalez-Beltran / 
  Maintainer: Alejandra Gonzalez-Beltran, ISA Team <isatools@googlegroups.com> 
}
\seealso{ \code{isatab2bioc} }
\examples{
   
}