\name{updateAssayMetadata}
\alias{updateAssayMetadata}
\docType{data}
\title{Update metadata into the assay file}
\description{
  \code{updateAssayMetadata} Updates metadata into a particular assay file. See an example of use for a metabolite dataset at \url{https://github.com/sneumann/mtbls2}.
}
\usage{
updateAssayMetadata(isa, assay.filename, col.name, values)
}
\arguments{
  \item{isa}{ An isatab object, as retrieved by the \code{readISAtab} function.}
  \item{assay.filename}{ the filename of the assay file to be augmented/modified}
  \item{col.name}{ the name of the column of the assay file to be modified}
  \item{values}{ the values to be added to the column of the assay file: it could be a single value, and in this case the value is repeated across the column, or it could be a list of values (whose length must match the number of rows of the assay file) }
  }
\value{
  The updated ISA-Tab object.
}
\author{
  Alejandra Gonzalez-Beltran / 
  Maintainer: Alejandra Gonzalez-Beltran, ISA Team <isatools@googlegroups.com> 
}
\seealso{ \code{readISAtab}, \url{https://github.com/sneumann/mtbls2}}
\examples{
 ### This example shows how to add values to the column "Derived Spectral Data File", assuming that the results are stored in the file "faahkoDSDF.txt"
 faahkoISA = readISAtab(find.package("faahKO"))
 assay.filename <- faahkoISA["assay.filenames"][[1]]
 updateAssayMetadata(faahkoISA, assay.filename,"Derived Spectral Data File","faahkoDSDF.txt" )
}