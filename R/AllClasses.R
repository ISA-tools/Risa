### class to hold ISATab information
#' An S4 class that stores ISATab information.
#' ISATab-class
#' 
#' @export
ISATab <- setClass("ISATab",
                   representation(
                     path="character",
                     investigation.filename="character",
                     investigation.file="data.frame",
                     investigation.identifier="character",
                     study.identifiers="character",
                     study.titles="character",
                     study.descriptions="character",
                     study.contacts="character",
                     study.contacts.affiliations="character",
                     study.filenames="character",                 
                     study.files="list",
                     assay.filenames="character",
                     assay.filenames.per.study="list",
                     assay.files="list",
                     assay.files.per.study="list",
                     assay.names="list",
                     assay.technology.types="character",
                     assay.measurement.types="character",
                     data.filenames="list",
                     samples="character",
                     samples.per.study="list",
                     samples.per.assay.filename="list",
                     assay.filenames.per.sample="list",
                     sample.to.rawdatafile="list",
                     sample.to.assayname="list",
                     rawdatafile.to.sample="list",
                     assayname.to.sample="list",
                     factors="list", 
                     treatments="list", 
                     groups="list",
                     assay.tabs="list")) 


AssayTab <- setClass("AssayTab",
                   representation(
                     path="character",
                     study.filename="character",
                     study.identifier="character",
                     assay.filename="character",
                     assay.file="data.frame",
                     assay.technology.type="character",
                     assay.measurement.type="character",
                     assay.names="data.frame",
                     data.filenames="data.frame"
                     ))

validMSAssayTabObject <- function(object) {
  if (object@assay.technology.type=="mass spectrometry") TRUE
  else paste("Technology type is not 'mass spectrometry' for ", 
             object, sep="")
}

MSAssayTab <- setClass("MSAssayTab",
                    representation(),
                    prototype(assay.technology.type="mass spectrometry"),  
                    contains="AssayTab",
                    validity=validMSAssayTabObject)

validMicroarrayAssayTabObject <- function(object) {
  if (object@assay.technology.type==technology.types$microarray) TRUE
  else paste("Technology type is not 'DNA microarray' for ", 
             object, sep="")
}

MicroarrayAssayTab <- setClass("MicroarrayAssayTab",
                       representation(),
                       prototype(assay.technology.type="DNA microarray"),  
                       contains="AssayTab",
                       validity=validMicroarrayAssayTabObject)

SeqAssayTab <- setClass("SeqAssayTab",
                               representation(),
                               prototype(assay.technology.type="nucleotide sequencing"),  
                               contains="AssayTab"
                               )

NMRAssayTab <- setClass("NMRAssayTab",
                        representation(),
                        prototype(assay.technology.type="NMR spectroscopy"),  
                        contains="AssayTab"
)

