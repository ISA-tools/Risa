### class to hold ISATab information
#' An S4 class that stores ISAtab information.
#' ISAtab-class
#' 
#' @export
ISAtab <- setClass("ISAtab",
                   representation(
                     path="character",
                     investigation.filename="character",
                     investigation.file="data.frame",
                     study.identifiers="factor",
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
                     assay.filenames.per.sample="character",
                     sample.to.rawdatafile="list",
                     sample.to.assayname="list",
                     rawdatafile.to.sample="list",
                     assayname.to.sample="list",
                     factors="list", 
                     treatments="list", 
                     groups="list")) 
