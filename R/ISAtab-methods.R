## getters
#' extract slots of ISAtab class
#'
#' @name [
#' @aliases [,ISAtab-method
#' @docType methods
#' @rdname extract-methods
setMethod(f="[",signature="ISAtab", definition=function(x, i,j, drop) {
  if (i=="path") { return(x@path) } else {}
  if (i=="investigation.filename") { return(x@investigation.filename) } else {}
  if (i=="investigation.file") { return(x@investigation.file) } else {}
  if (i=="study.identifiers") { return(x@study.identifiers) } else {}
  if (i=="study.filenames") { return(x@study.filenames) } else {}
  if (i=="study.files") { return(x@study.files) } else {}   
  if (i=="assay.filenames") { return(x@assay.filenames) } else {}
  if (i=="assay.filenames.per.study") { return(x@assay.filenames.per.study) } else {}
  if (i=="assay.files") { return(x@assay.files) } else {}
  if (i=="assay.files.per.study") { return(x@assay.files.per.study) } else {}
  if (i=="assay.technology.types") { return(x@assay.technology.types) } else {}
  if (i=="assay.measurement.types") { return(x@assay.measurement.types) } else {}
  if (i=="data.filenames") { return(x@data.filenames) } else {}
  if (i=="samples") { return(x@samples) } else {}
  if (i=="samples.per.study") { return(x@samples.per.study) } else {}
  if (i=="samples.per.assay.filename") { return(x@samples.per.assay.filename) } else {}
  if (i=="assay.filenames.per.sample") { return(x@assay.filenames.per.sample) } else {}
  if (i=="sample.to.rawdatafile") { return(x@sample.to.rawdatafile) } else {}
  if (i=="sample.to.assayname") { return(x@sample.to.assayname) } else {}
  if (i=="rawdatafile.to.sample") { return(x@rawdatafile.to.sample) } else {}
  if (i=="assayname.to.sample") { return(x@assayname.to.sample) } else {}
}
) 

## setter methods
#' @name [
#' @aliases [<-,ISAtab-method
#' @docType methods
#' @rdname extract-methods
setReplaceMethod(f="[",signature="ISAtab", definition=function(x,i,j,value){
  if (i=="path") { x@path<-value } else {}
  if (i=="investigation.filename") { x@investigation.filename<-value } else {}
  if (i=="investigation.file") { x@investigation.file<-value } else {}
  if (i=="study.identifiers") { x@study.identifiers<-value } else {}
  if (i=="study.filenames") { x@study.filenames<-value } else {}
  if (i=="study.files") { x@study.files<-value } else {}   
  if (i=="assay.filenames") { x@assay.filenames<-value } else {}
  if (i=="assay.filenames.per.study") { x@assay.filenames.per.study<-value } else {}
  if (i=="assay.files") { x@assay.files<-value } else {}
  if (i=="assay.files.per.study") { x@assay.files.per.study<-value} else {}
  if (i=="assay.technology.types") { x@assay.technology.types<-value} else {}
  if (i=="assay.measurement.types") { x@assay.measurement.types<-value } else {}
  if (i=="data.filenames") { x@data.filenames<-value } else {}
  if (i=="samples") { x@samples<-value } else {}
  if (i=="samples.per.study") { x@samples.per.study<-value } else {}
  if (i=="samples.per.assay.filename") { x@samples.per.assay.filename<-value } else {}
  if (i=="assay.filenames.per.sample") { x@assay.filenames.per.sample<-value } else {}
  if (i=="sample.to.rawdatafile") { x@sample.to.rawdatafile<-value } else {}
  if (i=="sample.to.assayname") { x@sample.to.assayname<-value } else {}
  if (i=="rawdatafile.to.sample") { x@rawdatafile.to.sample<-value } else {}
  if (i=="assayname.to.sample") { x@assayname.to.sample<-value } else {} 
  return (x)
  }
)

setMethod(
  f="initialize",
  signature="ISAtab",
  definition=function(.Object,path){
    
    # Assignment of the slots
    .Object["path"] <- path
    
    #### Parse ISATab files
    d = dir(path)
    
    ## Investigation filename
    ifilename = grep(isatab.syntax$investigation.prefix, d, value=TRUE)
    if (length(ifilename)==0)
      stop("Did not find any investigation file at folder ", path)
    else if (!file.exists(file.path(path, ifilename)))
      stop("Did not find investigation file: ", ifilename)
    
    .Object["investigation.filename"] <- ifilename
    
    ## Reading in investigation file into a data frame
    ifile = read.table(file.path(path, ifilename), sep="\t", fill=TRUE, na.strings = "NA")
    
    .Object["investigation.file"] <- ifile
    
    
    ## Study Identifiers  - as a list of strings
    sidentifiers = ifile[grep(isatab.syntax$study.identifier, ifile[,1], useBytes=TRUE),][2][[1]]
    
    .Object["study.identifiers"] <- sidentifiers
    
    ## Study filenames (one or more)
    sfilenames = unlist(sapply(ifile[grep(isatab.syntax$study.file.name, ifile[,1], useBytes=TRUE),], function(i) grep(isatab.syntax$study.prefix, i, value=TRUE, useBytes=TRUE)))
    if (length(sidentifiers)!=length(sfilenames))
      stop("There are study files with no identifier assigned")
    ## Assign sidentifiers as names of the list sfilenames
    names(sfilenames) <- sidentifiers
    
    .Object["study.filenames"] <- sfilenames
    
    ## TODO pretty printing sfilenames
    ## Validation of existance of study files
    if (!all(sapply(sfilenames, function(i) file.exists(file.path(path, i)))))
      stop("Did not find some of the study files: ", sfilenames)
    
    ## Reading study files into a list of data frames
    sfiles = lapply(sfilenames, function(i) read.table(file.path(path, i), sep="\t", header=TRUE, stringsAsFactors=FALSE, check.names=FALSE))
    
    .Object["study.files"] <- sfiles
    
    ## List of assay filenames 
    #afilenames is a list with all the assay filenames (without association to studies)
    afilenames = unlist(sapply(ifile[grep(isatab.syntax$study.assay.file.name, ifile[,1], useBytes=TRUE),], function(i) grep(isatab.syntax$assay.prefix, i, value=TRUE, useBytes=TRUE)))
    
    .Object["assay.filenames"] <- afilenames
    
    
    #getting afilenames associated with studies
    afilenames.df = ifile[grep(isatab.syntax$study.assay.file.name, ifile[,1], useBytes=TRUE),]
    afilenames.matrix = apply(afilenames.df,c(1,2),function(row) grep(isatab.syntax$assay.prefix,row, value=TRUE))  
    afilenames.lists = split(afilenames.matrix, row(afilenames.matrix, as.factor=TRUE))
    afilenames.per.study = lapply(seq_len(length(afilenames.lists)), function(i) Filter(function(j) !identical(character(0), j), afilenames.lists[[i]]))
    names(afilenames.per.study) <- sidentifiers
    
    .Object["assay.filenames.per.study"] <- afilenames.per.study
    
    
    ## Reading in assay files 
    # afiles is a list of data frames (containing all the assay files)
    afiles <- lapply(afilenames, function(i) read.table(file.path(path, i), sep="\t", header=TRUE, stringsAsFactors=FALSE,  check.names=FALSE))
    names(afiles) <- afilenames
    
    .Object["assay.files"] <- afiles
    
    # afiles.per.study is a list (one element per study) of lists (one element per assay) 
    afiles.per.study = lapply(seq_len(length(afilenames.per.study)), 
                              function(j) (lapply(seq_len(length(afilenames.per.study[[j]])),
                                                  function(i) read.table(file.path(path,afilenames.per.study[[j]][[i]]), sep="\t", header=TRUE, stringsAsFactors=FALSE, check.names=FALSE))))
    names(afiles.per.study) <- sidentifiers
    
    .Object["assay.files.per.study"] < afiles.per.study
    
    
    ## Assay technology types
    #data frame with types
    assay.tech.types = ifile[which(ifile[[1]]==isatab.syntax$study.assay.technology.type),] 
    #remove empty types - results in a list of types
    assay.tech.types = na.omit(assay.tech.types[assay.tech.types != ""])
    #remove headers
    assay.tech.types = assay.tech.types[ assay.tech.types != isatab.syntax$study.assay.technology.type]
    
    ## Validate number of assay technology types == number of afiles
    if (length(assay.tech.types)!=length(afiles)){
      stop("The number of assay files mismatches the number of assay types")
    }
    
    .Object["assay.technology.types"] <- assay.tech.types
    
    
    ## Assay measurement types
    assay.meas.types = ifile[which(ifile[[1]]==isatab.syntax$study.assay.measurement.type),] 
    assay.meas.types = na.omit(assay.meas.types[assay.meas.types != ""])
    assay.meas.types = assay.meas.types[ assay.meas.types != isatab.syntax$study.assay.measurement.type]
    
    .Object["assay.measurement.types"] <- assay.meas.types
    
    
    ## Identifying what sample is studied in which assay
    ## assays is a list of data frames (one for each assay file)
    assays = lapply(seq_len(length(sfiles)), 
                    function(j) (lapply(seq_len(length(afiles)), 
                                        function(i) sfiles[[j]]$Sample.Name %in% afiles[[i]]$Sample.Name)))
    
    
    samples = unlist(lapply(sfiles, function(i) i[,grep(isatab.syntax$sample.name, colnames(i))]))
    
    .Object["samples"] <- samples
    
    samples.per.assay.filename = lapply(seq_len(length(afiles)), 
                                        function(i) afiles[[i]][[isatab.syntax$sample.name]])
    names(samples.per.assay.filename) <- afilenames
    
    .Object["samples.per.assay.filename"] <- samples.per.assay.filename
    
    
    samples.per.study <- lapply(seq_len(length(sfiles)),
                                function(i) sfiles[[i]][[isatab.syntax$sample.name]])
    names(samples.per.study) <- sidentifiers
    
    .Object["samples.per.study"] <- samples.per.study
    
    assay.filenames.per.sample <- unlist(lapply(seq_len(length(samples)), 
                                                function(j) lapply(seq_len(length(afilenames)), 
                                                                   function(i)   if (samples[[j]] %in% afiles[[i]][[isatab.syntax$sample.name]]) {
                                                                     afilenames[[i]]
                                                                   }
                                                )))
    
    .Object["assay.filenames.per.sample"] <- assay.filenames.per.sample
   
    .Object <- setAssayDependentSlots(.Object)
    
    return(.Object) # return of the object
    }
  )

setGeneric("setAssayFile",function(.Object,assay.filename,assay.file){standardGeneric("setAssayFile")})
setMethod("setAssayFile",
          signature(.Object = "ISAtab", assay.filename = "character", assay.file = "data.frame"),
          function (.Object, assay.filename, assay.file) 
          {
            .Object["assay.files"][[assay.filename]] <- assay.file
            .Object <- setAssayDependentSlots(.Object)
            return(.Object)
          }
)

setGeneric("setAssayDependentSlots",function(.Object){standardGeneric("setAssayDependentSlots")})
setMethod("setAssayDependentSlots",
          signature(.Object = "ISAtab"),
          function (.Object) 
          { 
            afiles <- .Object["assay.files"]
            
            ## List of data filenames with assay filenames as keys
            dfilenames.per.assay = lapply(afiles, function(i) i[,grep(isatab.syntax$data.file, colnames(i))])
            
            .Object["data.filenames"] <- dfilenames.per.assay
            
            
            data.col.names = lapply(seq_len(length(afiles)),
                                    function(i) if (isatab.syntax$raw.data.file %in% colnames(afiles[[i]])){
                                      isatab.syntax$raw.data.file
                                    }else if (isatab.syntax$free.induction.decay.data.file %in% colnames(afiles[[i]])){
                                      isatab.syntax$free.induction.decay.data.file
                                    }else if (isatab.syntax$array.data.file %in% colnames(afiles[[i]])){
                                      isatab.syntax$array.data.file
                                    }else if (isatab.syntax$raw.spectral.data.file %in% colnames(afiles[[i]])){
                                      isatab.syntax$raw.spectral.data.file
                                    })
            
            
            sample.to.rawdatafile <- lapply( seq_len(length(afiles)), 
                                             function(i) afiles[[i]][,c(isatab.syntax$sample.name,data.col.names[[i]])] )
            sample.to.rawdatafile <- lapply(seq_len(length(afiles)), function(i)
              merge(sample.to.rawdatafile[[i]][ !duplicated(sample.to.rawdatafile[[i]] [[isatab.syntax$sample.name]]), ], sample.to.rawdatafile[[i]][ duplicated(sample.to.rawdatafile[[i]][[isatab.syntax$sample.name]]), ], all=TRUE))  
            
            .Object@sample.to.rawdatafile <- sample.to.rawdatafile
            
            sample.to.assayname <-lapply( afiles,
                                          function(i) i[,c(isatab.syntax$sample.name,grep(isatab.syntax$assay.name, colnames(i), value=TRUE))])
            sample.to.assayname <- lapply(seq_len(length(afiles)), function(i)
              merge(sample.to.assayname[[i]][ !duplicated(sample.to.assayname[[i]][[isatab.syntax$sample.name]]), ], sample.to.assayname[[i]][ duplicated(sample.to.assayname[[i]][[isatab.syntax$sample.name]]), ], all=TRUE))
            
            .Object["sample.to.assayname"] <- sample.to.assayname
            
            rawdatafile.to.sample <- lapply( seq_len(length(afiles)), 
                                             function(i) afiles[[i]][,c(data.col.names[[i]],isatab.syntax$sample.name)] )
            rawdatafile.to.sample <- lapply(seq_len(length(afiles)), function(i)
              merge(rawdatafile.to.sample[[i]][ !duplicated(rawdatafile.to.sample[[i]][[data.col.names[[i]]]]), ], rawdatafile.to.sample[[i]][ duplicated(rawdatafile.to.sample[[i]][[data.col.names[[i]]]]), ], all=TRUE))
            
            .Object@rawdatafile.to.sample <- rawdatafile.to.sample
            
            assayname.to.sample <- lapply( afiles,
                                           function(i) i[,c(grep(isatab.syntax$assay.name, colnames(i), value=TRUE),isatab.syntax$sample.name)])
            assayname.to.sample <- lapply(seq_len(length(afiles)), function(i)
              merge(assayname.to.sample[[i]][ !duplicated(assayname.to.sample[[i]][,c(grep(isatab.syntax$assay.name, colnames(assayname.to.sample[[i]]), value=TRUE))]), ], 
                    assayname.to.sample[[i]][  duplicated(assayname.to.sample[[i]][,c(grep(isatab.syntax$assay.name, colnames(assayname.to.sample[[i]]), value=TRUE))]), ], 
                    all=TRUE))
            
            .Object["assayname.to.sample"] <- assayname.to.sample
            return(.Object)
          }
)


