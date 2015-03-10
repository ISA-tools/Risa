## getters
#' extract slots of ISATab class
#'
#' @name [
#' @aliases [,ISATab-method
#' @docType methods
#' @rdname extract-methods
setMethod(f="[",signature="ISATab", definition=function(x, i,j, drop) {
  if (i=="path") { return(x@path) } else {}
  if (i=="investigation.filename") { return(x@investigation.filename) } else {}
  if (i=="investigation.file") { return(x@investigation.file) } else {}
  if (i=="investigation.identifier") { return(x@investigation.identifier) } else {}
  if (i=="study.identifiers") { return(x@study.identifiers) } else {}
  if (i=="study.titles") { return(x@study.titles) } else {}
  if (i=="study.descriptions") { return(x@study.descriptions) } else {}
  if (i=="study.contacts") { return(x@study.contacts) } else {}
  if (i=="study.contacts.affiliations") { return(x@study.contacts.affiliations) } else {}
  if (i=="study.filenames") { return(x@study.filenames) } else {}
  if (i=="study.files") { return(x@study.files) } else {}   
  if (i=="assay.filenames") { return(x@assay.filenames) } else {}
  if (i=="assay.filenames.per.study") { return(x@assay.filenames.per.study) } else {}
  if (i=="assay.files") { return(x@assay.files) } else {}
  if (i=="assay.names") { return(x@assay.names) } else {}
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
  if (i=="factors") { return(x@factors) } else {}
  if (i=="treatments") { return(x@treatments) } else {}
  if (i=="groups") { return(x@groups) } else {}
  if (i=="assay.tabs") { return(x@assay.tabs) } else {}
}
) 

## setter methods
#' @name [
#' @aliases [<-,ISATab-method
#' @docType methods
#' @rdname extract-methods
setReplaceMethod(f="[",signature="ISATab", definition=function(x,i,j,value){
  if (i=="path") { x@path<-value } else {}
  if (i=="investigation.filename") { x@investigation.filename<-value } else {}
  if (i=="investigation.file") { x@investigation.file<-value } else {}
  if (i=="investigation.identifier") { x@investigation.identifier<-value } else {}
  if (i=="study.identifiers") { x@study.identifiers<-value } else {}
  if (i=="study.titles") { x@study.titles<-value } else {}
  if (i=="study.descriptions") { x@study.descriptions<-value } else {}  
  if (i=="study.contacts") { x@study.contacts<-value } else {}
  if (i=="study.contacts.affiliations") { x@study.contacts.affiliations<-value } else {}
  if (i=="study.filenames") { x@study.filenames<-value } else {}
  if (i=="study.files") { x@study.files<-value } else {}   
  if (i=="assay.filenames") { x@assay.filenames<-value } else {}
  if (i=="assay.filenames.per.study") { x@assay.filenames.per.study<-value } else {}
  if (i=="assay.files") { x@assay.files<-value } else {}
  if (i=="assay.names") { x@assay.names<-value } else {}
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
  if (i=="factors") { x@factors<-value } else {} 
  if (i=="treatments") { x@treatments<-value } else {} 
  if (i=="groups") { x@groups<-value } else {} 
  if (i=="assay.tabs") { x@assay.tabs<-value } else {} 
  return (x)
  }
)


## getters
#' extract slots of AssayTab class
#'
#' @name [
#' @aliases [,AssayTab-method
#' @docType methods
#' @rdname extract-methods
setMethod(f="[",signature="AssayTab", definition=function(x, i,j, drop) {
  if (i=="path") { return(x@path) } else {}
  if (i=="data.filenames") { return(x@data.filenames) } else {}
}
) 

## setter methods
#' @name [
#' @aliases [<-,AssayTab-method
#' @docType methods
#' @rdname extract-methods
setReplaceMethod(f="[",signature="AssayTab", definition=function(x,i,j,value){
  if (i=="path") { x@path<-value } else {}
  if (i=="data.filenames") { x@data.filenames<-value } else {}
}
)


setMethod(
  f="initialize",
  signature="ISATab",
  definition=function(.Object,path){
    
    # Assignment of the slots
    .Object["path"] <- path
    
    #### Parse ISATab files
    d = dir(path)
    
    ## Investigation filename, and avoid editor autosave-files
    ifilename = grep(paste("^", Risa:::isatab.syntax$investigation.prefix, ".*[a-zA-Z]$", sep=""), d, value=TRUE, perl=TRUE)
    if (length(ifilename)==0)
      stop("Did not find any investigation file at folder ", path)
    else if (length(ifilename)>1)
      stop("Found too many possible investigation files: ", ifilename)
    else if (!file.exists(file.path(path, ifilename)))
      stop("Did not find investigation file: ", ifilename)
          
    .Object["investigation.filename"] <- ifilename
    
    ## Reading in investigation file into a data frame
    number.columns <- max(count.fields(file.path(path, ifilename), sep = "\t", comment.char = "#", quote = "\"'", blank.lines.skip = TRUE), na.rm = TRUE)
    
    ifile = read.table(file.path(path, ifilename), sep="\t", fill=TRUE, na.strings = "NA", comment.char = "#", blank.lines.skip = TRUE , col.names = paste0("V",seq_len(number.columns)))
    
    .Object["investigation.file"] <- ifile
    
    iidentifier <- as.character(ifile[grep(isatab.syntax$investigation.identifier, ifile[,1], useBytes=TRUE),][2][[1]])
    
    .Object["investigation.identifier"] <- iidentifier
        
    ## Study Identifiers  - as a list of strings
    sidentifiers = ifile[grep(isatab.syntax$study.identifier, ifile[,1], useBytes=TRUE),][2][[1]]
    
    .Object["study.identifiers"] <- as.character(sidentifiers)
    
    stitles = ifile[grep(isatab.syntax$study.title, ifile[,1], useBytes=TRUE),][2][[1]]
    
    .Object["study.titles"] <- as.character(stitles)
    
    sdescriptions = ifile[grep(isatab.syntax$study.description, ifile[,1], useBytes=TRUE),][2][[1]]
    
    .Object["study.descriptions"] <- as.character(sdescriptions)
    
    spersonfirstnames <- trim(as.character(ifile[grep(isatab.syntax$study.person.first.name, ifile[,1], useBytes=TRUE),][2][[1]]))
    spersonlastnames <- trim(as.character(ifile[grep(isatab.syntax$study.person.last.name, ifile[,1], useBytes=TRUE),][2][[1]]))
    spersonmidinitials <- trim(as.character(ifile[grep(isatab.syntax$study.person.mid.initial, ifile[,1], useBytes=TRUE),][2][[1]]))
    
    .Object["study.contacts"] <- apply(cbind(spersonfirstnames, spersonmidinitials, spersonlastnames), 1, paste, collapse=" " )
   
    spersonaffiliations <- trim(as.character(ifile[grep(isatab.syntax$study.person.affiliation, ifile[,1], useBytes=TRUE),][2][[1]]))
    
    .Object["study.contacts.affiliations"] <- spersonaffiliations
    
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
    
    .Object["assay.files.per.study"] <- afiles.per.study
    
  
    ###assay.names
    assay.names <- lapply( afiles, function(i) i[ grep(isatab.syntax$assay.name, colnames(i) ) ])  
    .Object["assay.names"] <- assay.names
        
    ## Assay technology types
    #data frame with types
    assay.tech.types = ifile[which(ifile[[1]]==isatab.syntax$study.assay.technology.type),] 
    #remove empty types - results in a list of types
    assay.tech.types = na.omit(assay.tech.types[assay.tech.types != ""])
    #remove headers
    assay.tech.types = assay.tech.types[ assay.tech.types != isatab.syntax$study.assay.technology.type]
    
    ## Validate number of assay technology types == number of afiles
    if (length(assay.tech.types)!=length(afiles)){
      message("The number of assay files mismatches the number of assay types")
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
        
    samples = unique(unlist(lapply(sfiles, function(i) i[,grep(isatab.syntax$sample.name, colnames(i))])))
    
    .Object["samples"] <- samples
    
    samples.per.assay.filename = lapply(seq_len(length(afiles)), 
                                        function(i) afiles[[i]][[isatab.syntax$sample.name]])
    names(samples.per.assay.filename) <- afilenames
    
    .Object["samples.per.assay.filename"] <- samples.per.assay.filename
        
    samples.per.study <- lapply(seq_len(length(sfiles)),
                                function(i) sfiles[[i]][[isatab.syntax$sample.name]])
    names(samples.per.study) <- sidentifiers
    
    .Object["samples.per.study"] <- samples.per.study
      
    assay.filenames.per.sample <- list()
    for (k in seq_len(length(samples))){
      for (j in seq_len(length(afilenames))){
        if (samples[[k]] %in% afiles[[j]][[isatab.syntax$sample.name]]) {
          if (length(assay.filenames.per.sample) < k)
            assay.filenames.per.sample[[k]] <- list()
          assay.filenames.per.sample[[k]] <- c(assay.filenames.per.sample[[k]],afilenames[[j]])
        }
      }     
    }                                                                                                                                                                                                                               
    
    
    if (is.null(assay.filenames.per.sample)){
     message("assay.filenames.per.sample not assigned") 
    }else{      
       .Object["assay.filenames.per.sample"] <- assay.filenames.per.sample
    }
        
      
    .Object <- setAssayDependentSlots(.Object)
    
    .Object <- setFactors(.Object)
        
    .Object <- setTreatments(.Object)
    
    .Object <- setGroups(.Object)
    
    .Object <- setAssayTabs(.Object)
    
    return(.Object) 
    }
  )



setGeneric("setAssayFile",function(.Object,assay.filename,assay.file){standardGeneric("setAssayFile")})
setMethod("setAssayFile",
          signature(.Object = "ISATab", assay.filename = "character", assay.file = "data.frame"),
          function (.Object, assay.filename, assay.file) 
          {
            .Object["assay.files"][[assay.filename]] <- assay.file
            .Object <- setAssayDependentSlots(.Object)
            return(.Object)
          }
)

setGeneric("setAssayDependentSlots",function(.Object){standardGeneric("setAssayDependentSlots")})
setMethod("setAssayDependentSlots",
          signature(.Object = "ISATab"),
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

setGeneric("setFactors",function(.Object){standardGeneric("setFactors")})
setMethod("setFactors",
          signature(.Object = "ISATab"),
          function (.Object) 
          { 
            study.files <- .Object["study.files"]  
            factors.list <- list()
            for(i in seq_len(length(study.files))){
              if (length(grep("Factor.Value", colnames(study.files[[i]]))) != 0) {
                factor.values  <-  study.files[[i]][ grep("Factor.Value", colnames(study.files[[i]]))]
                factors.list[[i]] <- lapply(factor.values, factor)  
              }else{
                study.filenames <- .Object["study.filenames"]
                factors.list[[i]] <- list()                
                message("No 'Factor Value' column defined in study file ",study.filenames[[i]], ". Factors slot will be an empty list for that study.")
              }
            }           
            .Object["factors"] <- factors.list
            return(.Object)
          }
          )

setGeneric("setTreatments",function(.Object){standardGeneric("setTreatments")})
setMethod("setTreatments",
          signature(.Object = "ISATab"),
          function (.Object) 
          { 
            factors <- .Object["factors"] 
            if (length(factors) == 0){
              treatments <- list()      
              message("Treatments slot will be an empty list")
            } else {
              factors.df.list <- lapply(factors, as.data.frame)
              for (i in seq(factors.df.list)){
                colnames(factors.df.list[[i]]) <-  names(factors[[i]])
              } 
           
              treatments <- lapply(factors.df.list, function(factors.df) factors.df[!duplicated(factors.df),])           
              names(treatments) <- sapply(factors.df.list, function(x) paste(colnames(x), collapse=' '))
            }
            
            .Object["treatments"] <- treatments            
            return(.Object)
          }
)


setGeneric("setGroups",function(.Object){standardGeneric("setGroups")})
setMethod("setGroups",
          signature(.Object = "ISATab"),
          function (.Object) 
          {                     
            treatments <- .Object["treatments"]            
                        
            study.files <- .Object["study.files"]
            
            samples.per.study <- .Object["samples.per.study"]
                                    
            groups <- list()
            
            if (length(treatments) != 0){
            
            for (j in seq_len(length(study.files))){
              subgroups <- list()
             
              if (class(treatments[[j]]) == "factor"){
                               
                for(i in seq_len(length(levels(treatments[[j]])))){
                  treatment <- treatments[[j]][[i]]  
                  list <-  rep(treatment, each = length(samples.per.study[[j]]))
                  subgroups[[i]] = samples.per.study[[j]][ apply(study.files[[j]][ names(treatments)[[j]] ] == as.data.frame(list) , 1, all) ]
                  groups[[j]] <- subgroups
                }
                                
              }else{
                n <- nrow(treatments[[j]])
                
                for(i in seq_len(n)){              
                  treatment <- data.frame(treatments[[j]][i,])
                  df <- data.frame(lapply(treatment, function(x) rep(x, each = length(samples.per.study[[j]]))))            
                  subgroups[[i]] = samples.per.study[[j]][ apply(study.files[[j]][ names(treatments[[j]])] == df, 1, all) ]
                }
                groups[[j]] <- subgroups
              }            
                
              }
            }else{
              message("Groups slot will be an empty list")
            }    
            
           .Object["groups"] <- groups
            return(.Object)
          }
)


setGeneric("setAssayTabs",function(.Object){standardGeneric("setAssayTabs")})
setMethod("setAssayTabs",
          signature(.Object = "ISATab"),
          function (.Object) 
          {
            
            atabs <- list()
            
            afiles <- .Object["assay.files"]
            assay.names <- lapply( afiles, function(i) i[ grep(isatab.syntax$assay.name, colnames(i) ) ]) 
            assay.filenames <- .Object["assay.filenames"]
            assay.tech.types <- .Object["assay.technology.types"]
                             
            assay.meas.types <- .Object["assay.measurement.types"]
            data.filenames <- .Object["data.filenames"]
            
            study.filenames <- .Object["study.filenames"]
            study.identifiers <- .Object["study.identifiers"]            
                       
            for(i in seq_len(length(assay.filenames))) {
              
              k <- getStudyFilenameIndex(.Object, assay.filenames[[i]])
              study.filename <- study.filenames[[k]]
              study.identifier <- study.identifiers[[k]]
                     
              if (class(data.filenames[[i]])!="data.frame"){
                data.filenames[[i]] <- as.data.frame(data.filenames[[i]])
                index <- grep("Data File", colnames(afiles[[i]]))
                colnames(data.filenames[[i]]) <- colnames(afiles[[i]])[[index]]
              }
              
              
              if (assay.tech.types[[i]]=="mass spectrometry"){
                
                atabs[[i]] <- new("MSAssayTab",
                                  path=.Object["path"],
                                  study.filename=study.filename,
                                  study.identifier=study.identifier,
                                  assay.filename=assay.filenames[[i]],
                                  assay.file=afiles[[i]],
                                  assay.technology.type=assay.tech.types[[i]],
                                  assay.measurement.type=assay.meas.types[[i]],
                                  assay.names=assay.names[[i]],
                                  data.filenames=data.filenames[[i]]
                )
                
                
              }else if (assay.tech.types[[i]]=="DNA microarray"){
                                
                
                atabs[[i]] <- new("MicroarrayAssayTab",
                                  path=.Object["path"],
                                  study.filename=study.filename,
                                  study.identifier=study.identifier,
                                  assay.filename=assay.filenames[[i]],
                                  assay.file=afiles[[i]],
                                  assay.technology.type=assay.tech.types[[i]],
                                  assay.measurement.type=assay.meas.types[[i]],
                                  assay.names=assay.names[[i]],
                                  data.filenames=data.filenames[[i]]
                )
                
                
              }else if (assay.tech.types[[i]]=="NMR spectroscopy"){
                
                
                atabs[[i]] <- new("NMRAssayTab",
                                  path=.Object["path"],
                                  study.filename=study.filename,
                                  study.identifier=study.identifier,
                                  assay.filename=assay.filenames[[i]],
                                  assay.file=afiles[[i]],
                                  assay.technology.type=assay.tech.types[[i]],
                                  assay.measurement.type=assay.meas.types[[i]],
                                  assay.names=assay.names[[i]],
                                  data.filenames=data.filenames[[i]]
                )
                
                
              }else{
                                
                
                atabs[[i]] <- new("AssayTab",
                                  path=.Object["path"],
                                  study.filename=study.filename,
                                  study.identifier=study.identifier,
                                  assay.filename=assay.filenames[[i]],
                                  assay.file=afiles[[i]],
                                  assay.technology.type=assay.tech.types[[i]],
                                  assay.measurement.type=assay.meas.types[[i]],
                                  assay.names=assay.names[[i]],
                                  data.filenames=data.filenames[[i]]                                                    
                )               
                
              }
              
            } #for
            
            .Object["assay.tabs"] <- atabs
            
            return(.Object)
          }
)



############################################################################################################
##### getRawDataFilenames
############################################################################################################
# generic method called 'getAssayRawDataFilenames' that
# dispatches on the type of object it's applied to
setGeneric(
  "getRawDataFilenames",
  function(.Object, full.path=TRUE) {
    standardGeneric("getRawDataFilenames")
  }
)

setMethod(
  "getRawDataFilenames",
  signature(.Object = "ISATab", full.path = "logical"),
  function(.Object, full.path=TRUE) {    
    assay.tabs <- .Object["assay.tabs"]    
    raw.data.filenames <- sapply(assay.tabs, function(x) getAssayRawDataFilenames(x, full.path))
    return(raw.data.filenames)
}
)

############################################################################################################
##### getAssayRawDataFilenames
############################################################################################################

# generic method called 'getAssayRawDataFilenames' that
# dispatches on the type of object it's applied to
setGeneric(
    "getAssayRawDataFilenames",
    function(.Object, full.path) {
      standardGeneric("getAssayRawDataFilenames")
    }
)
  
setMethod(
  "getAssayRawDataFilenames",
  signature(.Object = "MSAssayTab", full.path = "logical"),
  function(.Object, full.path=TRUE) {
    
    msfiles <- as.list(.Object["data.filenames"][isatab.syntax$raw.spectral.data.file])
    if (full.path)
      msfiles <- sapply(msfiles, function(x) sapply(x, function(y) paste(.Object["path"], y, sep=.Platform$file.sep)))  
    return(msfiles)
  }
)


setMethod(
  "getAssayRawDataFilenames",
  signature(.Object = "MicroarrayAssayTab", full.path ="logical"),
  function(.Object, full.path=TRUE) {
    
    #if (!isMicroarrayAssay(isa, assay.filename))
    #  stop("The ", assay.filename, " is not a microarray assay")
    
    microarray.files <- as.list(.Object["data.filenames"][isatab.syntax$array.data.file])
    if (full.path)
      microarray.files <- sapply(microarray.files, function(x) sapply(x, function(y) paste(.Object["path"], y, sep=.Platform$file.sep)))  
    return(microarray.files) 
  }
)

setMethod(
  "getAssayRawDataFilenames",
  signature(.Object = "SeqAssayTab", full.path = "logical"),
  function(.Object, full.path=TRUE) {
    
    msfiles <- as.list(.Object["data.filenames"][isatab.syntax$raw.data.file])
    if (full.path)
      msfiles <- sapply(msfiles, function(x) sapply(x, function(y) paste(.Object["path"], y, sep=.Platform$file.sep)))  
    return(msfiles)
  }
)

setMethod(
  "getAssayRawDataFilenames",
  signature(.Object = "AssayTab", full.path = "logical"),
  function(.Object, full.path=TRUE) {
    
    raw.files <- as.list(.Object["data.filenames"][isatab.syntax$raw.data.file])
    if (full.path)
      msfiles <- sapply(raw.files, function(x) sapply(x, function(y) paste(.Object["path"], y, sep=.Platform$file.sep)))  
    return(msfiles)
  }
)

setMethod(
  "getAssayRawDataFilenames",
  signature(.Object = "NMRAssayTab", full.path = "logical"),
  function(.Object, full.path=TRUE) {
    
    raw.files <- as.list(.Object["data.filenames"][isatab.syntax$free.induction.decay.data.file])
    if (full.path)
      msfiles <- sapply(raw.files, function(x) sapply(x, function(y) paste(.Object["path"], y, sep=.Platform$file.sep)))  
    return(msfiles)
  }
)

############################################################################################################
##### getAssayFilenames
############################################################################################################

setGeneric(
  "getAssayNames",
  function(.Object, full.path) {
    standardGeneric("getAssayNames")
  }
)

setMethod(
  "getAssayNames",
  signature(.Object = "MicroarrayAssayTab", full.path ="logical"),
  function(.Object, full.path=TRUE) {
    #  assay.filenames <- isa["assay.filenames"]
    #  assay.files <- isa["assay.files"]
    #  microarray.assay.filenames <- assay.filenames[ sapply(assay.files, function(x) isatab.syntax$hybridization.assay.name %in% names(x)) ]
    #  return(microarray.assay.filenames)
    
    assay.filename <- .Object["assay.filename"]
    assay.file <- .Object["assay.file"]
    assay.names <- assay.file[ isatab.syntax$hybridization.assay.name]
}
)



#getMicroarrayAssayFilenames <- function(isa){
#  assay.filenames <- isa["assay.filenames"]
#  assay.files <- isa["assay.files"]
#  microarray.assay.filenames <- assay.filenames[ sapply(assay.files, function(x) isatab.syntax$hybridization.assay.name %in% names(x)) ]
#  return(microarray.assay.filenames)
#}

############################################################################################################
##### getDerivedDataFilenames
############################################################################################################

setGeneric(
  "getDerivedDataFilenames",
  function(.Object, full.path) {
    standardGeneric("getDerivedDataFilenames")
  }
)

