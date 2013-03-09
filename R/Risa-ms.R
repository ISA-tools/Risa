### Methods to deal with assays whose technology type is mass spectrometry 

### specific function to deal with assays whose technology type is mass spectrometry using the xcms package
### it returns an xcmsSet
processAssayXcmsSet.1factor = function(isa, assay.filename, ...){
  
  i <- which(isa["assay.filenames"]==assay.filename)
  
  #if 'Raw Spectral Data File' is one of the columns in the assay file = it is a mass spectrometry assay
  if (isatab.syntax$raw.spectral.data.file %in% colnames(isa["data.filenames"][[i]]))
  {
    #mass spectrometry files
    msfiles = isa["data.filenames"][[i]][[ isatab.syntax$raw.spectral.data.file ]]
    
    #the assay file as an AnnotatedDataFrame
    pd = try(read.AnnotatedDataFrame(file.path(isa["path"], isa["assay.filenames"][i]),
                                     row.names = NULL, blank.lines.skip = TRUE, fill = TRUE,
                                     varMetadata.char = "$", quote="\""))
    
    #Adding the raw spectral data files as the row names
    sampleNames(pd) = pd$Raw.Spectral.Data.File
    
    if (length(grep(isatab.syntax$factor.value, colnames(isa["assay.files"][[i]]))) != 0) {
      ## If there are explicit factors, use the first one of them
      sclass = isa["assay.files"][[i]][ which(isa["assay.files"][[i]][[isatab.syntax$sample.name]] %in% pd$Sample.Name), grep(isatab.syntax$factor.value, colnames(isa["assay.files"][[i]]))[1]]
      
      wd <- getwd()
      setwd(isa["path"])
      xset = xcmsSet(files=msfiles, sclass=sclass, ...)
      setwd(wd)
      
    } else {
      wd <- getwd()
      setwd(isa["path"])
      ## Otherwise just use what was there
      xset = try(xcmsSet(msfiles, phenoData=pData(pd), ...))
      setwd(wd)
    }
    return(xset)
  }#if
   
}#processAssayXcmsSet



### specific function to deal with assays whose technology type is mass spectrometry using the xcms package
### it returns an xcmsSet, 
processAssayXcmsSet = function(isa, assay.filename, ...){
      
  i <- which(isa["assay.filenames"]==assay.filename)
  
  if (isatab.syntax$raw.spectral.data.file %in% colnames(isa["data.filenames"][[i]]))
  {    
              
    #mass spectrometry files
    msfiles = isa["data.filenames"][[i]][[ isatab.syntax$raw.spectral.data.file ]]        
    
    pd = try(read.AnnotatedDataFrame(file.path(isa["path"], isa["assay.filenames"][i]),
                                     row.names = NULL, blank.lines.skip = TRUE, fill = TRUE,
                                     varMetadata.char = "$", quote="\""))
    
    sampleNames(pd) = pd$Raw.Spectral.Data.File
    
    if (length(grep(isatab.syntax$factor.value, colnames(isa["assay.files"][[i]]))) != 0) {
      ## If there are explicit factors, use them
      sclass <- as.data.frame(isa["factors"])
      
      wd <- getwd()
      setwd(isa["path"])
      xset = xcmsSet(files=msfiles, sclass=sclass, ...)
          
      setwd(wd)
      
    } else {
      wd <- getwd()
      setwd(isa["path"])
      ## Otherwise just use what was there
      xset = try(xcmsSet(msfiles, phenoData=pData(pd), ...))
      setwd(wd)
    }
    return(xset)
  }#if
}#processAssayXcmsSet