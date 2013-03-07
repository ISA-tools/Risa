### Methods to deal with assays whose technology type is mass spectrometry 

### specific function to deal with assays whose technology type is mass spectrometry using the xcms package
### it returns an xcmsSet
processAssayXcmsSet.1factor = function(isa, assay.filename, ...){
  for(i in seq_len(length(isa["assay.filenames"]))){
    
    if (isa["assay.filenames"][[i]]==assay.filename){
      
      if (isatab.syntax$raw.spectral.data.file %in% colnames(isa["data.filenames"][[i]]))
      {
        #mass spectrometry files
        msfiles = isa["data.filenames"][[i]][[ isatab.syntax$raw.spectral.data.file ]]
        
        pd = try(read.AnnotatedDataFrame(file.path(isa["path"], isa["assay.filenames"][i]),
                                         row.names = NULL, blank.lines.skip = TRUE, fill = TRUE,
                                         varMetadata.char = "$", quote="\""))
        
        sampleNames(pd) = pd$Raw.Spectral.Data.File
        
        if (length(grep("Factor.Value", colnames(isa["assay.files"][[i]]))) != 0) {
          ## If there are explicit factors, use them
          sclass = isa["assay.files"][[i]][ which(isa["assay.files"][[i]][[isatab.syntax$sample.name]] %in% pd$Sample.Name), grep("Factor.Value", colnames(isa["assay.files"][[i]]))[1]]
          
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
      
    }#if 
  }#for
}#processAssayXcmsSet



### specific function to deal with assays whose technology type is mass spectrometry using the xcms package
### it returns an xcmsSet
### TODO - change implementation to include all factors
processAssayXcmsSet = function(isa, assay.filename, ...){
  
  
  phenodata.data.frame <- as.data.frame(isa["factors"])
  
  
  
  for(i in seq_len(length(isa["assay.filenames"]))){
    
    if (isa["assay.filenames"][[i]]==assay.filename){
      
      if (isatab.syntax$raw.spectral.data.file %in% colnames(isa["data.filenames"][[i]]))
      {
        #mass spectrometry files
        msfiles = isa["data.filenames"][[i]][[ isatab.syntax$raw.spectral.data.file ]]
                
        pd = try(read.AnnotatedDataFrame(file.path(isa["path"], isa["assay.filenames"][i]),
                                         row.names = NULL, blank.lines.skip = TRUE, fill = TRUE,
                                         varMetadata.char = "$", quote="\""))
        
        sampleNames(pd) = pd$Raw.Spectral.Data.File
        
        if (length(grep("Factor.Value", colnames(isa["assay.files"][[i]]))) != 0) {
          ## If there are explicit factors, use them
          sclass = isa["assay.files"][[i]][ which(isa["assay.files"][[i]][[isatab.syntax$sample.name]] %in% pd$Sample.Name), grep("Factor.Value", colnames(isa["assay.files"][[i]]))[1]]
          
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
      
    }#if 
  }#for
}#processAssayXcmsSet