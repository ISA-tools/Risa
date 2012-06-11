isatab2bioczip = function(zip, path = getwd())
{
  # TODO the paths only work if the zip file does not contain a directory
  d = unzip(zipfile = zip, exdir = extract <- path)
  isaobj = isatab2bioc(path)
  return(isaobj)
}##end function isatab2bioczip


isatab2bioc = function(path = getwd())
{
  d = dir(path)
  
  ## Investigation filename
  ifilename = grep("i_", d, value=TRUE)
  ## Reading in investigation file into a data frame
  ifile = read.table(file.path(path, ifilename), sep="\t", fill=TRUE)

  ## Assay filenames
  afilenames = unlist(sapply(ifile[grep("Study Assay File Name", ifile[,1], useBytes=TRUE),], function(i) grep("a_", i, value=TRUE, useBytes=TRUE)))

  ## Reading in assay files into data frames
  afiles = lapply(afilenames, function(i) read.table(file.path(path, i), sep="\t", header=TRUE, stringsAsFactors=FALSE))

  ## Assay technology types
  types = ifile[grep("Study Assay Technology Type$", ifile[,1], useBytes=TRUE),]    
  types = na.omit(types[types != ""])
  types = types[-1]

  ## Data filenames
  dfiles = lapply(afiles, function(i) i[,grep("Data.File", colnames(i))])

  ## Study filenames (one or more)
  sfilenames = unlist(sapply(ifile[grep("Study File Name", ifile[,1], useBytes=TRUE),], function(i) grep("s_", i, value=TRUE, useBytes=TRUE)))
  
  ## Reading study files
  sfiles = lapply(sfilenames, function(i) read.table(file.path(path, i), sep="\t", header=TRUE, stringsAsFactors=FALSE))
  
  ## Identifying what sample is studied in which assay
  assays = lapply(seq_len(length(afiles)), function(i) sfile$Sample.Name %in% afiles[[i]]$Sample.Name)
  assays = do.call(cbind, assays)
  for(i in seq_len(ncol(assays)))
    {
      assays[assays[,i]==TRUE,i] = paste("isa",i, sep="")
      assays[assays[,i]==FALSE,i] = ""
    }      

  ## Adding the study file content to the isa object
  metadata = cbind(sfile, assays)
	
  isa = list()

  for(i in seq_len(length(dfiles)))
    {
      #############################################################################
      if(types[i] == "DNA microarray")
        {
          ## Raw and processed data filenames
          rawfilenames = if("Array.Data.File" %in% colnames(dfiles[[i]])) dfiles[[i]][,"Array.Data.File"] else NULL
          procfilenames = if("Derived.Array.Data.File" %in% colnames(dfiles[[i]])) dfiles[[i]][,"Derived.Array.Data.File"] else NULL
			   
          ## URL for ADF (Array Design Format) file
          urladf = paste("http://www.ebi.ac.uk/microarray-as/ae/files/", unique(afiles[[i]][,"Array.Design.REF"]), "/", unique(afiles[[i]][,"Array.Design.REF"]), ".adf.txt", sep="")
          adffilename = file.path(path,unique(afiles[[i]][,"Array.Design.REF"]))
          adf_download = download.file(urladf, adffilename, mode="wb")

          ## List containing rawfiles, sdrf, idf, adf & directory containing the files
          ## as required by {ArrayExpress} magetab2bioc function
          files = list(path = path,
            rawfiles = rawfilenames,
            procfile = procfilenames,
            sdrf = afilenames[[i]],
            idf = ifilename,
            adf = basename(adffilename))
			   
          if (is.null(dim(dfiles[[i]])[2]))
            ## No processed files
            isa[[i]] = try(magetab2bioc(files)) 
          else {
              raw = try(magetab2bioc(files))
              ## TODO more testing
              ## The following issues an R CMD check warning,
              ## no visible binding for global variable ‘procol’
              ## likely to be a true issue:
              cn = getcolproc(files)
              procol = cn[1]
              proc = try(procset(files, procol = procol))
              isa[[i]] = list(raw=raw, processed=proc)}
        }## end microarray
      #############################################################################
      
      #############################################################################
      if(types[i] == "flow cytometry")
        {
          pd = try(read.AnnotatedDataFrame(file.path(path, afilenames[i]),row.names = NULL, blank.lines.skip = TRUE, fill = TRUE, varMetadata.char = "$", quote="\""))
          sampleNames(pd) = pd$Raw.Data.File
                                        #if(all(dfiles[[i]] %in% dir(path)) && all(dfiles[[i]] %in% sampleNames(pd)))
          isa[[i]] = try(suppressWarnings(read.flowSet(dfiles[[i]], phenoData=pd, path=path)))
                                        #if(!all(dfiles[[i]] %in% dir(path)))
                                        #isa[[i]] = try(suppressWarnings(read.flowSet(dfiles[[i]][dfiles[[i]] %in% dir(path)], phenoData=pd, path=path)))				
                                        #if(!all(dfiles[[i]] %in% sampleNames(pd)))
                                        #isa[[i]] = try(suppressWarnings(read.flowSet(dfiles[[i]][dfiles[[i]] %in% sampleNames(pd)], phenoData=pd, path=path)))				
        }## end flow cytometry
      #############################################################################
      
      #############################################################################
      if(types[i] == "mass spectrometry")
        {
          if("Raw.Spectral.Data.File" %in% colnames(dfiles[[i]]))
            {
              msfiles = dfiles[[i]]$Raw.Spectral.Data.File
              pd = try(read.AnnotatedDataFrame(file.path(path, afilenames[i]),
                row.names = NULL, blank.lines.skip = TRUE, fill = TRUE,
                varMetadata.char = "$", quote="\""))
              sampleNames(pd) = pd$Raw.Spectral.Data.File

              if(length(grep("Factor.Value", colnames(metadata))) != 0) {
                ## If there are explicit factors, use them
                sclass=metadata[which(metadata$Sample.Name %in% pd$Sample.Name),grep("Factor.Value", colnames(metadata))[1]]
                  isa[[i]] = xcmsSet(files=msfiles, sclass=sclass)
                } else {
                  ## Otherwise just use what was there
                  isa[[i]] = try(xcmsSet(msfiles, phenoData=pData(pd)))
                }

            }			
        }## end mass spectrometry
      #############################################################################
      
    }## end for on dfiles

		
  names(isa) = do.call(paste, list("isa", seq_len(length(isa)), sep=""))
  isaobj = list(metadata, isa)

  names(isaobj) = c("metadata","data")
  return(isaobj)

}##end function isatab2bioc
