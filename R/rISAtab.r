risazip = function(zip, path = getwd())
{
  d = unzip(zipfile = zip, exdir = extract <- path)
  isaobj = risatab(path)
  return(isaobj)
}##end function risazip


risatab = function(path = getwd())
{
  d = dir(path)

  ifi = grep("i_", d, value=TRUE)
  ifile = read.table(file.path(path, ifi), sep="\t", fill=TRUE)

  af = unlist(sapply(ifile[grep("Study Assay File Name", ifile[,1]),], function(i) grep("a_", i, value=TRUE)))

  afiles = lapply(af, function(i) read.table(file.path(path, i), sep="\t", header=TRUE, stringsAsFactors=FALSE))

  types = ifile[grep("Study Assay Technology Type$", ifile[,1]),]      
  types = na.omit(types[types != ""])
  types = types[-1]

  dfiles = lapply(afiles, function(i) i[,grep("Data.File", colnames(i))])

  ## Reading study file
  sfile = read.table(file.path(path, grep("^s_", d, value=TRUE)), sep="\t", header=TRUE)

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
      if(types[i] == "DNA microarray")
        {
          rawfiles = if("Array.Data.File" %in% colnames(dfiles[[i]])) dfiles[[i]][,"Array.Data.File"] else NULL
          procfile = if("Derived.Array.Data.File" %in% colnames(dfiles[[i]])) dfiles[[i]][,"Derived.Array.Data.File"] else NULL
			   
          urladf = paste("http://www.ebi.ac.uk/microarray-as/ae/files/", unique(afiles[[i]][,"Array.Design.REF"]), "/", unique(afiles[[i]][,"Array.Design.REF"]), ".adf.txt", sep="")
          adffile = file.path(path,unique(afiles[[i]][,"Array.Design.REF"]))
          adf = download.file(urladf, adffile, mode="wb")


          files = list(path = path,
            rawfiles = rawfiles,
            procfile = procfile,
            sdrf = af[[i]],
            idf = ifi,
            adf = basename(adffile))
			   
          if(is.null(dim(dfiles[[i]])[2]))
            isa[[i]] = try(magetab2bioc(files)) else {
              raw = try(magetab2bioc(files))
              proc = try(procset(files, procol = procol))
              isa[[i]] = list(raw=raw, processed=proc)}
        }## end microarray

      if(types[i] == "flow cytometry")
        {
          pd = try(read.AnnotatedDataFrame(file.path(path, af[i]),row.names = NULL, blank.lines.skip = TRUE, fill = TRUE, varMetadata.char = "$", quote="\""))
          sampleNames(pd) = pd$Raw.Data.File
                                        #if(all(dfiles[[i]] %in% dir(path)) && all(dfiles[[i]] %in% sampleNames(pd)))
          isa[[i]] = try(suppressWarnings(read.flowSet(dfiles[[i]], phenoData=pd, path=path)))
                                        #if(!all(dfiles[[i]] %in% dir(path)))
                                        #isa[[i]] = try(suppressWarnings(read.flowSet(dfiles[[i]][dfiles[[i]] %in% dir(path)], phenoData=pd, path=path)))				
                                        #if(!all(dfiles[[i]] %in% sampleNames(pd)))
                                        #isa[[i]] = try(suppressWarnings(read.flowSet(dfiles[[i]][dfiles[[i]] %in% sampleNames(pd)], phenoData=pd, path=path)))				
        }## end flow cytometry

      if(types[i] == "mass spectrometry")
        {
          if("Raw.Spectral.Data.File" %in% colnames(dfiles[[i]]))
            {
              msfiles = dfiles[[i]]$Raw.Spectral.Data.File
              pd = try(read.AnnotatedDataFrame(file.path(path, af[i]),row.names = NULL, blank.lines.skip = TRUE, fill = TRUE, varMetadata.char = "$", quote="\""))
              sampleNames(pd) = pd$Raw.Spectral.Data.File
              if(length(grep("Factor.Value", colnames(metadata))) != 0)
                {
                  isa[[i]] = try(xcmsSet(file.path(path,msfiles), phenoData=pData(pd), sclass= metadata[which(metadata$Sample.Name %in% pd$Sample.Name),grep("Factor.Value", colnames(metadata))[1]]))
                } else isa[[i]] = try(xcmsSet(file.path(path,msfiles), phenoData=pData(pd)))

            }			
        }## end mass spectrometry

    }## end for on dfiles

		
  names(isa) = do.call(paste, list("isa", seq_len(length(isa)), sep=""))
  isaobj = list(metadata, isa)

  names(isaobj) = c("metadata","data")
  return(isaobj)

}##end function risatab
