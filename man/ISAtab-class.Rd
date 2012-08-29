\name{ISAtab-class}
\docType{class}
\alias{ISAtab-class}
\alias{[,ISAtab-method]}
\title{Class ISAtab, a class for maintaing ISAtab information}
\description{
  This class maintains informatino from an ISAtab dataset, including an investigation file, one or more study files, and one or more assay files for each study file. 
}
\section{Objects from the Class}{
  Objects can be created with the \code{ISAtab} constructor. Objects can also be
  created by calls of the form \code{new("ISAtab", ...)}.
}
\section{Slots}{
  \describe{
    \item{\code{path}:}{
      the path of the ISA-Tab dataset
    }
    \item{\code{investigation.filename}:}{
      the name of the Investigation file
    }
    \item{\code{investigation.file}:}{
       a data frame with the contents of the Investigation file
    }
    \item{\code{study.identifiers}:}{
       the list of study identifiers
    }
    \item{\code{study.filenames}:}{
       the names of the study files
    }
    \item{\code{study.files}:}{
       a list of data frames wiht the contents of the study files
    }
    \item{\code{assay.filenames}:}{
      the names of the assay files
    }
    \item{\code{assay.filenames.per.study}:}{
      the names of the assay files according to the study they belong to
    }
  \item{\code{assay.files}}{
	    a list of data frames with the contents of the assay files
	}
    \item{\code{assay.files.per.study}:}{
     a list of data frames with the contents of the assay files divided per study they belong to
    }
    \item{\code{assay.technology.types}:}{   
     a list with the technology types corresponding to each assay
    }
    \item{\code{assay.measurement.types}:}{
       a list with the names of the data files
    }
    \item{\code{data.filenames}:}{
     a list with the names of the samples
    }
     \item{\code{samples}:}{
     the samples classified according to the assay filename they belong to
    }
    \item{\code{samples.per.study}:}{
     the names of the assay files classified per sample name
    }
    \item{\code{samples.per.assay.filename}:}{
    
    }
    \item{\code{assay.filenames.per.sample}:}{
    
    }
     \item{\code{sample.to.rawdatafile}:}{
       a list of data frames with samples associated with raw data files
    }
    \item{\code{sample.to.assayname}:}{
       maintains the association between samples and assay names
    }
      \item{\code{rawdatafile.to.sample}:}{
       mantains the association raw data file and samples
    }
    \item{\code{assayname.to.sample}:}{
       mantains the associations between assay names and samples
    }
  }
}


    
    
    
    
  
  
    

\section{Methods}{
  \describe{
    \item{\link[xcms:c.xcmsSet]{c}}{
      \code{signature("xcmsSet")}: combine objects together
    }
    \item{filepaths<-}{
      \code{signature(object = "xcmsSet")}: set \code{filepaths} slot
    }
    \item{filepaths}{
      \code{signature(object = "xcmsSet")}: get \code{filepaths} slot
    }
    \item{\link{diffreport}}{
      \code{signature(object = "xcmsSet")}: create report of
      differentially regulated ions including EICs
    }
    \item{\link{fillPeaks}}{
      \code{signature(object = "xcmsSet")}: fill in peak data for
      groups with missing peaks
    }
    \item{\link{getEIC}}{
      \code{signature(object = "xcmsSet")}: get list of EICs for
      each sample in the set
    }
    \item{groupidx<-}{
      \code{signature(object = "xcmsSet")}: set \code{groupidx} slot
    }
    \item{groupidx}{
      \code{signature(object = "xcmsSet")}: get \code{groupidx} slot
    }
    \item{\link{groupnames}}{
      \code{signature(object = "xcmsSet")}: get textual names for
      peak groups
    }
    \item{groups<-}{
      \code{signature(object = "xcmsSet")}: set \code{groups} slot
    }
    \item{groups}{
      \code{signature(object = "xcmsSet")}: get \code{groups} slot
    }
    \item{\link{groupval}}{
      \code{signature(object = "xcmsSet")}: get matrix of values
      from peak data with a row for each peak group
    }
    \item{\link{group}}{
      \code{signature(object = "xcmsSet")}: find groups of peaks
      across samples that share similar m/z and retention times
    }
    \item{peaks<-}{
      \code{signature(object = "xcmsSet")}: set \code{peaks} slot
    }
    \item{peaks}{
      \code{signature(object = "xcmsSet")}: get \code{peaks} slot
    }
    \item{\link{plotrt}}{
      \code{signature(object = "xcmsSet")}: plot retention time
      deviation profiles
    }
    \item{profinfo<-}{
      \code{signature(object = "xcmsSet")}: set \code{profinfo} slot
    }
    \item{profinfo}{
      \code{signature(object = "xcmsSet")}: get \code{profinfo} slot
    }
    \item{\link{retcor}}{
      \code{signature(object = "xcmsSet")}: use initial grouping
      of peaks to do nonlinear loess retention time correction
    }
    \item{sampclass<-}{
      \code{signature(object = "xcmsSet")}: \strong{DEPRECATED}. If used, the
      experimental design will be replaced with a data frame with a
      single column matching the supplied factor.
    }
    \item{sampclass}{
      \code{signature(object = "xcmsSet")}: get the interaction of the
      experimental design factors
    }
    \item{phenoData<-}{
      \code{signature(object = "xcmsSet")}: set the \code{phenoData} slot
    }
    \item{phenoData}{
      \code{signature(object = "xcmsSet")}: get the \code{phenoData} slot
    }
    \item{progressCallback<-}{
      \code{signature(object = "xcmsSet")}: set the \code{progressCallback} slot
    }
    \item{progressCallback}{
      \code{signature(object = "xcmsSet")}: get the \code{progressCallback} slot
    }
    \item{sampnames<-}{
      \code{signature(object = "xcmsSet")}: set rownames in the
      \code{phenoData} slot
    }
    \item{\link{sampnames}}{
      \code{signature(object = "xcmsSet")}: get rownames in the
      \code{phenoData} slot
    }
    \item{\link[xcms:split.xcmsSet]{split}}{
      \code{signature("xcmsSet")}: divide into a list of objects
    }
  }
}
\references{
  A parallel effort in metabolite profiling data sharing:
  \url{http://metlin.scripps.edu/}
}
\author{Colin A. Smith, \email{csmith@scripps.edu}}
\note{
  No notes yet.
}
\seealso{
  \code{\link{xcmsSet}}
}
\keyword{classes}
