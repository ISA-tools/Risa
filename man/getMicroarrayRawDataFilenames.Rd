\name{getMicroarrayRawDataFilenames}
\alias{getMicroarrayRawDataFilenames}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Retrieves the raw data filenames of microarray assays.
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
getMicroarrayRawDataFilenames(isa, full.path = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{isa}{
%%     ~~Describe \code{isa} here~~
}
  \item{full.path}{
%%     ~~Describe \code{full.path} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (isa, full.path = TRUE) 
{
    microarray.assay.filenames <- getMicroarrayAssayFilenames(isa)
    microarray.files <- lapply(isa["data.filenames"][microarray.assay.filenames], 
        function(x) x[isatab.syntax$array.data.file])
    if (full.path) 
        microarray.files <- sapply(microarray.files, function(x) sapply(x, 
            function(y) paste(isa["path"], y, sep = .Platform$file.sep)))
    return(microarray.files)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
