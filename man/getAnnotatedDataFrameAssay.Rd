\name{getAnnotatedDataFrameAssay}
\alias{getAnnotatedDataFrameAssay}
\title{
Retrieves an AnnotatedDataFrameAssay.
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
getAnnotatedDataFrameAssay(isa, assay.filename)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{isa}{
%%     ~~Describe \code{isa} here~~
}
  \item{assay.filename}{
%%     ~~Describe \code{assay.filename} here~~
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
function (isa, assay.filename) 
{
    i <- which(names(isa["assay.files"]) == assay.filename)
    dataf <- as.data.frame(isa@factors[[i]])
    colnames(dataf) <- names(isa@factors[[i]])
    row.names(dataf) <- isa@samples
    dataf.desc <- as.data.frame(names(isa@factors[[i]]))
    pdata <- new("AnnotatedDataFrame", data = dataf, varMetadata = dataf.desc)
    return(pdata)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
