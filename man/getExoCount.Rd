\docType{methods}
\name{getExoCount}
\alias{getExoCount}
\alias{getExoCount,Experiment-method}
\alias{getExoCount,ExperimentLoaded-method}


\title{Get the number of reads aligned to the exogenous reference genome}

\description{
    Access and returns the number of reads of an experiment that were aligned 
    to the exogenous reference genome.
}

\usage{
getExoCount(theObject)
            
            
\S4method{getExoCount}{Experiment}(theObject)

\S4method{getExoCount}{ExperimentLoaded}(theObject)
}

\arguments{
    \item{theObject}{An \code{Experiment} or \code{ExperimentLoaded} object}
}

\value{
  A numeric of the number of reads aligned to the exogenous reference genome
}

\examples{
data(result_extractBinding);
getExoCount(csds[[1]][[1]]);
}

\seealso{
    \code{\link{getCount}}
}

\author{Nicolas Descostes}
