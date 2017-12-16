\docType{methods}
\name{getExogenousBam}
\alias{getExogenousBam}
\alias{getExogenousBam,Experiment-method}
\alias{getExogenousBam,ExperimentLoaded-method}


\title{Get the path to an exoogenous experiment bam file}

\description{
    Accesses and returns the path to the bam file containing the reads of an 
    experiment aligned to the exogenous reference genome.
}

\usage{
getExogenousBam(theObject)
            
            

\S4method{getExogenousBam}{Experiment}(theObject)

\S4method{getExogenousBam}{ExperimentLoaded}(theObject)
}

\arguments{
    \item{theObject}{An \code{Experiment} or \code{ExperimentLoaded} object}
}

\value{
  A string of the path to the exogenous bam file
}

\examples{
data(result_extractBinding)
getExogenousBam(csds[[1]][[1]])
}

\seealso{
    \code{\link{getBam}}
}

\author{Nicolas Descostes}
