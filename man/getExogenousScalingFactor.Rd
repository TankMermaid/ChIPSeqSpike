\docType{methods}
\name{getExogenousScalingFactor}
\alias{getExogenousScalingFactor}
\alias{getExogenousScalingFactor,Experiment-method}
\alias{getExogenousScalingFactor,ExperimentLoaded-method}

\title{Get the exogenous scaling factor}

\description{
    Accesses and returns the experiment exogenous scaling factor.
}

\usage{
getExogenousScalingFactor(theObject)
            
            
\S4method{getExogenousScalingFactor}{Experiment}(theObject)

\S4method{getExogenousScalingFactor}{ExperimentLoaded}(theObject)
}

\arguments{
    \item{theObject}{An \code{Experiment} or \code{ExperimentLoaded} object}
}

\value{
  A numeric of the exogenous scaling factor
}

\examples{
data(result_extractBinding);
getExogenousScalingFactor(csds[[1]][[1]]);
}

\seealso{
    \code{\link{getScalingFactor}}
}

\author{Nicolas Descostes}
