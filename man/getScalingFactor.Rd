\docType{methods}
\name{getScalingFactor}
\alias{getScalingFactor}
\alias{getScalingFactor,ChIPSeqSpikeDataset-method}
\alias{getScalingFactor,ChIPSeqSpikeDatasetBoost-method}
\alias{getScalingFactor,Experiment-method}
\alias{getScalingFactor,ExperimentLoaded-method}
\alias{getScalingFactor,ChIPSeqSpikeCore-method}


\title{Get the endogenous scaling factor}

\description{
    Accesses and returns the experiment or input DNA experiment endogenous 
    scaling factor.
}

\usage{
getScalingFactor(theObject)
            
            
\S4method{getScalingFactor}{ChIPSeqSpikeDataset}(theObject)

\S4method{getScalingFactor}{ChIPSeqSpikeDatasetBoost}(theObject)

\S4method{getScalingFactor}{Experiment}(theObject)

\S4method{getScalingFactor}{ExperimentLoaded}(theObject)

\S4method{getScalingFactor}{ChIPSeqSpikeCore}(theObject)
}

\arguments{
    \item{theObject}{A \code{ChIPSeqSpikeCore}, \code{ChIPSeqSpikeDataset}, 
                     \code{ChIPSeqSpikeDatasetBoost}, \code{Experiment} or 
                     \code{ExperimentLoaded} object}
}

\value{
  A numeric of the endogenous scaling factor
}
\details{
  If the object is ChIPSeqSpikeCore, ChIPSeqSpikeDataset or 
  ChIPSeqSpikeDatasetBoost, getScalingFactor returns the input DNA experiment 
  endogenous scaling factor.
  
  If the object is Experiment or ExperimentLoaded, getScalingFactor returns 
  the experiment endogenous scaling factor.
}

\examples{
data(result_extractBinding)
getScalingFactor(csds[[1]])
}

\seealso{
    \code{\link{getExogenousScalingFactor}}
}

\author{Nicolas Descostes}
