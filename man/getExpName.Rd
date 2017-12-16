\docType{methods}
\name{getExpName}
\alias{getExpName}
\alias{getExpName,ChIPSeqSpikeDataset-method}
\alias{getExpName,ChIPSeqSpikeDatasetBoost-method}
\alias{getExpName,Experiment-method}
\alias{getExpName,ExperimentLoaded-method}


\title{Get the experiment name}

\description{
    Accesses and returns the experiment names associated with a dataset.
}

\usage{
getExpName(theObject)
            
            
\S4method{getExpName}{ChIPSeqSpikeDataset}(theObject)

\S4method{getExpName}{ChIPSeqSpikeDatasetBoost}(theObject)

\S4method{getExpName}{Experiment}(theObject)

\S4method{getExpName}{ExperimentLoaded}(theObject)
}

\arguments{
    \item{theObject}{A \code{ChIPSeqSpikeDataset}, 
                     \code{ChIPSeqSpikeDatasetBoost}, \code{Experiment} or 
                     \code{ExperimentLoaded} object}
}

\value{
  A string or character vector of the names of the experiments 
}

\examples{
data(result_extractBinding)
getExpName(csds[[1]])
}

\author{Nicolas Descostes}
