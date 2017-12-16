\docType{methods}
\name{getExperimentList}
\alias{getExperimentList}
\alias{getExperimentList,ChIPSeqSpikeDataset-method}
\alias{getExperimentList,ChIPSeqSpikeDatasetBoost-method}


\title{Get all Experiment or ExperimentLoaded objects associated with a 
dataset}

\description{
    Accesses and returns all Experiment or ExperimentLoaded objects associated 
    with a dataset.
}

\usage{
getExperimentList(theObject)
            
            
\S4method{getExperimentList}{ChIPSeqSpikeDataset}(theObject)

\S4method{getExperimentList}{ChIPSeqSpikeDatasetBoost}(theObject)
}

\arguments{
    \item{theObject}{A \code{ChIPSeqSpikeDataset} or 
                     \code{ChIPSeqSpikeDatasetBoost} object}
}

\value{
  A list of Experiment or ExperimentLoaded objects
}

\examples{
data(result_extractBinding)
getExperimentList(csds[[1]])
}

\seealso{
    \code{\link{experimentList}}
}

\author{Nicolas Descostes}
