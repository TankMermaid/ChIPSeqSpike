\docType{methods}
\name{getExperimentList}
\alias{getExperimentList}
\alias{getExperimentList,ChIPSeqSpikeDataset-method}
\alias{getExperimentList,ChIPSeqSpikeDatasetBoost-method}


\title{Get all paths to the bigwig files associated with a dataset}

\description{
    Accesses and returns all paths to the bigwig files associated with a 
    dataset.
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
  A character vector of all paths to the bigwig files associated with a dataset
}

\examples{
data(result_extractBinding);
getExperimentList(csds[[1]]);
}

\author{Nicolas Descostes}
