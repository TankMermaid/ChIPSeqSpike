\docType{methods}
\name{getExperimentListBigWigs}
\alias{getExperimentListBigWigs}
\alias{getExperimentListBigWigs,ChIPSeqSpikeDataset-method}
\alias{getExperimentListBigWigs,ChIPSeqSpikeDatasetBoost-method}


\title{Get all paths to the bigwig files associated with a dataset}

\description{
    Accesses and returns all paths to the bigwig files associated with a 
    dataset.
}

\usage{
getExperimentListBigWigs(theObject)
            
            
\S4method{getExperimentListBigWigs}{ChIPSeqSpikeDataset}(theObject)

\S4method{getExperimentListBigWigs}{ChIPSeqSpikeDatasetBoost}(theObject)
}

\arguments{
    \item{theObject}{A \code{ChIPSeqSpikeDataset} or 
                     \code{ChIPSeqSpikeDatasetBoost} object}
}

\value{
  A character vector of all paths to the bigwig files associated with a dataset
}

\examples{
data(result_extractBinding)
getExperimentListBigWigs(csds[[1]])
}

\author{Nicolas Descostes}
