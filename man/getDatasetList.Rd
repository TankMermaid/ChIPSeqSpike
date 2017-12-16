\docType{methods}
\name{getDatasetList}
\alias{getDatasetList}
\alias{getDatasetList,ChIPSeqSpikeDatasetList-method}
\alias{getDatasetList,ChIPSeqSpikeDatasetListBoost-method}


\title{Get the list of ChIPSeqSpike objects}

\description{
    Access and returns the list of ChIPSeqSpikeDataset or 
    ChIPSeqSpikeDatasetBoost objects contained in the structure
}

\usage{
getDatasetList(theObject)
            
            
\S4method{getDatasetList}{ChIPSeqSpikeDatasetList}(theObject)

\S4method{getDatasetList}{ChIPSeqSpikeDatasetListBoost}(theObject)
}

\arguments{
    \item{theObject}{A \code{ChIPSeqSpikeDatasetList} or 
    \code{ChIPSeqSpikeDataseListBoost} object}
}

\value{
  A list of \code{ChIPSeqSpikeDataset} or \code{ChIPSeqSpikeDatasetBoost} 
  objects
}

\examples{
data(result_extractBinding)
result <- getDatasetList(csds)
is(result)
is(result[[1]])
}

\seealso{
    \code{\link{datasetList}}
}

\author{Nicolas Descostes}
