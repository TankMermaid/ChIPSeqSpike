\docType{methods}
\name{datasetList}
\alias{datasetList}
\alias{datasetList<-}
\alias{datasetList<-,ChIPSeqSpikeDatasetList-method}
\alias{datasetList<-,ChIPSeqSpikeDatasetListBoost-method}


\title{Set the list of dataset representing experiments}

\description{
    Set the list of ChIPSeqSpikeDataset or ChIPSeqSpikeDatasetBoost objects.
}

\usage{
datasetList(theObject) <- value
            
            
\S4method{datasetList}{ChIPSeqSpikeDatasetList}(theObject) <- value

\S4method{datasetList}{ChIPSeqSpikeDatasetListBoost}(theObject) <- value
}

\arguments{
    \item{theObject}{A ChIPSeqSpikeDatasetList or ChIPSeqSpikeDatasetListBoost 
    object}
 
    \item{value}{A list of ChIPSeqSpikeDataset or ChIPSeqSpikeDatasetBoost 
    objects}
}

\value{
  The modified object is returned
}

\examples{
data(result_extractBinding)
datasetList(csds) <- getDatasetList(csds)
}

\seealso{
    \code{\link{getDatasetList}}
}

\author{Nicolas Descostes}
