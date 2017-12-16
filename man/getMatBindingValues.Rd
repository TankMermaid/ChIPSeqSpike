\docType{methods}
\name{getMatBindingValues}
\alias{getMatBindingValues}
\alias{getMatBindingValues,ChIPSeqSpikeDataset-method}
\alias{getMatBindingValues,ChIPSeqSpikeDatasetBoost-method}
\alias{getMatBindingValues,ChIPSeqSpikeCore-method}

\title{Get the list of matrices of binding scores}

\description{
    Accesses and returns the list of binding scores matrices of all experiments
    associated with a dataset. These matrices are used to plot boxplots and 
    correlation plots.
}

\usage{
getMatBindingValues(theObject)
            
            
\S4method{getMatBindingValues}{ChIPSeqSpikeDataset}(theObject)

\S4method{getMatBindingValues}{ChIPSeqSpikeDatasetBoost}(theObject)

\S4method{getMatBindingValues}{ChIPSeqSpikeCore}(theObject)
}

\arguments{
    \item{theObject}{A \code{ChIPSeqSpikeCore}, \code{ChIPSeqSpikeDataset} or 
                     \code{ChIPSeqSpikeDatasetBoost} object.}
}

\value{
  A list of matrices containing binding scores
}

\examples{
data(result_extractBinding)
getMatBindingValues(csds[[1]])
}

\seealso{
    \code{\link{getAverageBindingValues}}
    \code{\link{boxplotSpike}}
    \code{\link{plotCor}}
}

\author{Nicolas Descostes}
