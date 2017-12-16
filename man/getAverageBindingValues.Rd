\docType{methods}
\name{getAverageBindingValues}
\alias{getAverageBindingValues}
\alias{getAverageBindingValues,ChIPSeqSpikeDataset-method}
\alias{getAverageBindingValues,ChIPSeqSpikeDatasetBoost-method}
\alias{getAverageBindingValues,ChIPSeqSpikeCore-method}


\title{Get the average binding values associated to a dataset}

\description{
    Accessor returning the average binding values associated to a dataset.
}

\usage{
getAverageBindingValues(theObject)
            
            
\S4method{getAverageBindingValues}{ChIPSeqSpikeDataset}(theObject)

\S4method{getAverageBindingValues}{ChIPSeqSpikeDatasetBoost}(theObject)

\S4method{getAverageBindingValues}{ChIPSeqSpikeCore}(theObject)
}

\arguments{
    \item{theObject}{A \code{ChIPSeqSpikeCore}, \code{ChIPSeqSpikeDataset} or 
                     \code{ChIPSeqSpikeDatasetBoost} object.}
 }

\value{
  A list of plotSetArray objects.
}

\details{
  Each element of the list contains the binding values to perform
meta-profile (see ?plotProfile); transformation profiles if not in boost mode 
(see ?plotTransform) and heatmaps (see ?plotHeatmaps). These values are stored 
in a plotSetArray object. This object is created by the method getPlotSetArray 
of the 'seqplots' package.
}

\examples{
data(result_extractBinding)
getAverageBindingValues(csds[[1]])
}

\seealso{
    \code{\link{plotProfile}}
    \code{\link{plotTransform}}
    \code{\link{plotHeatmaps}}
    \code{\link{getMatBindingValues}}
    \code{\link[seqplots]{getPlotSetArray}}
}

\author{Nicolas Descostes}
