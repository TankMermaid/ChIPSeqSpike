\docType{methods}
\name{getCount}
\alias{getCount}
\alias{getCount,ChIPSeqSpikeDataset-method}
\alias{getCount,ChIPSeqSpikeDatasetBoost-method}
\alias{getCount,Experiment-method}
\alias{getCount,ExperimentLoaded-method}
\alias{getCount,ChIPSeqSpikeCore-method}


\title{Get the number of reads aligned to the endogenous reference genome}

\description{
    Access and returns the number of reads of an experiment or input DNA 
    experiment that were aligned to the endogenous reference genome.
}

\usage{
getCount(theObject)
            
            
\S4method{getCount}{ChIPSeqSpikeDataset}(theObject)

\S4method{getCount}{ChIPSeqSpikeDatasetBoost}(theObject)

\S4method{getCount}{Experiment}(theObject)

\S4method{getCount}{ExperimentLoaded}(theObject)

\S4method{getCount}{ChIPSeqSpikeCore}(theObject)
}

\arguments{
    \item{theObject}{A \code{ChIPSeqSpikeCore}, \code{ChIPSeqSpikeDataset}, 
                     \code{ChIPSeqSpikeDatasetBoost}, \code{Experiment} or 
                     \code{ExperimentLoaded} object}
}

\value{
  A numeric of the number of reads aligned to the endogenous reference genome
}
\details{
  If the object is ChIPSeqSpikeCore, ChIPSeqSpikeDataset or 
  ChIPSeqSpikeDatasetBoost, getCount returns the number of reads of the input 
  DNA experiment that were aligned to the endogenous reference genome.
  
  If the object is Experiment or ExperimentLoaded, getCount returns the number 
  of reads of the experiment that were aligned to the endogenous reference 
  genome.
}

\examples{
data(result_extractBinding);
getCount(csds[[1]]);
}

\seealso{
    \code{\link{getExoCount}}
}

\author{Nicolas Descostes}
