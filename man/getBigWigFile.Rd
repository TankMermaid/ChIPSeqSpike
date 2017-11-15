\docType{methods}
\name{getBigWigFile}
\alias{getBigWigFile}
\alias{getBigWigFile,ChIPSeqSpikeDataset-method}
\alias{getBigWigFile,ChIPSeqSpikeDatasetList-method}
\alias{getBigWigFile,ChIPSeqSpikeDatasetBoost-method}
\alias{getBigWigFile,Experiment-method}
\alias{getBigWigFile,ExperimentLoaded-method}
\alias{getBigWigFile,ChIPSeqSpikeCore-method}


\title{Get the path to an endogenous experiment bigwig file}

\description{
    Access and returns the path to the endogenous bigwig file of an experiment
}

\usage{
getBigWigFile(theObject)
            
            
\S4method{getBigWigFile}{ChIPSeqSpikeDataset}(theObject)

\S4method{getBigWigFile}{ChIPSeqSpikeDatasetList}(theObject)

\S4method{getBigWigFile}{ChIPSeqSpikeDatasetBoost}(theObject)

\S4method{getBigWigFile}{Experiment}(theObject)

\S4method{getBigWigFile}{ExperimentLoaded}(theObject)

\S4method{getBigWigFile}{ChIPSeqSpikeCore}(theObject)
}

\arguments{
    \item{theObject}{A \code{ChIPSeqSpikeCore}, \code{ChIPSeqSpikeDataset}, 
                    \code{ChIPSeqSpikeDatasetList}, 
                    \code{ChIPSeqSpikeDatasetBoost}, \code{Experiment} or 
                    \code{ExperimentLoaded} object}
}

\value{
  A string of the path to the endogenous bigwig file
}
\details{
  If the object is ChIPSeqSpikeCore, ChIPSeqSpikeDataset or 
  ChIPSeqSpikeDatasetBoost, getBigWigFile returns the path to the endogenous 
  input DNA experiment bigwig file.
  
  If the object is Experiment or ExperimentLoaded, getBigWigFile returns the 
  path to the endogenous experiment bigwig file.
  
  If the object is ChIPSeqSpikeDatasetList, all bigwig files associated with 
  all dataset are returned at once.
}

\examples{
data(result_extractBinding);
getBigWigFile(csds);
}

\author{Nicolas Descostes}
