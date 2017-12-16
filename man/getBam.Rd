\docType{methods}
\name{getBam}
\alias{getBam}
\alias{getBam,ChIPSeqSpikeDataset-method}
\alias{getBam,ChIPSeqSpikeDatasetBoost-method}
\alias{getBam,Experiment-method}
\alias{getBam,ExperimentLoaded-method}
\alias{getBam,ChIPSeqSpikeCore-method}


\title{Get the path to an endogenous experiment bam file}

\description{
    Access and returns the path to the bam file containing the reads of an 
    experiment aligned to the endogenous reference genome.
}

\usage{
getBam(theObject)
            
            
\S4method{getBam}{ChIPSeqSpikeDataset}(theObject)

\S4method{getBam}{ChIPSeqSpikeDatasetBoost}(theObject)

\S4method{getBam}{Experiment}(theObject)

\S4method{getBam}{ExperimentLoaded}(theObject)

\S4method{getBam}{ChIPSeqSpikeCore}(theObject)
}

\arguments{
    \item{theObject}{A \code{ChIPSeqSpikeCore}, \code{ChIPSeqSpikeDataset}, 
                    \code{ChIPSeqSpikeDatasetBoost}, \code{Experiment} or 
                    \code{ExperimentLoaded} object}
}

\value{
  A string of the path to the endogenous bam file
}
\details{
  If the object is ChIPSeqSpikeCore, ChIPSeqSpikeDataset or 
  ChIPSeqSpikeDatasetBoost, getBam returns the path to the input DNA experiment
  bam file containing reads aligned to the endogenous reference genome.
  
  If the object is Experiment or ExperimentLoaded, getBam returns the path to 
  the experiment bam file containing reads aligned to the endogenous reference 
  genome.
}

\examples{
data(result_extractBinding)
getBam(csds[[1]])
}

\seealso{
    \code{\link{getExogenousBam}}
}

\author{Nicolas Descostes}
