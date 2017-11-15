\docType{methods}
\name{spikeSummary}
\alias{spikeSummary}
\alias{spikeSummary,Experiment-method}
\alias{spikeSummary,ExperimentLoaded-method}
\alias{spikeSummary,ChIPSeqSpikeDataset-method}
\alias{spikeSummary,ChIPSeqSpikeDatasetBoost-method}
\alias{spikeSummary,ChIPSeqSpikeDatasetList-method}
\alias{spikeSummary,ChIPSeqSpikeDatasetListBoost-method}

\title{Output dataset summary information}

\description{
    Output a table giving endogenous/exogenous scaling scores and read counts 
for each experiment contained in a dataset 
}

\usage{

    spikeSummary(theObject)

\S4method{spikeSummary}{Experiment}(theObject)
\S4method{spikeSummary}{ExperimentLoaded}(theObject)
\S4method{spikeSummary}{ChIPSeqSpikeDataset}(theObject)
\S4method{spikeSummary}{ChIPSeqSpikeDatasetList}(theObject)
\S4method{spikeSummary}{ChIPSeqSpikeDatasetBoost}(theObject)
\S4method{spikeSummary}{ChIPSeqSpikeDatasetListBoost}(theObject)
}

\arguments{
    \item{theObject}{\code{ChIPSeqSpike} dataset (see ?spikeDataset)}
}

\value{
    a matrix of information about each experiment
}

\details{
  The rows represent an experiment (with corresponding input below) and the 
columns represent the endogenous scaling factor (endoScalFact), the exogenous 
scaling factor (exoScalFact), the number of reads aligned to the reference 
genome (endoCount) and the number of reads aligned to the exogenous genome 
(exoCount).
}

\section{Methods (by class)}{
\itemize{
\item \code{Experiment}: Method for signature theObject = 'Experiment'

\item \code{ExperimentLoaded}: Method for signature theObject = 
                                                             'ExperimentLoaded'

\item \code{ChIPSeqSpikeDataset}: Method for signature theObject=
                                                          'ChIPSeqSpikeDataset'

\item \code{ChIPSeqSpikeDatasetBoost}: Method for signature theObject=
                                                     'ChIPSeqSpikeDatasetBoost'

\item \code{ChIPSeqSpikeDatasetList}: Method for signature theObject=
                                                      'ChIPSeqSpikeDatasetList'

\item \code{ChIPSeqSpikeDatasetListBoost}: Method for signature theObject=
                                                 'ChIPSeqSpikeDatasetListBoost'
}}

\examples{
data("result_estimateScalingFactors");
spikeSummary(csds);
}

\seealso{
    \code{\link{getRatio}}
}

\author{Nicolas Descostes}
