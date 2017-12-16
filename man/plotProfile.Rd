\docType{methods}
\name{plotProfile}
\alias{plotProfile}
\alias{plotProfile,ChIPSeqSpikeDataset-method}
\alias{plotProfile,ChIPSeqSpikeDatasetBoost-method}
\alias{plotProfile,ChIPSeqSpikeDatasetList-method}
\alias{plotProfile,ChIPSeqSpikeDatasetListBoost-method}

\title{Plots average profiles of ChIP-seq experiments}

\description{
    Plots average profiles of all experiments contained in a 'spikeDataset' at 
different locations
}

\usage{

    plotProfile(theObject, legends = FALSE, colVec = NULL, notScaled = FALSE)

\S4method{plotProfile}{ChIPSeqSpikeDataset}(theObject, legends = FALSE, 
colVec = NULL, notScaled = FALSE)

\S4method{plotProfile}{ChIPSeqSpikeDatasetBoost}(theObject, legends = FALSE, 
colVec = NULL)

\S4method{plotProfile}{ChIPSeqSpikeDatasetList}(theObject, legends = FALSE, 
colVec = NULL, notScaled = FALSE)

\S4method{plotProfile}{ChIPSeqSpikeDatasetListBoost}(theObject, legends = 
 FALSE, colVec = NULL)
}

\arguments{
    \item{theObject}{\code{ChIPSeqSpike} dataset (see ?spikeDataset)}

    \item{legends}{Logical indicating if legend should be indicated. Default is
                  FALSE.}

    \item{colVec}{Character vector indicating the colors to use. Default is 
                  NULL}

    \item{notScaled}{If not in boost mode, logical indicating if untransformed 
                    data should be plotted. Default is FALSE.}
}

\value{
    Nothing
}

\details{
  plotProfile plots average signal over annotations that were given to the 
'extractBinding' method. plotProfile is using the 'plotAverage' method of 
the 'seqplots' package.

  By default, the spiked signal is plotted. If the 'notScaled' parameter is
set to TRUE and the object is of type ChIPSeqSpikeDatasetList or 
ChIPSeqSpikeDataset (no boost mode), the RPM scaled and input subtracted values
are also plotted. This option enables to visualize the effect of spike-in 
scaling.

  The signal is plotted at four different annotation locations: 'start', 'end',
'midpoint' and 'composite' (pf, ef, mf and af options of 'plotAverage' method).
}

\section{Methods (by class)}{
\itemize{

\item \code{ChIPSeqSpikeDataset}: Method for signature theObject=
                                                          'ChIPSeqSpikeDataset'

\item \code{ChIPSeqSpikeDatasetBoost}: Method for signature theObject=
                                                     'ChIPSeqSpikeDatasetBoost'

\item \code{ChIPSeqSpikeDatasetList}: Method for signature theObject=
                                                      'ChIPSeqSpikeDatasetList'

\item \code{ChIPSeqSpikeDatasetListBoost}: Method for signature theObject=
                                                  ChIPSeqSpikeDatasetListBoost'
}}

\examples{
data("result_extractBinding")
plotProfile(csds)
plotProfile(csds, TRUE)
plotProfile(csds, TRUE, notScaled = TRUE)
}

\seealso{
    \code{\link{spikeDataset}}
    \code{\link{extractBinding}}
    \code{\link[seqplots]{plotAverage}}
    \code{\link{plotTransform}}
    \code{\link{boxplotSpike}}
    \code{\link{plotHeatmaps}}
    \code{\link{plotCor}}
}

\author{Nicolas Descostes}
