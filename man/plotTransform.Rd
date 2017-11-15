\docType{methods}
\name{plotTransform}
\alias{plotTransform}
\alias{plotTransform,ChIPSeqSpikeDataset-method}
\alias{plotTransform,ChIPSeqSpikeDatasetList-method}

\title{Plots average profiles of steps of the spike-in normalization}

\description{
    For each experiment, plots average profiles of RPM scaled, input 
subtracted, RPM reversed and spiked data at different annotation locations (not
 available in boost mode)
}

\usage{

  plotTransform(theObject, legends = FALSE, colVec = NULL, 
                separateWindows = FALSE)

\S4method{plotTransform}{ChIPSeqSpikeDataset}(theObject, legends = FALSE, 
colVec = NULL, separateWindows = FALSE)

\S4method{plotTransform}{ChIPSeqSpikeDatasetList}(theObject, legends = FALSE, 
colVec = NULL, separateWindows = FALSE)
}

\arguments{
    \item{theObject}{A \code{ChIPSeqSpikeDatasetList} or 
                                         \code{ChIPSeqSpikeDataset} object}

    \item{legends}{Logical indicating if legend should be indicated. Default is
                   FALSE.}

    \item{colVec}{Character vector indicating the colors to use. Default is 
                  NULL}

    \item{separateWindows}{Plot each experiment in a separate window. Default 
                           is FALSE.}
}

\value{
    Nothing
}

\details{
    plotTransform plots average signal over annotations that were given to the 
'extractBinding' method. plotTransform is using the 'plotAverage' method of 
the 'seqplots' package.

   The signal is plotted at four different annotation locations: 'start', 
   'end', 'midpoint' and 'composite' (pf, ef, mf and af options of 
   'plotAverage' method).

  As objects created in boost mode only hold the binding values in GRanges 
objects, the previously applied transformations are not kept in memory. 
Therefore, this method does not work with ChIPSeqSpikeDatasetBoost and 
ChIPSeqSpikeDatasetListBoost objects.
}

\section{Methods (by class)}{
\itemize{

\item \code{ChIPSeqSpikeDataset}: Method for signature theObject=
                                                          'ChIPSeqSpikeDataset'

\item \code{ChIPSeqSpikeDatasetList}: Method for signature theObject=
                                                      'ChIPSeqSpikeDatasetList'
}}

\examples{
data("result_extractBinding");
plotTransform(csds, TRUE);
plotTransform(csds, TRUE, separateWindows=TRUE);
}

\seealso{
    \code{\link{spikeDataset}}
    \code{\link{extractBinding}}
    \code{\link[seqplots]{plotAverage}}
    \code{\link{plotProfile}}
    \code{\link{boxplotSpike}}
    \code{\link{plotHeatmaps}}
    \code{\link{plotCor}}
}

\author{Nicolas Descostes}
