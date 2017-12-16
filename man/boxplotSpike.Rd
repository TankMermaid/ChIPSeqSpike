\docType{methods}
\name{boxplotSpike}
\alias{boxplotSpike}
\alias{boxplotSpike,ChIPSeqSpikeDataset-method}
\alias{boxplotSpike,ChIPSeqSpikeDatasetBoost-method}
\alias{boxplotSpike,ChIPSeqSpikeDatasetList-method}
\alias{boxplotSpike,ChIPSeqSpikeDatasetListBoost-method}

\title{plot boxplots of ChIP-seq experiments}

\description{
    plot boxplots of the mean values of ChIP-seq experiments on the annotations
    given to the extractBinding method
}

\usage{
  boxplotSpike(theObject, col = NULL, rawFile = FALSE, 
                rpmFile = FALSE, bgsubFile = FALSE, revFile = FALSE, 
                spiked = TRUE, ylab = NULL, outline = TRUE, violinPlot = FALSE,
                notch = TRUE, mean_with_sd = FALSE, mean = FALSE, 
                median = FALSE, boxplot = FALSE, jitter = FALSE, plot = TRUE,
                verbose = FALSE)

\S4method{boxplotSpike}{ChIPSeqSpikeDataset}(theObject, col = NULL, rawFile = 
FALSE, rpmFile = FALSE, bgsubFile = FALSE, revFile = FALSE, spiked = TRUE, ylab
= NULL, outline = TRUE, violinPlot = FALSE, notch = TRUE, mean_with_sd = FALSE,
 mean = FALSE, median = FALSE, boxplot = FALSE, jitter = FALSE, plot = TRUE,
verbose = FALSE)

\S4method{boxplotSpike}{ChIPSeqSpikeDatasetBoost}(theObject, col = NULL, 
rawFile = FALSE, rpmFile = FALSE, bgsubFile = FALSE, revFile = FALSE, 
spiked = TRUE, ylab = NULL, outline = TRUE, violinPlot = FALSE, notch = TRUE, 
mean_with_sd = FALSE, mean = FALSE, median = FALSE, boxplot = FALSE, 
jitter = FALSE, plot = TRUE, verbose = FALSE)

\S4method{boxplotSpike}{ChIPSeqSpikeDatasetList}(theObject, col = NULL, rawFile
= FALSE, rpmFile = FALSE, bgsubFile = FALSE, revFile = FALSE, spiked = TRUE, 
ylab = NULL, outline = TRUE, violinPlot = FALSE, notch = TRUE, mean_with_sd = 
FALSE, mean = FALSE, median = FALSE, boxplot = FALSE, jitter = FALSE, 
plot = TRUE, verbose = FALSE)

\S4method{boxplotSpike}{ChIPSeqSpikeDatasetListBoost}(theObject, col = NULL, 
rawFile = FALSE, rpmFile = FALSE, bgsubFile = FALSE, revFile = FALSE, 
spiked = TRUE, ylab = NULL, outline = TRUE, violinPlot = FALSE, notch = TRUE, 
mean_with_sd = FALSE, mean = FALSE, median = FALSE, boxplot = FALSE, 
jitter = FALSE, plot = TRUE, verbose = FALSE)
}

\arguments{
    \item{theObject}{\code{ChIPSeqSpike} dataset (see ?spikeDataset)}

    \item{col}{Vector of colors for each experiment. Default is NULL.}

    \item{rawFile}{If TRUE, use the untransformed data. Not available in boost 
                  mode. Default is FALSE. (see details)}

    \item{rpmFile}{If TRUE, use the RPM scaled data. Not available in boost 
                   mode. Default is FALSE. (see details)}

    \item{bgsubFile}{If TRUE, use the input subtracted data. Not available in 
                     boost mode. Default is FALSE. (see details)}

    \item{revFile}{If TRUE, use the RPM reverted data. Not available in boost 
                   mode. Default is FALSE. (see details)}

    \item{spiked}{If TRUE, use the spiked data. Default is TRUE. (see details)}

    \item{ylab}{Character string of the name of the y-axis. Default is NULL.}

    \item{outline}{Logical indicating if outliers are shown. Default is TRUE.}

    \item{violinPlot}{Logical indicating if a violin plot representation is 
                      used. Default is FALSE.}

    \item{notch}{Logical indicating if confidence intervals are shown. Default 
                 is TRUE.}

    \item{mean_with_sd}{Logical indicating if the mean and standard deviation 
                        are shown on the violin plot. Default is FALSE.}

    \item{mean}{Logical indicating if the mean is shown on the violin plot. 
                Default is FALSE.}

    \item{median}{Logical indicating if the median is shown on the violin plot.
                  Default is FALSE.}

    \item{boxplot}{Logical indicating if boxplot is shown on the violin plot. 
                   Default is FALSE}

    \item{jitter}{Logical indicating if each mean values is represented as a 
                  point on the violin plot. Default is FALSE.}

    \item{plot}{Logical indicating if the boxplot should be plotted. Default is
               TRUE.}

    \item{verbose}{Logical indicating if processing messages are shown. Default
                   is FALSE.}
}

\value{
    Return a list with the components 'stats', 'n', 'conf', 'out', 'group' and 
'names'. See ?boxplot for details.
}

\details{
  The 'rawFile', 'rpmFile', 'bgsubFile', 'revFile', 'spiked' parameters 
indicate if the untransformed, RPM scaled, input DNA subtracted, RPM reversed 
or spiked data should be plotted. This option is only available if not in 
boost mode (see ?spikePipe for details on what these steps are).
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
                                                 'ChIPSeqSpikeDatasetListBoost'
}}

\examples{
data("result_extractBinding")
boxplotSpike(csds)
boxplotSpike(csds, outline = FALSE, violinPlot = TRUE)
boxplotSpike(csds, outline = FALSE, violinPlot = TRUE, mean = TRUE, jitter = 
TRUE)
boxplotSpike(csds, rawFile = TRUE, rpmFile = TRUE, bgsubFile = TRUE, revFile = 
TRUE, spiked = TRUE, outline = FALSE, violinPlot = TRUE)
}

\seealso{
    \code{\link{spikeDataset}}
    \code{\link{spikePipe}}
    \code{\link{plotProfile}}
    \code{\link{plotTransform}}
    \code{\link{plotHeatmaps}}
    \code{\link{plotCor}}
}

\author{Nicolas Descostes}
