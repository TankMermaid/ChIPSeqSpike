\docType{methods}
\name{plotCor}
\alias{plotCor}
\alias{plotCor,ChIPSeqSpikeDataset-method}
\alias{plotCor,ChIPSeqSpikeDatasetBoost-method}
\alias{plotCor,ChIPSeqSpikeDatasetList-method}
\alias{plotCor,ChIPSeqSpikeDatasetListBoost-method}

\title{Plot the correlation between ChIP-seq experiments}

\description{
    Plot the correlation between ChIP-seq experiments using heatscatter plot or,
 if heatscatterplot = FALSE, correlation tables.
}

\usage{

  plotCor(theObject, rawFile = FALSE, rpmFile = FALSE, bgsubFile = FALSE, 
revFile = FALSE, spiked = TRUE, main = "", add_contour = FALSE, method_cor = 
"spearman", nlevels = 10, color_contour = "black", show_cor = TRUE, 
allOnPanel = TRUE, method_scale = "none", method_corrplot = "circle", 
heatscatterplot = TRUE, type_corrplot = "upper", diag_corrplot = FALSE, 
separateWindows = FALSE, verbose = FALSE, ...)

\S4method{plotCor}{ChIPSeqSpikeDataset}(theObject, rawFile = FALSE, rpmFile = 
FALSE, bgsubFile = FALSE, revFile = FALSE, spiked = TRUE, main = "", 
add_contour = FALSE, method_cor = "spearman", nlevels = 10, color_contour = 
"black", show_cor = TRUE, allOnPanel = TRUE, method_scale = "none", 
method_corrplot = "circle", heatscatterplot = TRUE, type_corrplot = "upper", 
diag_corrplot = FALSE, separateWindows = FALSE, verbose = FALSE, ...)

\S4method{plotCor}{ChIPSeqSpikeDatasetBoost}(theObject, rawFile = FALSE, 
rpmFile = FALSE, bgsubFile = FALSE, revFile = FALSE, spiked = TRUE, main = "", 
add_contour = FALSE, method_cor = "spearman", nlevels = 10, 
color_contour = "black", show_cor = TRUE, allOnPanel = TRUE, method_scale = 
"none", method_corrplot = "circle", heatscatterplot = TRUE, type_corrplot = 
"upper", diag_corrplot = FALSE, separateWindows = FALSE, verbose = FALSE, ...)

\S4method{plotCor}{ChIPSeqSpikeDatasetList}(theObject, rawFile = FALSE, rpmFile
= FALSE, bgsubFile = FALSE, revFile = FALSE, spiked = TRUE, main = "", 
add_contour = FALSE, method_cor = "spearman", nlevels = 10, color_contour = 
"black", show_cor = TRUE, allOnPanel = TRUE, method_scale = "none",
 method_corrplot = "circle", heatscatterplot = TRUE, type_corrplot = "upper", 
 diag_corrplot = FALSE, separateWindows = FALSE, verbose = FALSE, ...)

\S4method{plotCor}{ChIPSeqSpikeDatasetListBoost}(theObject, rawFile = FALSE, 
rpmFile = FALSE, bgsubFile = FALSE, revFile = FALSE, spiked = TRUE, main = "", 
add_contour = FALSE, method_cor = "spearman", nlevels = 10, 
color_contour = "black", show_cor = TRUE, allOnPanel = TRUE, 
method_scale = "none", method_corrplot = "circle", heatscatterplot = TRUE, 
type_corrplot = "upper", diag_corrplot = FALSE, separateWindows = FALSE, 
verbose = FALSE, ...)
}

\arguments{
    \item{theObject}{\code{ChIPSeqSpike} dataset (see ?spikeDataset)}

    \item{rawFile}{If TRUE, use the untransformed data. Not available in boost 
                   mode. Default is FALSE. (see details)}

    \item{rpmFile}{If TRUE, use the RPM scaled data. Not available in boost 
                   mode. Default is FALSE. (see details)}

    \item{bgsubFile}{If TRUE, use the input subtracted data. Not available in 
                     boost mode. Default is FALSE. (see details)}

    \item{revFile}{If TRUE, use the RPM reverted data. Not available in boost 
                   mode. Default is FALSE. (see details)}

    \item{spiked}{If TRUE, use the spiked data. Default is TRUE. (see details)}

    \item{main}{Main title of the plot. No title is displayed by default.}

    \item{add_contour}{If TRUE, contours are added to the heatscatter. Default 
                       is FALSE.}

    \item{method_cor}{A character string indicating which correlation 
                     coefficient is to be computed.  One of 'pearson' (default)
                      , 'kendall' or 'spearman'.}

    \item{nlevels}{An integer giving the number of levels of the contour lines.
                   Not used if heatscatteplot is FALSE. Default is 10}

    \item{color_contour}{Character string defining the color of the contour 
                         line. Not used if heatscatteplot is FALSE. Default is 
                         'black'.}

    \item{show_cor}{Logical indicating if the correlation is added to the 
                    title. Not used if heatscatteplot is FALSE. Default is 
                    TRUE.}

    \item{allOnPanel}{Logical indicating if all correlations should be on the 
                     same panel. Not used if heatscatteplot is FALSE. Default 
                     is TRUE.}

    \item{method_scale}{Character string indicating the scaling to be applied 
                        to the data. Possible values are 'none', 'log', 
                        'asinh', 'cuberoot' or 'zscore'. Not used if 
                        heatscatteplot is FALSE.}

    \item{method_corrplot}{If heatscatterplot is FALSE, define the graphical 
                     representation used for the correlation table. Possible 
                     values are 'circle', 'square', 'ellipse', 'number', 'pie',
                     'shade' and 'color'. See ?corrplot::corrplot for more 
                     details. Default is 'circle'.}

    \item{heatscatterplot}{If TRUE, use a heatscatter representation instead of
                     correlation table. see ?LSD::heatscatter for more details.
                     Default is TRUE.}

    \item{type_corrplot}{If heatscatter is FALSE, define if the full ('full'), 
                         the lower triangular ('lower') or upper triangular 
                         matrix is displayed. Default is 'upper'.}

    \item{diag_corrplot}{If heatscatter is FALSE, logical indicating if the 
                   correlation coefficients are displayed on the principal 
                   diagonal. Default is FALSE.}

    \item{separateWindows}{If heatscatterplot is TRUE, Logical indicating if 
                          each plot is output to a separate window. Default is 
                          FALSE.} 

    \item{verbose}{Logical indicating if processing messages are displayed. 
                   Default is FALSE}

    \item{...}{Additional parameter to pass to the LSD::heatscatter or 
              corrplot::corrplot functions.}
}

\value{
    If heatscatterplot is FALSE, return the correlation matrix.
}

\details{
  The 'rawFile', 'rpmFile', 'bgsubFile', 'revFile', 'spiked' parameters 
indicate if the untransformed, RPM scaled, input DNA subtracted, RPM reversed 
or spiked data should be used plotted. This option is only available if not in 
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
## Sub-sample
data("result_extractBinding");

## Heatscatter of spiked data using Spearman correlation
plotCor(csds)

## Pearson correlation of log transformed raw data
plotCor(csds,  rawFile = TRUE, spiked = FALSE,  main = "heatscatter",  
method_cor = "pearson", method_scale = "log");
 
## Correlation table of all transformation steps with circle representation
plotCor(csds, rawFile = TRUE, rpmFile = TRUE, bgsubFile = TRUE, revFile = TRUE,
 spiked = TRUE, heatscatterplot = FALSE, verbose = TRUE);
 
## Correlation table of all transformation steps with number representation
plotCor(csds, rawFile = TRUE, rpmFile = TRUE, bgsubFile = TRUE, revFile = TRUE,
 spiked = TRUE, heatscatterplot = FALSE, verbose = TRUE, method_corrplot = 
 "number")
}

\seealso{
    \code{\link{spikeDataset}}
    \code{\link[LSD]{heatscatter}}
    \code{\link[corrplot]{corrplot}}
    \code{\link{spikePipe}}
    \code{\link{plotTransform}}
    \code{\link{boxplotSpike}}
    \code{\link{plotHeatmaps}}
    \code{\link{plotProfile}}
}

\author{Nicolas Descostes}
