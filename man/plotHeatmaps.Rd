\docType{methods}
\name{plotHeatmaps}
\alias{plotHeatmaps}
\alias{plotHeatmaps,ChIPSeqSpikeDataset-method}
\alias{plotHeatmaps,ChIPSeqSpikeDatasetBoost-method}
\alias{plotHeatmaps,ChIPSeqSpikeDatasetList-method}
\alias{plotHeatmaps,ChIPSeqSpikeDatasetListBoost-method}

\title{plot heatmaps of ChIP-seq experiments}

\description{
    Clusters and output binding values signal in the form of heatmaps
}

\usage{

  plotHeatmaps(theObject, location = "start", 
                transformType = "spiked", legend = TRUE, plot_scale = "no",
                sort_rows = "decreasing", nb_of_groups = 1, 
                clustering_method = "none", include_exp_vec = NULL, 
                auto_scale = FALSE, raster_value = TRUE, col_value = "blue",
                ...)

\S4method{plotHeatmaps}{ChIPSeqSpikeDataset}(theObject, location = "start", 
transformType = "spiked", legend = TRUE, plot_scale = "no", sort_rows = 
"decreasing", nb_of_groups = 1, clustering_method = "none", include_exp_vec = 
NULL, auto_scale = FALSE, raster_value = TRUE, col_value = "blue", ...)

\S4method{plotHeatmaps}{ChIPSeqSpikeDatasetBoost}(theObject, location = 
"start", transformType = "spiked", legend = TRUE, plot_scale = "no", 
sort_rows = "decreasing", nb_of_groups = 1, clustering_method = "none", 
include_exp_vec = NULL, auto_scale = FALSE, raster_value = TRUE, 
col_value = "blue", ...)

\S4method{plotHeatmaps}{ChIPSeqSpikeDatasetList}(theObject, location = "start",
transformType = "spiked", legend = TRUE, plot_scale = "no", sort_rows = 
"decreasing", nb_of_groups = 1, clustering_method = "none", include_exp_vec = 
NULL, auto_scale = FALSE, raster_value = TRUE, col_value = "blue", ...)

\S4method{plotHeatmaps}{ChIPSeqSpikeDatasetListBoost}(theObject, location = 
"start", transformType = "spiked", legend = TRUE, plot_scale = "no", 
sort_rows = "decreasing", nb_of_groups = 1, clustering_method = "none", 
include_exp_vec = NULL, auto_scale = FALSE, raster_value = TRUE, 
col_value = "blue", ...)
}

\arguments{
    \item{theObject}{\code{ChIPSeqSpike} dataset (see ?spikeDataset)}

    \item{location}{Define the location at which heatmaps should be plotted. 
                    Possible values are 'start', 'midpoint', 'end' and 
                    'composite' (see details). Default is 'start'.}

    \item{transformType}{Indicate what step of chip-seq spike-in normalization 
                         is used. Possible values are 'spiked', 'reverse', 
                         'BGSub' 'RPM' and 'raw' (see details). Default is 
                         'spiked'. Not available in boost mode.}

    \item{legend}{Logical indicating if legend should be indicated. Default is 
                  FALSE.}

    \item{plot_scale}{Indicate the transformation that is applied to the data 
                     before plotting. Possible values are 'no' (default), 
                     'linear', 'log2' or 'zscore'. See ?seqplots::plotHeatmap 
                     for more details.}

    \item{sort_rows}{Indicate how rows of the heatmap should be sorted. 
                     Possible values are 'increasing', 'decreasing' or FALSE. 
                     See ?seqplots::plotHeatmap for more details. Default is 
                     'decreasing'.}

    \item{nb_of_groups}{If clustering_method is different than 'none', define 
                       the number of groups of the cluster. Default is 1.}

    \item{clustering_method}{Determine the heatmap clustering algorithm. 
                            Possible values are 'k-means', 'hclust', 'ssom', 
                            'bed_scores' and 'none'. See 
                            ?seqplots::plotHeatmap for more details. Default is
                             'none'.}

    \item{include_exp_vec}{Logical vector indicating on which experiments 
                         clustering is performed. NULL meaning all experiments.
                         See ?seqplots::plotHeatmap for more details. Default 
                         is NULL.}

    \item{auto_scale}{Logical indicating if color scaled should be specific to 
                      each experiment. Default is FALSE.}

    \item{raster_value}{Logical indicating if the bitmap raster is used. See 
                       ?seqplots::plotHeatmap for more details. Default is 
                       FALSE.}

    \item{col_value}{The vector or list of colour values used generate 
                     sub-heatmaps colorspaces. See ?seqplots::plotHeatmap for 
                     more details. Default is 'blue'.}

    \item{...}{Additional parameter to pass to seqplots::plotHeatmap.}
}

\value{
    See ?seqplots::plotHeatmap for details.
}

\details{
  This method relies on the 'plotHeatmap' method of the 'seqplots' package and 
enables to generate different heatmaps on the chIPSeqSpike objects (see 
?spikeDataset). 

  The 'transformType' indicates if the untransformed (raw), RPM scaled (RPM), 
input subtracted (BGSub), RPM reversed (reverse) or spiked (spiked) data should
 be used to generate the heatmaps. This option is only available if not in 
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
data("result_extractBinding");

## Spiked data in decreasing order
plotHeatmaps(csds);

## Raw data in decreasing order
plotHeatmaps(csds, transformType = "raw");

## K-means clustering with 3 groups
plotHeatmaps(csds, nb_of_groups = 3, clustering_method = "kmeans");
}

\seealso{
    \code{\link{spikeDataset}}
    \code{\link{spikePipe}}
    \code{\link[seqplots]{plotHeatmap}}
    \code{\link{plotProfile}}
    \code{\link{plotTransform}}
    \code{\link{boxplotSpike}}
    \code{\link{plotCor}}
}

\author{Nicolas Descostes}
