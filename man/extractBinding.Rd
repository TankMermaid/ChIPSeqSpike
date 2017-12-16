\docType{methods}
\name{extractBinding}
\alias{extractBinding}
\alias{extractBinding,ChIPSeqSpikeDataset-method}
\alias{extractBinding,ChIPSeqSpikeDatasetBoost-method}
\alias{extractBinding,ChIPSeqSpikeDatasetList-method}
\alias{extractBinding,ChIPSeqSpikeDatasetListBoost-method}

\title{Extract binding values for graphical representations}

\description{
    Extracts and formats binding scores for each experiment into structures 
    adapted to performing different graphical representations.
}

\usage{
extractBinding(theObject, gff_vec, genome, binsize = 50, before = 2000, 
               after=2000, mean_or_median = "mean", interpolation_number = 100,
               interpolation_average = 10000, ignore_strand = FALSE, 
               verbose = FALSE)
            
\S4method{extractBinding}{ChIPSeqSpikeDataset}(theObject, gff_vec, genome, 
                binsize = 50, before = 2000, after=2000, 
                mean_or_median = "mean", interpolation_number = 100, 
                interpolation_average = 10000, ignore_strand = FALSE, 
                verbose = FALSE)

\S4method{extractBinding}{ChIPSeqSpikeDatasetBoost}(theObject, gff_vec, genome,
                binsize = 50, before = 2000, after=2000, 
                mean_or_median = "mean", interpolation_number = 100, 
                interpolation_average = 10000, ignore_strand = FALSE, 
                verbose = FALSE)

\S4method{extractBinding}{ChIPSeqSpikeDatasetList}(theObject, gff_vec, genome, 
                binsize = 50, before = 2000, after=2000, 
                mean_or_median = "mean", interpolation_number = 100, 
                interpolation_average = 10000, ignore_strand = FALSE, 
                verbose = FALSE)

\S4method{extractBinding}{ChIPSeqSpikeDatasetListBoost}(theObject, gff_vec, 
                genome, binsize = 50, before = 2000, after=2000, 
                mean_or_median = "mean", interpolation_number = 100, 
                interpolation_average = 10000, ignore_strand = FALSE, 
                verbose = FALSE)
}

\arguments{
    \item{theObject}{\code{ChIPSeqSpike} dataset (see ?spikeDataset)}

    \item{gff_vec}{File in GFF format containing annotations used to plot 
                   information}

    \item{genome}{The UCSC code of reference genome, e.g. 'hg19' for Homo 
                  sapiens (see details)}

    \item{binsize}{Binning size used to create bigwig files. Default is 50.}

    \item{before}{Length in bp of the interval upstream annotation. Default is 
                  2000.}

    \item{after}{Length in bp of the interval downstream annotation. Default is
                 2000.}

    \item{mean_or_median}{For average profiles, should the 'mean' or 'median' 
                          values be used. Default is 'mean'.}

    \item{interpolation_number}{Number of interpolated points to create 
                                matrices
                               (see details). Default is 100.}

    \item{interpolation_average}{Number of interpolated points of profiles and 
                           heatmaps (see details). Default is 10000.}

    \item{ignore_strand}{If TRUE, the directionality is ignored, that is all 
                         features strands, regardless of annotation in GFF 
                         file, are treated as undetermined ("*"). Default is 
                         FALSE.}

    \item{verbose}{If TRUE, output processing messages. Default is FALSE.}
}

\value{
    Returns the same object with binding values in the form of plotSetArray and
    matrices (see details).
}

\details{
    This method should be called before performing any graphical analysis. It 
  updates two slots of theObject:
  \itemize{
             \item \code{SetArrayList}: Contains the binding values to perform
meta-profile (see ?plotProfile); transformation profiles if not in boost mode 
(see ?plotTransform) and heatmaps (see ?plotHeatmaps). These values are stored 
in a plotSetArray object. This object is created by the method getPlotSetArray 
of the 'seqplots' package.

             \item \code{matBindingList}: Contains list of matrices for each 
experiment. Each row correspond to an annotation given by gff_vec and the 
number of columns is defined by the interpolation_number parameter. These 
matrices are used to perform boxplots (see ?boxplotSpike) and correlation plots
 (see ?plotCor).
  }


    The SetArrayList contains values for 4 kind of representations (profiles 
and heatmaps): Representation at the start of the annotation (-before/
+after parameters); at the midpoint of the annotation; at the end of the 
annotation (-before/+after) or at the entire annotation (-before/+after).
    For representations using the entire annotations and upstream (before)/
downstream intervals, the number of points used for the within annotation 
interpolation is defined by the interpolation_average parameter.

    For details on installing reference genomes, see details of the function 
    'getPlotSetArray' of the 'seqplots' package.
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

## Mock example on a restricted number of reads
info_file_csv <- system.file("extdata/info.csv", package="ChIPSeqSpike")
bam_path <- system.file(c("extdata/bam_files"), package="ChIPSeqSpike")
bigwig_path <- system.file(c("extdata/bigwig_files"), package="ChIPSeqSpike")
gff_vec <- system.file("extdata/test_coord.gff", package="ChIPSeqSpike")
genome_name <- "hg19"

csds <- spikeDataset(infoFile = info_file_csv, bamPath = bam_path, 
                     bigWigPath = bigwig_path)

## Copying test files to the current folder
originalBW_vec <- as.character(getBigWigFile(csds))
dir.create("./test_chipseqspike")
result <- file.copy(originalBW_vec, "test_chipseqspike")

csds <- estimateScalingFactors(csds)

## Apply RPM scaling
csds <- scaling(csds, outputFolder = "test_chipseqspike")

## Perform input subtraction
csds <- inputSubtraction(csds)

## Reverse RPM scaling after input subtraction
csds <- scaling(csds, reverse = TRUE)

## Apply exogenous scaling factors
csds <- scaling(csds, type = "exo")

## Extract binding values
csds <- extractBinding(csds, gff_vec, genome_name)

## Delete all files generated in this example
unlink("test_chipseqspike/", recursive = TRUE)
}

\seealso{
    \code{\link{spikeDataset}}
    \code{\link{plotProfile}}
    \code{\link{plotTransform}}
    \code{\link{plotHeatmaps}}
    \code{\link{boxplotSpike}}
    \code{\link{plotCor}}
    \code{\link[seqplots]{getPlotSetArray}}
}

\author{Nicolas Descostes}
