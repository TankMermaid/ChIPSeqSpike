\name{spikePipe}
\alias{spikePipe}

\title{ChIP-seq spike-in normalization wrapper function}

\description{
  This function performs all steps of spike-in normalization: Dataset creation,
  RPM scaling, input DNA subtraction, RPM scaling reversal, exogenous DNA 
  scaling (spike) and binding values extraction.
}

\usage{
    spikePipe(infoFile, bamPath, bigWigPath, anno, genome_version, 
              paired = FALSE, binsize = 50, profile_length_before = 2000, 
              profile_length_after= 2000, mean_or_median = "mean", 
              interpolation_number = 100, interpolation_average = 10000,
              ignore_strand = FALSE, verbose = FALSE, boost = FALSE, 
              outputFolder = NULL)
}

\arguments{
    \item{infoFile}{csv or tab separated txt file containing information about 
                    files (see details)}

    \item{bamPath}{Path to the folder containing bam files}

    \item{bigWigPath}{Path to the folder containing bigwig files}

    \item{anno}{File in GFF format containing annotations used to plot 
                information}

    \item{genome_version}{The UCSC code of reference genome, e.g. 'hg19' for 
                          Homo sapiens (see details)}

    \item{paired}{Indicate if sequences are single- or paired-ended. Default is
                 FALSE}

    \item{binsize}{Binning size used to create bigwig files. Default is 50.}

    \item{profile_length_before}{Length in bp of the interval upstream 
                                 annotation (see details). Default is 2000.}

    \item{profile_length_after}{Length in bp of the interval downstream 
                                annotation (see details). Default is 2000.}

    \item{mean_or_median}{For average profiles, should the 'mean' or 'median' 
                          values be used. Default is 'mean'.}

    \item{interpolation_number}{Number of interpolated points to create 
                               matrices (see details). Default is 100.}

    \item{interpolation_average}{Number of interpolated points of profiles and 
                           heatmaps (see details). Default is 10000.}

    \item{ignore_strand}{If TRUE the directionality is ignored, that is all 
                         features strands, regardless of annotation in GFF 
                         file, are treated as undetermined ("*"). default is 
                         FALSE.}

    \item{verbose}{If TRUE, output processing messages. Default is FALSE.}

    \item{boost}{If TRUE, the object created enables to perform the analysis in
                 boost mode (see details). Default is FALSE}

    \item{outputFolder}{Define the folder where scaled bigwig are output. 
                        Default is NULL (see details).}
}


\details{
    'infoFile' should be a csv or a tab separated txt file. The column names 
should be: expName, endogenousBam, exogenousBam, inputBam, bigWigEndogenous and
bigWigInput. These columns indicate the experiment names; the bam file names of
 data aligned to the reference genome; the bam file names of data aligned to 
the exogenous genome; the input DNA bam file names corresponding to each 
experiment; the bigwig file names of data aligned to the reference genome and 
the bigwig file names of input DNA experiments.

    If 'infoFile' contains only one input file (specified for each experiment),
 a ChIPSeqSpikeDataset (or ChIPSeqSpikeDatasetBoost) object is created. if 
'infoFile' contains different input DNA files, an object of type 'list' is 
created (ChIPSeqSpikeDatasetList or ChIPSeqSpikeDatasetListBoost). Each 
element of the list will contain all experiments corresponding to a given 
input DNA one.

    This function calls different processing steps that overall perform 
ChIP-seq spike-in normalization. The steps and functions are called in the 
following order: Dataset creation (see ?spikeDataset), RPM scaling (see 
?scaling), input DNA subtraction (see ?inputSubtraction), RPM scaling reversal 
(see ?scaling), exogenous DNA scaling (see ?scaling) and binding values 
extraction (see ?extractBinding). 

    For details on installing reference genomes, see details of the function 
'?getPlotSetArray' of the 'seqplots' package.

    For more details on parameters profile_length_before, profile_length_after,
mean_or_median, interpolation_number, interpolation_average and ignore_strand, 
see ?extractBinding.

    If boost = TRUE, either a ?ChIPSeqSpikeDatasetBoost or 
?ChIPSeqSpikeDatasetListBoost object is created. The boost mode enables to 
store the binding values in the form of a GRanges object and avoid reading/
writing files at each processing step. Even if faster, this mode however 
consumes much more memory and should be used with caution.

  If outputFolder is not NULL, the original bigwig files should be copied to 
this folder before performing the analysis. This parameter was created to test
the package with the provided files in extdata/.
}

\value{
    Returns a spike-in normalized object with extracted binding values that can
be used to perform graphical representations (see ?plotProfile, ?plotTransform,
?plotHeatmaps, ?boxplotSpike and ?plotCor).

    According to the files provided in 'infoFile', different objects are 
returned:
\itemize{

\item  A ChIPSeqSpikeDataset object if only one input DNA experiment is 
provided.
\item  A ChIPSeqSpikeDatasetBoost object if only one input DNA experiment is 
    provided and boost = TRUE.
\item   A ChIPSeqSpikeDatasetList object if several input DNA experiments are 
provided.
\item A ChIPSeqSpikeDatasetListBoost object if several input DNA experiments 
are provided and boost = TRUE.
}
}

\examples{
info_file_csv <- system.file("extdata/info.csv", package="ChIPSeqSpike");
bam_path <- system.file("extdata/bam_files", package="ChIPSeqSpike");
bigwig_path <- system.file("extdata/bigwig_files", package="ChIPSeqSpike");
gff_vec <- system.file("extdata/test_coord.gff", package="ChIPSeqSpike");
genome_name <- "hg19";
output_folder <- "test_chipseqspike";
bigwig_files <- system.file("extdata/bigwig_files", 
                            c("H3K79me2_0-filtered.bw",
                              "H3K79me2_100-filtered.bw",
                              "H3K79me2_50-filtered.bw",
                              "input_0-filtered.bw",
                              "input_100-filtered.bw",
                              "input_50-filtered.bw"), package="ChIPSeqSpike");

## Copying example files
dir.create("./test_chipseqspike");
result <- file.copy(bigwig_files, "test_chipseqspike");

csds <- spikePipe(info_file_csv, bam_path, bigwig_path, gff_vec, genome_name, 
              verbose = TRUE, outputFolder = output_folder);

csds2 <- spikePipe(info_file_csv, bam_path, bigwig_path, gff_vec, genome_name, 
              boost = TRUE, verbose = TRUE, outputFolder = output_folder);

unlink("test_chipseqspike/", recursive = TRUE);
is(csds);
is(csds2);
}

\seealso{
  \code{\link{ChIPSeqSpikeDataset}}
  \code{\link{ChIPSeqSpikeDatasetList}}
  \code{\link{ChIPSeqSpikeDatasetBoost}}
  \code{\link{ChIPSeqSpikeDatasetListBoost}}
  \code{\link{spikeDataset}}
  \code{\link{scaling}}
  \code{\link{inputSubtraction}}
  \code{\link{extractBinding}}
  \code{\link{plotProfile}}
  \code{\link{plotTransform}}
  \code{\link{plotHeatmaps}}
  \code{\link{boxplotSpike}}
  \code{\link{plotCor}}
  \code{\link[seqplots]{getPlotSetArray}}
}

\author{Nicolas Descostes}
