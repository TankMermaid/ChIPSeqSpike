\name{spikeDataset}
\alias{spikeDataset}

\title{spikeDataset constructors function}

\description{
  Function for calling the correct constructor using csv or txt input file
}


\usage{
  spikeDataset(infoFile, bamPath, bigWigPath, boost = FALSE, 
        verbose = TRUE)
}

\arguments{

    \item{infoFile}{csv or txt file containing information about files (see 
                    details).}
    
    \item{bamPath}{Path to the folder containing bam files.}
    
    \item{bigWigPath}{Path to the folder containing bigwig files.}
    
    \item{boost}{If TRUE, the object created enables to perform the analysis in
                 boost mode (see details). Default is FALSE.}
           
    \item{verbose}{If FALSE, do not print messages about object creation. 
                   Default is TRUE}
}


\details{
    'infoFile' should be a csv or a tab separated txt file. The column names 
should be: expName, endogenousBam, exogenousBam, inputBam, bigWigEndogenous and
 bigWigInput. These columns indicate the experiment names; the bam file names 
of data aligned to the reference genome; the bam file names of data aligned to 
the exogenous genome; the input DNA bam file names corresponding to each 
experiment; the bigwig file names of data aligned to the reference genome and 
the bigwig file names of input DNA experiments.

    If 'infoFile' contains only one input file (specified for each experiment),
 a ChIPSeqSpikeDataset (or ChIPSeqSpikeDatasetBoost) object is created. if 
'infoFile' contains different input DNA files, an object of type 'list' is 
created (ChIPSeqSpikeDatasetList or ChIPSeqSpikeDatasetListBoost). Each element
 of the list will contain all experiments corresponding to a given input DNA 
one.

    If boost = TRUE, either a ChIPSeqSpikeDatasetBoost or 
ChIPSeqSpikeDatasetListBoost object is created. The boost mode enables to store
the binding values in the form of a GRanges object and avoid reading/writing 
files at each processing step. Even if faster, this mode however consumes much 
more memory and should be used with caution.
}

\value{
    A ChIPSeqSpikeDataset object if only one input DNA experiment is provided.

    A ChIPSeqSpikeDatasetBoost object if only one input DNA experiment is 
    provided and boost = TRUE.

    A ChIPSeqSpikeDatasetList object if several input DNA experiments are 
provided.
    
    A ChIPSeqSpikeDatasetListBoost object if several input DNA experiments are 
provided and boost = TRUE.
}

\examples{
info_file_csv <- system.file("extdata/info.csv", package="ChIPSeqSpike")
bam_path <- system.file(c("extdata/bam_files"), package="ChIPSeqSpike")
bigwig_path <- system.file(c("extdata/bigwig_files"), package="ChIPSeqSpike")

## Returns ChIPSeqSpikeDatasetList
csds <- spikeDataset(infoFile = info_file_csv, bamPath = bam_path, 
                     bigWigPath = bigwig_path)
is(csds)
csds

## Returns ChIPSeqSpikeDatasetListBoost
csds <- spikeDataset(infoFile = info_file_csv, bamPath = bam_path, 
                     bigWigPath = bigwig_path, boost = TRUE)
is(csds)
csds
}

\seealso{
  \code{\link{ChIPSeqSpikeDataset}}
  \code{\link{ChIPSeqSpikeDatasetBoost}}
  \code{\link{ChIPSeqSpikeDatasetList}}
  \code{\link{ChIPSeqSpikeDatasetListBoost}}
  \code{\link{spikePipe}}
  
}


\author{Nicolas Descostes}
