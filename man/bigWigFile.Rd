\docType{methods}
\name{bigWigFile}
\alias{bigWigFile}
\alias{bigWigFile<-}
\alias{bigWigFile<-,ChIPSeqSpikeDataset-method}
\alias{bigWigFile<-,ChIPSeqSpikeDatasetBoost-method}
\alias{bigWigFile<-,Experiment-method}
\alias{bigWigFile<-,ExperimentLoaded-method}
\alias{bigWigFile<-,ChIPSeqSpikeCore-method}


\title{Set the path to the bigwig file}

\description{
    Set the path to the input or experiment bigwig file.
}

\usage{
bigWigFile(theObject) <- value
            
            
\S4method{bigWigFile}{ChIPSeqSpikeDataset}(theObject) <- value

\S4method{bigWigFile}{ChIPSeqSpikeDatasetBoost}(theObject) <- value

\S4method{bigWigFile}{Experiment}(theObject) <- value

\S4method{bigWigFile}{ExperimentLoaded}(theObject) <- value

\S4method{bigWigFile}{ChIPSeqSpikeCore}(theObject) <- value
}

\arguments{
    \item{theObject}{A \code{ChIPSeqSpikeCore}, \code{ChIPSeqSpikeDataset}, 
                    \code{ChIPSeqSpikeDatasetBoost}, \code{Experiment} or 
                    \code{ExperimentLoaded} object}
 
    \item{value}{A string representing a path to a bigwig file}
}

\value{
  The modified object is returned
}

\details{
  If the object is ChIPSeqSpikeCore, ChIPSeqSpikeDataset or 
  ChIPSeqSpikeDatasetBoost, bigWigFile will set the path to the input DNA 
  experiment bigwig file corresponding to all experiments defined by the 
  object.
  
  If the object is Experiment or ExperimentLoaded, bigWigFile will set the path
  to the experiment bigwig file.
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

write("toto", file="./test_chipseqspike/toto.bw");
getBigWigFile(csds[[1]][[1]]);
bigWigFile(csds[[1]][[1]]) <- "test_chipseqspike/toto.bw";
getBigWigFile(csds[[1]][[1]]);

unlink("test_chipseqspike/", recursive = TRUE);
}

\author{Nicolas Descostes}
