\name{ChIPSeqSpikeDatasetListBoost-class}
\docType{class}
\alias{ChIPSeqSpikeDatasetListBoost-class}
\alias{ChIPSeqSpikeDatasetListBoost}
\alias{[[,ChIPSeqSpikeDatasetListBoost,ANY,ANY-method}
\title{ChIPSeqSpikeDatasetListBoost Reference Class}

\description{
 Class containing a list of ChIPSeqSpikeDatasetBoost objects
}

\section{Fields}{
  \describe{
    datasetList: A list of ChIPSeqSpikeDatasetBoost objects
    }
}

\section{Constructor}{
  \describe{
    ChIPSeqSpikeDatasetListBoost(dataset_list, verbose)
  }
}

\section{Arguments}{
  \describe{
    \item{dataset_list}{A properly formatted list of information needed 
                               to create the object (see details).}

    \item{verbose}{Indicate if processing messages should be output.}
  }
}

\section{Getters}{
  \describe{
    \item{\code{exportBigWigs}}{Output all bigwig files corresponding to the 
                                previously performed transformations. The list
                                of bigwig files is given in info.csv (see 
                                details).}

    \item{\code{x[[i]]}}{Get the ChIPSeqSpikeDatasetBoost object at index i}
  }
}

\section{Details}{
This class enables to process datasets containing different input files in 
boost mode. It will creates a list of ChIPSeqSpikeDatasetBoost objects, each of
 them containing a different input.

The above indicated constructor should not be used directly. One would rather
create the object by calling the "meta-constructor" spikeDataset function, 
which takes a info.csv file as input, setting boost = TRUE. spikeDataset 
function formats properly the different information into a dataset_list which 
is submitted to the ChIPSeqSpikeDatasetListBoost constructor.
}


\examples{
info_file_csv <- system.file("extdata/info.csv", package="ChIPSeqSpike")
bam_path <- system.file(c("extdata/bam_files"), package="ChIPSeqSpike")
bigwig_path <- system.file(c("extdata/bigwig_files"), package="ChIPSeqSpike")

csds <- spikeDataset(infoFile = info_file_csv, bamPath = bam_path, 
                     bigWigPath = bigwig_path, boost = TRUE)
                     
is(csds)
csds
}

\seealso{
  \code{\link{ChIPSeqSpikeDatasetBoost-class}}
  \code{\link{spikeDataset}}
  \code{\link{spikeSummary}}
  \code{\link{getRatio}}
}

\author{Nicolas Descostes}
