\name{ChIPSeqSpikeDataset-class}
\docType{class}
\alias{ChIPSeqSpikeDataset-class}
\alias{ChIPSeqSpikeDataset}
\alias{[[,ChIPSeqSpikeDataset,ANY,ANY-method}
\alias{[[<-,ChIPSeqSpikeDataset,ANY,ANY-method}
\title{ChIPSeqSpikeDataset Reference Class}

\description{
 Main class containing file paths, scores and values for spike normalization. 
 It inherits from ChIPSeqSpikeCore.
}

\section{Fields}{
  \describe{
    experimentList: A list of Experiment objects
  }
}

\section{Constructor}{
  \describe{
    ChIPSeqSpikeDataset(endogenousBam_vec, exogenousBam_vec, 
        bigWigFile_endogenous_vec, inputBigWigFile, inputBamFile, expnames,
        inputSF = 0, inputNb = 0, SetArrayList = list(), 
        matBindingList = list())
  }
}

\section{Arguments}{
  \describe{
    \item{endogenousBam_vec}{Character vector of file paths to the BAM 
                                   files aligned to the reference genome.}

    \item{exogenousBam_vec}{Character vector of file paths to the BAM 
                                   files aligned to the exogenous genome.}

    \item{bigWigFile_endogenous_vec}{Character vector of file paths to 
                             the bigWig files aligned to the reference genome.}

    \item{inputBigWigFile}{String representing the file path to the input
                                  control bigWig file. (see details)}

    \item{inputBamFile}{String representing the file path to 
                               the input control BAM file.}

    \item{expnames}{Character vector of experiment names. (see details)}

    \item{inputSF}{Numeric scaling factor. Default is 0. (see details)}

    \item{inputNb}{Numeric read counts. Default is 0. (see details)}

    \item{SetArrayList}{List of PlotSetArray objects. 
                               Default is an empty list. (see details)}

    \item{matBindingList}{List of binding value matrices.
                                 Default is an empty list. (see details)}
  }
}

\section{Getters}{
  \describe{
    \item{\code{getBam}}{Returns the input BAM path}

    \item{\code{getBigWigFile}}{Returns the input bigWig path}

    \item{\code{getExperimentListBigWigs}}{Returns a character vector of paths 
                                         to the experiment bigWig files}

    \item{\code{getExpName}}{Returns a character vector of experiment names}

    \item{\code{getScalingFactor}}{Returns the input scaling factor}

    \item{\code{getCount}}{Returns the number of reads contained in the input 
                           BAM file}

    \item{\code{getAverageBindingValues}}{Returns a list of PlotSetArray 
                                          objects. (see details)}

    \item{\code{getMatBindingValues}}{Returns a list of matrices containing 
                                      binding values. (see details)}

    \item{\code{x[[i]]}}{Get the Experiment object at index i}
  }
}

\section{Setters}{
  \describe{
    \item{\code{scalingFactor}}{Modifies the input scaling factor value}

    \item{\code{count}}{Modifies the input count value}

    \item{\code{bigWigFile}}{Modifies the input bigWig file path}

    \item{\code{averageBindingValues}}{Modifies the PlotSetArray list. 
                                      (see details)}

    \item{\code{matBindingValues}}{Modifies the list of binding value matrices.
                                   (see details)}

    \item{\code{x[[i]] <- value}}{Set \code{value} to ChIPSeqSpikeDataset i}
  }
}

\section{Details}{'expnames' character vector is used to define the names of 
the experiment list and are used as labels in plotting, summary and getRatio 
functions.

'inputSF' is the scaling factor that will be applied to the input bigWigFile
before input subtraction of the different experiments. 'inputNb' which holds 
the number of aligned reads is used to calculate the aforementioned factor. 
Only the endogenous count and factor are needed for the input conversely to the
 experiments for which both the endogenous/exogenous scaling factors and counts
 are needed to perform spike-in normalization.

'SetArrayList' contains PlotSetArray objects. The PlotSetArray class is defined
in the Bioconductor package 'seqplots' and holds the values that are necessary 
to plot profiles and heatmaps. These values can be retrieved with the 
'getAverageBindingValues' function.

'matBindingList' contains matrices of binding values for each experiment. These
values are used to generate boxplots and correlations plots. They are retrieved
by calling the function 'BWGFile_summary' of the bioconductor package 
'rtracklayer'. 

If the dataset contains more than one input, one would want to use the 
ChIPSeqSpikeDatasetList class. Boost mode classes (ChIPSeqSpikeDatasetBoost and
ChIPSeqSpikeDatasetListBoost) can also be considered to speed up the analysis. 
}


\examples{

file_vec <- system.file("extdata", 
                        c("bam_files/H3K79me2_0_dm3-filtered.bam", 
                        "bam_files/H3K79me2_0_hg19-filtered.bam", 
                        "bigwig_files/H3K79me2_0-filtered.bw", 
                        "bigwig_files/input_0-filtered.bw", 
                        "bam_files/input_0_hg19-filtered.bam"),
                        package="ChIPSeqSpike")

csds <- ChIPSeqSpikeDataset(endogenousBam_vec = file_vec[2],
                         exogenousBam_vec = file_vec[1],
                         bigWigFile_endogenous_vec = file_vec[3], 
                         inputBigWigFile = file_vec[4], 
                         inputBamFile = file_vec[5], 
                         expnames = "H3K79me2_0")
csds
}

\seealso{
  \code{\link{Experiment-class}}
  \code{\link{ChIPSeqSpikeCore-class}}
  \code{\link{ChIPSeqSpikeDatasetBoost-class}}
  \code{\link{ChIPSeqSpikeDatasetList-class}}
  \code{\link{spikeDataset}}
  \code{\link[seqplots]{PlotSetArray}}
  \code{\link{spikeSummary}}
  \code{\link{getRatio}}
}

\author{Nicolas Descostes}
