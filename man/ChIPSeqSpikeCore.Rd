\name{ChIPSeqSpikeCore-class}
\docType{class}
\alias{ChIPSeqSpikeCore-class}
\alias{ChIPSeqSpikeCore}
\title{ChIPSeqSpikeCore Reference Class}

\description{
 Main class containing input DNA file paths, scores and binding values. The 
 classes ChIPSeqSpikeDataset and ChIPSeqSpikeDatasetBoost inherit from this 
 class. This class is not used directly in the package but is part of 
 ChIPSeqSpikeDataset and ChIPSeqSpikeDatasetBoost. The constructor should not 
 be used directly.
}

\section{Fields}{
  \describe{
    
    inputBam: File path to the input control BAM file
    
    inputBigWig: File path to the input control BigWig file
    
    inputScalingFactor: Input control scaling factor
    
    inputCount: Input control reads count
    
    plotSetArrayList: List of average binding values
    
    matBindingValList: List of binding values matrices
  }
}

\section{Constructor}{
  \describe{
    ChIPSeqSpikeCore(inputBamFile, inputBigWigFile, inputSF = 0, 
        inputNb = 0, SetArrayList = list(), matBindingList = list())
  }
}

\section{Arguments}{
  \describe{
    \item{inputBamFile}{String representing the file path to 
                               the input control BAM file.}

    \item{inputBigWigFile}{String representing the file path to the 
                                  input control bigWig file. (see details)}

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

    \item{\code{getScalingFactor}}{Returns the input scaling factor}

    \item{\code{getCount}}{Returns the number of reads contained in the input 
                           BAM file}

    \item{\code{getAverageBindingValues}}{Returns a list of PlotSetArray 
                                          objects. (see details)}

    \item{\code{getMatBindingValues}}{Returns a list of matrices containing 
                                      binding values. (see details)}
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
  }
}

\section{Details}{'inputSF' is the scaling factor that will be applied to the 
input bigWigFile before input subtraction of the different experiments. 
'inputNb' which holds the number of aligned reads is used to calculate the 
aforementioned factor. Only the endogenous count and factor are needed for the 
input conversely to the experiments for which both the endogenous/exogenous 
scaling factors and counts are needed to perform spike-in normalization.

'SetArrayList' contains PlotSetArray objects. The PlotSetArray class is defined
in the Bioconductor package 'seqplots' and holds the values that are necessary 
to plot profiles and heatmaps. These values can be retrieved with the 
'getAverageBindingValues' function.

'matBindingList' contains matrices of binding values for each experiment. These
values are used to generate boxplots and correlations plots. They are retrieved
by calling the function 'BWGFile_summary' of the bioconductor package 
'rtracklayer'. 
}


\seealso{
  \code{\link{Experiment-class}}
  \code{\link{ExperimentLoaded-class}}
  \code{\link{ChIPSeqSpikeDatasetBoost-class}}
  \code{\link{ChIPSeqSpikeDataset-class}}
  \code{\link{spikeDataset}}
  \code{\link[seqplots]{PlotSetArray}}
}

\author{Nicolas Descostes}
