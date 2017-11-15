\docType{data}
\name{ratio}
\alias{ratio}

\title{Result of method getRatio on the complete dataset}

\description{
  Result of method getRatio on the complete dataset.
}

\usage{
  data(ratio)
}

\format{matrix}

\details{
## Complete Dataset

The data used in this documentation represent a gold-standard example of the 
importance of using spike-in controls with ChIP-Seq experiments. It uses 
Drosophila Melanogaster chromatin as exogenous spike-in control to correct 
experimental biases. Without spike-in control and using only RPM normalization,
proper differences of H3K79me2 histone modification in human Jurkat cells upon 
EPZ5676 inhibitor treatment are not observed [1]. 

This dataset is made of bigwig and bam files of H3K79me2 ChIP-Seq data and 
corresponding input DNA controls.
Bam files contain data aligned to the Human reference genome Hg19 or to the 
Drosophila reference genome dm3. The latest is used to compute external 
spike-in scaling factors. All above mentioned data are available at 0, 50 and 
100 percent EPZ5676 inhibitor treatment (see vignette for data references).
}

\references{
  [1] Orlando et al, "Quantitative ChIP-Seq normalization reveals global 
  modulation of the epigenome", Cell Rep, 2014.
}
\keyword{datasets}