\name{callHomozygous}
\alias{callHomozygous}
\title{
Call homozygous regions from sequencing data
}
\description{
Function calls homozygous regions from whole genome sequencing data.
}
\usage{
callHomozygous(sample, chromosomes, caller = "speedseq", MA_Window = 2000, HMZ_length  = 100000, min_n_HMZ = 20)
}

\arguments{
  \item{sample}{A name of SNV sample file to be analyzed.}
  \item{chromosomes}{A vector of strings indicating chromosomes to be analyzed.}
  \item{caller}{A string indicating vcf caller. Default is "speedseq", supports "GATK"}
    \item{MA_Window}{A number indicating window size for moving average function. Recommended value for genome is 2000, for exome is 20. Default is 2000.}
  \item{HMZ_length}{Minimal length of homozygous region to be called. Default is 100000.}
  \item{min_n_HMZ}{Minimal number of variants necessary to call a region. Default is 20.}
}

\value{
    \item{comp1}{function calls homozygous regions from whole genome sequencing data and returns them in a tab separated txt file.}
}

\author{
Tomasz Stokowy
}

\examples{
sample = system.file("extdata", "CoriellIndex_S1.vcf.gz",
                     package = "DataRareVariantVis")
callHomozygous(sample=sample, chromosomes=c("19"))
}

\keyword{ ~call homozygous regions }

