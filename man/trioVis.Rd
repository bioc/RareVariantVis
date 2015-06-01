\name{trioVis}
\alias{trioVis}

\title{
Interactive visualization of rare variants on the chromosome, applicable for trio data
}
\description{
Reads files containing table of rare variants from one chromosome and provides adequate trio visualization. Input files can be obtained from function chromosomeVis.
Function outputs visualization html figure in current working directory. Figure depicts three samples in subfigures: mother, index and father. Subfigures illustrate variants (dots) in their genomic coordinates (x axis). Ratio of alternative reads and depth (y axis) gives information about type of variant: homozygous alternative (expected ratio 1) and heterozygous (expected ratio 0.5). Zoom to the figures is possible, by marking the region of interest with mouse left click. Right click induces zoom out and return to the original plot.
Pointing on variants provides basic information about the variant - gene name and position on chromosome.

}
\usage{
trioVis(fileMother, fileIndex, fileFather, sampleMother,
    sampleIndex, sampleFather, chromosome, centromeres)
}

\arguments{
        \item{fileMother}{a connection or a character string giving the
         name of the file with mother variants to load. File should include tab delimited text file
         with columns:
         Start.position - for location on chromosome, SNP.Frequency - for dbSNP frequency,
         DP - sequencing depth, AD - allelic depths for reference and alternative alleles,
         Gene.name - gene symbol.}
        \item{fileIndex}{a connection or a character string giving the
         name of the file with mother variants to load. File should include tab delimited text file
         with columns:
         Start.position - for location on chromosome, SNP.Frequency - for dbSNP frequency,
         DP - sequencing depth, AD - allelic depths for reference and alternative alleles,
         Gene.name - gene symbol.}
        \item{fileFather}{a connection or a character string giving the
         name of the file with mother variants to load. File should include tab delimited text file
         with columns:
         Start.position - for location on chromosome, SNP.Frequency - for dbSNP frequency,
         DP - sequencing depth, AD - allelic depths for reference and alternative alleles,
         Gene.name - gene symbol.}
        \item{sampleMother}{character string giving the mother sample name}
        \item{sampleIndex}{character string giving the index sample name}
        \item{sampleFather}{character string giving the father sample name}
        \item{chromosome}{a numeric value giving the
         analyzed chromosome number}
        \item{centromeres}{a data frame providing location of centromeres on chromosomes}
}

\value{
    \item{comp1}{function returns html visualization file for trio and txt file with
    inheritance annotated rare variants}
}

\author{
Tomasz Stokowy
}

\examples{
trioVis(MotherRareVariantsChr21, SonRareVariantsChr21, FatherRareVariantsChr21, "Mother", "Son",
"Father", 21, CentromeresHg19)
}

\keyword{ ~visualization }
\keyword{ ~rare varaints }
