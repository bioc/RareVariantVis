\name{multipleVis}
\alias{multipleVis}

\title{
Interactive visualization of rare variants on the chromosome, applicable for multiple files
}
\description{
Reads files containing table of rare variants from one chromosome and provides adequate multiple sample visualization. Input files can be obtained from function chromosomeVis.
Function outputs visualization html figure. Figure depicts samples in subfigures. Subfigures illustrate variants (dots) in their genomic coordinates (x axis). Ratio of alternative reads and depth (y axis) gives information about type of variant: homozygous alternative (expected ratio 1) and heterozygous (expected ratio 0.5). Zoom to the figures is possible, by marking the region of interest with mouse left click. Right click induces zoom out and return to the original plot.
Pointing on variants provides basic information about the variant - gene name and position on chromosome.
}

\usage{
multipleVis(inputFiles, outputFile, sampleNames, chromosome)
}

\arguments{
    \item{inputFiles}{Vector of strings containing input file names.}
    \item{outputFile}{Output file name (string).}
    \item{sampleNames}{Vector of sample names (strings).}
    \item{chromosome}{Name of the chromosome to be analyzed (string).}

}

\value{
    \item{comp1}{function returns html visualization file for specified samples}
}

\author{
Adam Gudys and Tomasz Stokowy
}


\examples{
file1 = system.file("extdata", "RareVariants_CoriellIndex_S1.txt",
                    package = "DataRareVariantVis")
file2 = system.file("extdata", "RareVariants_Coriell_S2.txt",
                    package = "DataRareVariantVis")
inputFiles = c(file1, file2)
sampleNames = c("CoriellIndex_S1", "Coriell_S2");
multipleVis(inputFiles, "CorielSamples.html", sampleNames, "19")
}

\keyword{ ~visualization }
\keyword{ ~rare varaints }
