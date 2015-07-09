\name{chromosomeVis}
\alias{chromosomeVis}

\title{
Visualization of variants on the chromosome
}
\description{
Reads file containing table of variants (or alternatively, vcf file) from one chromosome and provides adequate visualization.
Function outputs visualization png figure in current working directory. Figure illustrates variants (dots) in their genomic coordinates (x axis). Ratio of alternative reads and depth (y axis) gives information about type of variant: homozygous alternative (expected ratio 1) and heterozygous (expected ratio 0.5). Green dots represent rare variants that pass filters: coding/UTR, nonsynonymous variant with dbSNP frequency < 0.01. Orange vertical lines depict position of centromere.
Red curve illustrates moving average of alternative reads/depth ratio. High values of this curve (exceedint 0.75) can suggest potential homozygous/deleterious regions.
In addition, file containing table with rare variants only is generated. Table includes the same columns as input file.
}
\usage{
chromosomeVis(file, sampleName, chromosome, centromeres, pngWidth,
pngHeight, plot, frFilter, dpFilter, vcf)
}

\arguments{
    \item{file}{a connection or a character string giving the
         name of the text file/vcf to load. File should include following information:
         Start.position - for location on chromosome, SNP.Frequency - for dbSNP frequency,
         DP - sequencing depth, AD - allelic depths for reference and alternative alleles,
         Gene.name - gene symbol, Gene.component and Variant.type.
         Example of the table is provided in attached SonVariantsChr21.txt file.
         Example file can be loaded using command data(SonVariantsChr21) from library
         AshkenazimSonChr21.}
    \item{sampleName}{character string giving the
         sample name}
    \item{chromosome}{a numeric value giving the
         analyzed chromosome number}
    \item{centromeres}{a data frame providing location of centromeres on chromosomes}
    \item{pngWidth}{a numeric value setting width of output figure}
    \item{pngHeight}{a numeric value setting height of output figure}
    \item{plot}{logic value indicating if plot should be displayed or saved on the current working
         directory}
    \item{frFilter}{numeric value indicating dbSNP frequency filtering threshold (default 0.01).
         It is recommended not to change this value in the first phase of experiments.}
    \item{dpFilter}{numeric value indicating DP filtering threshold (default 10).
         It is recommended not to change this value in the first phase of experiments.}
    \item{vcf}{logic value indicating if input file is in vcf format. Vcf is supported, however
         data frame inputs are recommended.}
}

\value{
    \item{comp1}{function returns html file with visualization and txt file with rare variants}
}

\author{
Tomasz Stokowy
}

\examples{
## example data
library(AshkenazimSonChr21)

## example of data frame input
chromosomeVis(SonVariantsChr21, "1", 21, CentromeresHg19, 1600, 1200)

## example of vcf file input
fl <- system.file("extdata", "SonVariantsChr21.vcf.gz", package="AshkenazimSonChr21")
chromosomeVis(fl, "1", 21, CentromeresHg19, 1600, 1200, vcf = TRUE)
}

\keyword{ ~visualization }
\keyword{ ~varaints }
