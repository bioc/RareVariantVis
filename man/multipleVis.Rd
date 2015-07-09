\name{multipleVis}
\alias{multipleVis}

\title{
Interactive visualization of rare variants on the chromosome, applicable for multiple files
}
\description{
Reads files containing table of rare variants from one chromosome and provides adequate multiple sample visualization. Input files can be obtained from function chromosomeVis.
Function outputs visualization html figure in current working directory. Figure depicts samples in subfigures. Subfigures illustrate variants (dots) in their genomic coordinates (x axis). Ratio of alternative reads and depth (y axis) gives information about type of variant: homozygous alternative (expected ratio 1) and heterozygous (expected ratio 0.5). Zoom to the figures is possible, by marking the region of interest with mouse left click. Right click induces zoom out and return to the original plot.
Pointing on variants provides basic information about the variant - gene name and position on chromosome.
}

\usage{
multipleVis(filesList,filesOrder,sampleNames,chromosome,centromeres,inHouse,inHouseSampNr)
}

\arguments{
  \item{filesList}{a string indicating files names that includes input files}
  \item{filesOrder}{vector including integers specifying order of displayed files. Default is 1:n, where n is number of files}
  \item{sampleNames}{a vector of strings including sample names. Default is names of the files}
  \item{chromosome}{a numeric value giving the analyzed chromosome number}
  \item{centromeres}{a data frame providing location of centromeres on chromosomes}
  \item{inHouse}{a logical value, if TRUE, creating inhouse table in .txt file, if FALSE, case is ignored}
  \item{inHouseSampNr}{a numeric value giving the threshold of variant repetition. Samples above this value are stored in inhouse table}
}

\value{
    \item{comp1}{function returns html visualization file for trio and txt file with
    inheritance annotated rare variants}
}

\author{
Mateusz Garbulowski
}


\examples{
dummyMother=MotherRareVariantsChr21
dummyFather=FatherRareVariantsChr21
filesList=c(MotherRareVariantsChr21, SonRareVariantsChr21,
FatherRareVariantsChr21, dummyMother, dummyFather)
filesOrder=c(1, 2, 3, 4, 5)
sampleNames=c('Mother', 'Son', 'Father', 'dummyMother', 'dummyFather')
multipleVis(filesList, filesOrder, sampleNames, 21, CentromeresHg19, inHouse=TRUE, inHouseSampNr=5)
}

\keyword{ ~visualization }
\keyword{ ~rare varaints }
