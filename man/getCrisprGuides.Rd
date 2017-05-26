\name{getCrisprGuides}
\alias{getCrisprGuides}
\title{
Retrieve CRISPR/Cas9 guides.
}
\description{
Function checks whether a guideRNA can be found that overlaps given SNP. Returns
sequence of the guideRNA with the variant marked with the lowercase letters.
When multiple guideRNAs are possible for given SNP, guideRNA with the variant
closest to the PAM site is being selected.
}

\arguments{
  \item{df}{ A data frame, preferably out from chromosomeVis.}
  \item{genome}{A object of the BSGenome, by default
  BSgenome.Hsapiens.UCSC.hg19.}
  \item{gsize}{Prefered size of the guideRNA, by default, standard 23 is used.}
  \item{PAM}{Prefered Protospacer Adjecent Motif "PAM", short motif that has to
  be found on the 5' end of the guideRNA, by default, Cas9 is used "GG".}
  \item{PAM_rev}{Prefered Protospacer Adjecent Motif "PAM", short motif that has
  to be found on the reverse strand, by default, Cas9 is used "CC".
  This is checked only when no guideRNA is found on the forward strand.}
}

\value{
    \item{character vector}{Vector of guideRNAs, when no guideRNA was found
    for the forward strand, reverse strand is checked, when no guideRNA is found
    NA is returned.}
}

\author{
Kornel Labun
}

\examples{
file <- system.file("extdata", "RareVariants_CoriellIndex_S1.txt",
    package = "RareVariantVis")
df <- read.delim(file, stringsAsFactors = FALSE)
getCrisprGuides(df)
}

\keyword{ ~CRISPR/Cas9}

