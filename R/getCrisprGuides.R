.alternate <- function(sequences, alt_allele, ref_width, gsize) {
    paste0(substr(sequences, start = 1, stop = gsize - nchar(alt_allele)),
           alt_allele,
           substring(sequences, first = gsize + ref_width - nchar(alt_allele) + 1))
}

getCrisprGuides <- function(df,
                            genome = BSgenome.Hsapiens.UCSC.hg19::
                                BSgenome.Hsapiens.UCSC.hg19,
                            gsize = 23,
                            PAM = "GG",
                            PAM_rev = "CC") {

    ref_width <- nchar(df$ref_allele)
    alt_width <- nchar(df$alt_allele)

    sequences <- BSgenome::getSeq(
        genome,
        paste0("chr", df$chromosome),
        as.numeric(df$final_positions) - gsize + alt_width,
        as.numeric(df$final_positions) + ref_width + gsize - alt_width - 1,
        as.character = FALSE)

    sequences <- .alternate(sequences, tolower(df$alt_allele), ref_width, gsize)
    PAMs <- gregexpr(pattern = PAM, sequences)
    PAMs <- sapply(PAMs, function(x) x[x > gsize][1])
    guides <- substr(sequences,
                     start = PAMs + nchar(PAM) - gsize,
                     stop = PAMs + nchar(PAM) - 1)

    # for those with NA try the other strand
    PAMs_rev <- gregexpr(pattern = PAM_rev, sequences[is.na(PAMs)])
    PAMs_rev <- sapply(PAMs_rev, function(x) {
        x <- x[x < gsize & x > -1]
        if (length(x) == 0) NA else x[length(x)]
    })

    guides[is.na(PAMs)] <- substr(sequences[is.na(PAMs)],
                                  start = PAMs_rev,
                                  stop = PAMs_rev - nchar(PAM_rev) + gsize + 1)
    guides
}
