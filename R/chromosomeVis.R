chromosomeVis <- function(
    sample,
    sv_sample,
    dbSNP_file,
    Exac_file,
    chromosomes,
    pngWidth = 1600,
    pngHeight = 1200,
    caller = "speedseq",
    MA_Window = 1000,
    coding_regions_file = system.file("extdata", "nexterarapidcapture_exome_targetedregions_v1.2_chr19_9-10.bed",
       package = "RareVariantVis"),
    annotation_file = system.file("extdata", "UCSC_hg19_chr9_9-10_refSeq_160702.txt",
       package = "RareVariantVis"),
    uniprot_file = system.file("extdata", "uniprot-all_chr19_9-10.txt",
       package = "RareVariantVis")
    ) {


    sampleName = unlist(strsplit(basename(sample), '[.]'))[1]

    if (missing(dbSNP_file)) {
        dbSNP_file = system.file("extdata", "All_20160601_chr19_9-10.vcf.recode.vcf.gz", package =
                                          "RareVariantVis")
    }

    if (missing(Exac_file)) {
        Exac_file = system.file("extdata", "ExAC.r0.3.1.sites.vep_chr19_9-10.vcf.recode.vcf.gz", package =
                                          "RareVariantVis")
    }

    ###############################################
    # Part 0 - prepare adnnotations
    ###############################################

    # Centromeres for hg19

    url = "http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/cytoBand.txt.gz"
    cytobands = read.delim(
        text=readLines(gzcon(url(url))),
        header=FALSE, sep="\t", stringsAsFactors = FALSE
    )

    cytobands_acen = cytobands[cytobands$V5 == "acen",]
    cytobands_chr = cytobands_acen[1:dim(cytobands_acen)[1] %% 2 == 1,]
    cytobands_start = cytobands_acen[1:dim(cytobands_acen)[1] %% 2 == 1,]$V2
    cytobands_end =  cytobands_acen[1:dim(cytobands_acen)[1] %% 2 == 0,]$V3
    cytobands_summary = cbind.data.frame(as.character(cytobands_chr$V1),
                                         cytobands_start, cytobands_end, stringsAsFactors = FALSE)

    colnames(cytobands_summary) = c("Chrom", "Start", "Stop")

    url = "http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/chromInfo.txt.gz"
    chromosome_lengths = read.delim(
        text=readLines(gzcon(url(url))),
        header=FALSE, sep="\t", stringsAsFactors = FALSE
    )

    chr_lengths = cbind(chromosome_lengths[1:24,1:2],
                        ceiling(chromosome_lengths[1:24,2] / 10000000) * 10000000)
    colnames(chr_lengths) = c("Chrom", "ChromLen", "WinSize")

    temp_centromeres = merge(cytobands_summary, chr_lengths, by = "Chrom")

    order_rows_centromeres = match(c(1:22, "X", "Y"), substr(temp_centromeres$Chrom, 4,
                                                             nchar(as.character(temp_centromeres$Chrom))))

    centromeres = temp_centromeres[order_rows_centromeres,]

    # End of centromeres preparation



    # Annotation files
    # centromeres_file    = system.file("extdata", "CentromeresHg19.txt", package =
    #                                          "RareVariantVis")
    # coding_regions_file = system.file("extdata",
    #                                  "nexterarapidcapture_exome_targetedregions_v1.2_chr19_9-10.bed",
    #                                  package = "RareVariantVis")
    # annotation_file     = system.file("extdata", "UCSC_hg19_chr9_9-10_refSeq_160702.txt", package =
    #                                      "RareVariantVis")
    # uniprot_file        = system.file("extdata", "uniprot-all_chr19_9-10.txt", package =
    #                                      "RareVariantVis")
    # load annotations
    # centromeres = read.delim(centromeres_file,
    #                             header = TRUE,
    #                             stringsAsFactors = FALSE)
    anno_whole = read.delim(annotation_file,
                            header = TRUE,
                            stringsAsFactors = FALSE)
    uniprot = read.delim(uniprot_file,
                         header = TRUE,
                         stringsAsFactors = FALSE)
    exome = read.delim(coding_regions_file,
                       header = FALSE,
                       stringsAsFactors = FALSE)

    # create table with empty row (to make binding work properly)
    final_table = data.frame(
        chromosome = "",
        final_positions = NA,
        variant_type = "",
        ref_allele = "",
        alt_allele = "",
        final_variations = NA,
        conservation = NA,
        gene_name = "",
        in_exon = NA,
        stringsAsFactors = FALSE
    )

    final_table = cbind(final_table, uniprot[1, ])

    sv_table = data.frame(
        chromosome = "",
        start_position = NA,
        ID = "",
        REF = "",
        ALT = "",
        QUAL = NA,
        SVTYPE = "",
        GT = "",
        end_position = NA,
        gene_name = "",
        in_exon = NA,
        stringsAsFactors = FALSE
    )

    sv_table = cbind(sv_table, uniprot[1, ])

    complex_table = data.frame(
        chromosome = "",
        start_position = NA,
        gene_name = "",
        in_exon = NA,
        stringsAsFactors = FALSE
    )

    complex_table = cbind(complex_table, uniprot[1, ])

    # begin chromosome loop
    for (chromosome in chromosomes) {
        # determine chromosome index
        if (chromosome == "X") {
            chromosomeIndex = 23
        } else if (chromosome == "Y") {
            chromosomeIndex = 24
        } else {
            chromosomeIndex = as.integer(chromosome)
        }

        ###############################################
        # Part 1 - read 1 chromosome from a genome
        ###############################################

        param = ScanVcfParam(which = GRanges(chromosome, IRanges(1, 250000000)))
        vcf = readVcf(sample, genome = "hg19", param = param)

        if (caller == "speedseq"){
            firstADs = sapply(geno(vcf)[["AO"]], function(m) m[1])
            secondADs = sapply(geno(vcf)[["RO"]], function(m) m[1])
            positions = start(rowRanges(vcf))
            variation = firstADs / (firstADs + secondADs)
        } else {
            firstADs = sapply( geno(vcf)[["AD"]], function(m) m[1] )
            secondADs = sapply( geno(vcf)[["AD"]], function(n) n[2] )
            positions = start(rowRanges(vcf))
            variation = secondADs / (firstADs + secondADs)
        }

        ###############################################
        # Part 2 - Select coding variants
        ###############################################

        # check for proper chromosome
        chr_select = which(exome[, 1] == paste0("chr", chromosome))

        # binary interval search
        lo = exome$V2[chr_select]

        hi = exome$V3[chr_select]

        # find overlapping ranges
        overlapping = vector(mode = "numeric", length = 0)
        if (length(lo) >= 2) {
            for (i in 2:length(lo)) {
                if (hi[i - 1] >= lo[i]) {
                    overlapping = c(overlapping, i)

                }
            }
        }

        # remove boundaries corresponding to overlapping ranges
        if (length(overlapping) > 0) {
            lo = lo[-overlapping] - 0.1 	# 0.1 is to handle properly boundary values
            hi = hi[-(overlapping - 1)] + 0.1
        }

        coding_ranges = c(rbind(lo, hi))
        # interleaved vector of low and high exon borders

        coding_sel = findInterval(positions, coding_ranges)

        # odd indices correspond to positions inside intervals
        coding_sel = c(0, which(coding_sel %% 2 == 1))

        coding_pos = positions[coding_sel]
        coding_var = variation[coding_sel]

        ###############################################
        # Part 3 - Select rare dbSNP variants
        ###############################################

        # beware of strange CAF format!
        param2 = ScanVcfParam(which = GRanges(
            chromosome,
            IRanges(start = coding_pos, end = coding_pos)
        ))
        vcf2 = readVcf(dbSNP_file, genome = "hg19", param = param2)

        # binary interval search
        lo_vcf2 = start(rowRanges(vcf2))
        hi_vcf2 = end(rowRanges(vcf2))

        # find overlapping ranges to identify complex variants
        overlapping_vcf2 = vector(mode = "numeric", length = 0)
        if (length(lo_vcf2) >= 2) {
            for (i in 2:length(lo_vcf2)) {
                if (hi_vcf2[i - 1] >= lo_vcf2[i]) {
                    overlapping_vcf2 = c(overlapping_vcf2, i)

                }
            }
        }

        # indexes of starting records for new coding positions
        first_codpos_ind = setdiff(1:length(lo_vcf2), overlapping_vcf2)

        # simple and complex variants in dbSNP
        prev_present_ind = first_codpos_ind - 1
        simple_variants_dbsnp_ind = intersect(first_codpos_ind, prev_present_ind)

        # positions of simple and complex variants
        simple_variant_positions = lo_vcf2[simple_variants_dbsnp_ind]
        complex_variant_positions = coding_pos[which(is.na(match(
            coding_pos, simple_variant_positions
        )))]

        # handle allelic frequencies
        CAF_array = matrix(NA,
                           nrow = length(simple_variant_positions),
                           ncol = 30)
        rareDBsnp = rep(TRUE, length(simple_variant_positions))

        # check if we found anything
        if (length(simple_variant_positions) > 0) {
            for (i in 1:length(simple_variant_positions)) {
                currentCAF = unlist(info(vcf2)[["CAF"]][simple_variants_dbsnp_ind][i])

                if (length(currentCAF) > 0) {
                    CAF_array[i, 1:length(currentCAF)] = currentCAF
                }

                if (any(currentCAF == ".")) {
                    rareDBsnp[i] = TRUE
                } else {
                    rareDBsnp[i] = !all(as.numeric(currentCAF) > 0.01)
                }
            }

            rareDBsnp_with_complex_pos = sort(c(
                simple_variant_positions[rareDBsnp],
                complex_variant_positions
            ))
            rareDBsnp_with_complex_var = variation[match(rareDBsnp_with_complex_pos, positions)]

            ###############################################
            # Part 4 - Select rare Exac variants
            ###############################################

            # beware of strange AF format!
            param3 = ScanVcfParam(which = GRanges(
                chromosome,
                IRanges(start = rareDBsnp_with_complex_pos, end = rareDBsnp_with_complex_pos)
            ))
            vcf3 = readVcf(Exac_file, genome = "hg19", param = param3)

            # binary interval search
            lo_vcf3 = start(rowRanges(vcf3))
            hi_vcf3 = end(rowRanges(vcf3))

            # find overlapping ranges to identify complex variants
            overlapping_vcf3 = vector(mode = "numeric", length = 0)
            if (length(lo_vcf3) >= 2) {
                for (i in 2:length(lo_vcf3)) {
                    if (hi_vcf3[i - 1] >= lo_vcf3[i]) {
                        overlapping_vcf3 = c(overlapping_vcf3, i)

                    }
                }
            }

            # indexes of starting records for new coding positions
            first_codpos_ind_exac = setdiff(1:length(lo_vcf3), overlapping_vcf3)

            # simple and complex variants in dbSNP
            prev_present_ind_exac = first_codpos_ind_exac - 1
            simple_variants_dbsnp_ind_exac = intersect(first_codpos_ind_exac, prev_present_ind_exac)

            # positions of simple and complex variants
            simple_variant_positions_exac = lo_vcf3[simple_variants_dbsnp_ind_exac]
            complex_variant_positions_exac = rareDBsnp_with_complex_pos[which(is.na(
                match(
                    rareDBsnp_with_complex_pos,
                    simple_variant_positions_exac
                )
            ))]

            # handle allelic frequencies
            Exac_AF_array = matrix(NA,
                                   nrow = length(simple_variant_positions_exac),
                                   ncol = 30)
            rareExac = rep(TRUE, length(simple_variant_positions_exac))
            multipleAltVariant = rep(FALSE, length(simple_variant_positions_exac))

            for (i in 1:length(simple_variant_positions_exac)) {
                currentAF_exac = unlist(info(vcf3)[["AF"]][simple_variants_dbsnp_ind_exac][i])

                if (length(currentAF_exac) > 1) {
                    rareExac[i] = TRUE
                    multipleAltVariant[i] = TRUE
                }
                else if (length(currentAF_exac) < 1) {
                    rareExac[i] = TRUE
                    multipleAltVariant[i] = FALSE
                }
                else{
                    rareExac[i] = currentAF_exac < 0.01
                    multipleAltVariant[i] = FALSE
                }
            }

            # Exac filtering Summary
            # one variant in position and one alternative allele
            simple_nexac_positions = simple_variant_positions_exac[rareExac &
                                                                       !multipleAltVariant]
            simple_nexac_variation = variation[match(simple_nexac_positions, positions)]

            # one variant in position and many alternative alleles
            complex1_exac_positions = simple_variant_positions_exac[rareExac &
                                                                        multipleAltVariant]
            # many variants in one position
            complex2_exac_positions = complex_variant_positions_exac

            complex_positions = sort(c(complex1_exac_positions, complex2_exac_positions))

            ###############################################
            # Part 5 - Predict non-synonymous
            ###############################################

            txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene

            if(length(simple_nexac_positions) > 0){
                param4 = ScanVcfParam(which = GRanges(
                    chromosome,
                    IRanges(start = simple_nexac_positions, end = simple_nexac_positions)
                ))
                vcf4 = readVcf(sample, genome = "hg19", param = param4)

                if(caller == "GATK"){
                    seqlevels(vcf4) = chromosome
                }

                vcf4 <- renameSeqlevels(vcf4, paste0("chr", chromosome))
                synonymous <- predictCoding(vcf4, txdb, Hsapiens)

                UniqueNonsyn = unique(synonymous[mcols(synonymous)$CONSEQUENCE == "nonsynonymous"])
                UniqueFrames = unique(synonymous[mcols(synonymous)$CONSEQUENCE == "frameshift"])
                UniqueNonsen = unique(synonymous[mcols(synonymous)$CONSEQUENCE == "nonsense"])

                RareNonsynCoding = c(UniqueNonsyn, UniqueFrames, UniqueNonsen)

                final_positions = sort(start(ranges(RareNonsynCoding)))
                final_variations = variation[match(final_positions, positions)]
            } else {
                final_positions = NULL
                final_variations = NULL
            }

            ###############################################
            # Part 6 - Annotate selected and report
            ###############################################

            # filter annotations by chromosome,
            anno = anno_whole[anno_whole$chrom == paste0("chr", chromosome),]

            # find genes containing analysed positions
            ordering_start = order(anno$txStart)
            ordering_end = order(anno$txEnd)
            intervals_start = findInterval(final_positions, anno$txStart[ordering_start])

            intervals_end = findInterval(final_positions, anno$txEnd[ordering_end])

            if (length(final_positions) > 0) {
                # get conservation scores
                conservations =  scores(phastCons100way.UCSC.hg19, scores.only=TRUE,
                                        GRanges(
                                            seqnames = paste0("chr", chromosome) ,
                                            IRanges(start = final_positions, width = 1)
                                        ))

                # analyse all positions
                for (i in 1:length(final_positions)) {
                    row_left = list(
                        chromosome,
                        final_positions[i],
                        as.character(mcols(RareNonsynCoding)$CONSEQUENCE[match(final_positions[i],
                                                                               start(RareNonsynCoding))]),
                        as.character(mcols(RareNonsynCoding)$REF[match(final_positions[i],
                                                                       start(RareNonsynCoding))]),
                        as.character(unlist(
                            mcols(RareNonsynCoding)$ALT
                        ))[match(final_positions[i],
                                 start(RareNonsynCoding))],
                        final_variations[i],
                        conservations[i]
                    )

                    starting_before = ordering_start[1:intervals_start[i]]
                    ending_after = ordering_end[(intervals_end[i] + 1):length(ordering_end)]
                    common = intersect(starting_before, ending_after)

                    # proceed only when position intersects with some gene
                    if (length(common) > 0) {
                        uniqueMask = !duplicated(anno[common,]$name2)

                        uniqueCommon = common[uniqueMask]

                        # for all genes containing the position
                        for (ig in 1:length(uniqueCommon)) {
                            # extract gene name
                            gene_id = uniqueCommon[ig]

                            gene = anno$name2[gene_id]

                            exon_starts = as.numeric(unlist(
                                strsplit(anno$exonStarts[gene_id], split = ",")
                            ))
                            exon_ends = as.numeric(unlist(
                                strsplit(anno$exonEnds[gene_id], split = ",")
                            ))

                            gene_ranges = c(rbind(exon_starts, exon_ends)) # interleaved vector of low and high exon borders
                            interval = findInterval(final_positions[i], gene_ranges)

                            exon = (interval %% 2 == 1)
                            uniprot_row = uniprot[which(uniprot$Gene.names...primary.. == anno[uniqueCommon[ig], ]$name2),]

                            if (dim(uniprot_row)[1] == 0) {
                                uniprot_row = rep(NA, 11)
                            }

                            row = c(row_left, gene, exon, uniprot_row)
                            final_table = rbind(final_table, unlist(row))
                        }
                    } else {
                        row = c(row_left, rep(NA, dim(final_table)[2] - 2))
                        final_table = rbind(final_table, row)
                    }

                } # end position loop


                ###############################################
                # Part 6b - complex variants
                ###############################################

                # # find genes containing analysed positions
                # intervals_start = findInterval(complex_positions, anno$txStart[ordering_start])
                # intervals_end = findInterval(complex_positions, anno$txEnd[ordering_end])
                #
                # if (length(complex_positions) > 0) {
                #
                #     # analyse all positions
                #     for (i in 1:length(complex_positions)) {
                #         row_left = list(
                #             chromosome,
                #             complex_positions[i]
                #         )
                #
                #         starting_before = ordering_start[1:intervals_start[i]]
                #         ending_after = ordering_end[(intervals_end[i] + 1):length(ordering_end)]
                #         common = intersect(starting_before, ending_after)
                #
                #         # proceed only when position intersects with some gene
                #         if (length(common) > 0) {
                #             uniqueMask = !duplicated(anno[common,]$name2)
                #
                #             uniqueCommon = common[uniqueMask]
                #
                #             # for all genes containing the position
                #             for (ig in 1:length(uniqueCommon)) {
                #                 # extract gene name
                #                 gene_id = uniqueCommon[ig]
                #
                #                 gene = anno$name2[gene_id]
                #
                #                 exon_starts = as.numeric(unlist(
                #                     strsplit(anno$exonStarts[gene_id], split = ",")
                #                 ))
                #                 exon_ends = as.numeric(unlist(
                #                     strsplit(anno$exonEnds[gene_id], split = ",")
                #                 ))
                #
                #                 gene_ranges = c(rbind(exon_starts, exon_ends)) # interleaved vector of low and high exon borders
                #                 interval = findInterval(complex_positions[i], gene_ranges)
                #
                #                 exon = (interval %% 2 == 1)
                #                 uniprot_row = uniprot[which(uniprot$Gene.names...primary.. == anno[uniqueCommon[ig], ]$name2),]
                #
                #                 if (dim(uniprot_row)[1] == 0) {
                #                     uniprot_row = rep(NA, 11)
                #                 }
                #
                #                 row = c(row_left, gene, exon, uniprot_row)
                #                 complex_table = rbind(complex_table, unlist(row))
                #             }
                #         } else {
                #             row = c(row_left, rep(NA, 13))
                #             complex_table = rbind(complex_table, row)
                #         }
                #
                #     } # end position loop
                # }
            }
        }

        ###############################################
        # Part 7 - Load and filter structural variants
        ###############################################

        if (!missing(sv_sample)) {

            sv_local_table = read.table(
                text = "",
                col.names = colnames(sv_table),
                colClasses = sapply(sv_table, class)
            )

            param = ScanVcfParam(which = GRanges(chromosome, IRanges(1, 250000000)))
            vcf_sv = readVcf(sv_sample, genome = "hg19", param = param)

            dupdel_variant_mask = (info(vcf_sv)[["SVTYPE"]] == "DUP") |
                (info(vcf_sv)[["SVTYPE"]] == "DEL")
            point_variant_mask = info(vcf_sv)[["SVTYPE"]] == "BND"

            dupdel_variant_positions_all = start(rowRanges(vcf_sv[dupdel_variant_mask]))
            dupdel_variant_ends_all = info(vcf_sv[dupdel_variant_mask])[["END"]]

            point_variant_positions_all = start(rowRanges(vcf_sv[point_variant_mask]))

            fitv_begins = findInterval(dupdel_variant_positions_all, coding_ranges)
            fitv_ends = findInterval(dupdel_variant_ends_all, coding_ranges)

            dupdel_variant_sel = which((fitv_begins %% 2 == 1) |
                                           # variant begins inside coding region
                                           (fitv_ends %% 2 == 1) | # variant ends inside coding region
                                           (fitv_ends - fitv_begins >= 1)) # variant spans over multiple intervals

            dupdel_variant_positions = dupdel_variant_positions_all[dupdel_variant_sel]
            dupdel_variant_ends = dupdel_variant_ends_all[dupdel_variant_sel]

            point_variant_sel = which(findInterval(point_variant_positions_all, coding_ranges) %% 2 == 1)
            point_variant_positions = point_variant_positions_all[point_variant_sel]

            matched_dupdel = match(dupdel_variant_positions, start(rowRanges(vcf_sv)))
            matched_point = match(point_variant_positions, start(rowRanges(vcf_sv)))

            ###############################################
            # Part 8 - Annotate duplications/deletions
            ###############################################

            # find genes containing analysed positions
            ordering_start = order(anno$txStart)
            ordering_end = order(anno$txEnd)

            if (length(dupdel_variant_positions)) {
                variantStart.geneStart = findInterval(dupdel_variant_positions, anno$txStart[ordering_start])

                variantStart.geneEnd = findInterval(dupdel_variant_positions, anno$txEnd[ordering_end])


                variantEnd.geneStart = findInterval(dupdel_variant_ends, anno$txStart[ordering_start])

                variantEnd.geneEnd = findInterval(dupdel_variant_ends, anno$txEnd[ordering_end])


                # analyse begining positions of duplications/deletions
                for (i in 1:length(dupdel_variant_positions)) {
                    row_left = list(
                        chromosome,
                        start(vcf_sv[matched_dupdel])[i],
                        rownames(vcf_sv[matched_dupdel])[i],
                        "N",
                        unlist(rowRanges(vcf_sv)$ALT)[matched_dupdel][i],
                        unlist(rowRanges(vcf_sv)$QUAL)[matched_dupdel][i],
                        info(vcf_sv)$SVTYPE[matched_dupdel][i],
                        geno(vcf_sv)$GT[matched_dupdel][i],
                        info(vcf_sv)$END[matched_dupdel][i]
                    )

                    # genes overlapping with variant:
                    # - start before variant end
                    # - end after variant begin
                    startingBeforeVariantEnd = ordering_start[1:variantStart.geneStart[i]]
                    endingAfterVariantStart = ordering_end[(variantStart.geneEnd[i] +
                                                                1):length(ordering_end)]
                    overlappingGenes = intersect(startingBeforeVariantEnd,
                                                 endingAfterVariantStart)


                    # proceed only when position intersects with some gene
                    if (length(overlappingGenes) > 0) {
                        uniqueMask = !duplicated(anno[overlappingGenes,]$name2)

                        uniqueCommon = overlappingGenes[uniqueMask]

                        # for all genes containing the position
                        for (ig in 1:length(uniqueCommon)) {
                            # extract gene name
                            gene_id = uniqueCommon[ig]

                            gene = anno$name2[gene_id]

                            exon_starts = as.numeric(unlist(
                                strsplit(anno$exonStarts[gene_id], split = ",")
                            ))
                            exon_ends = as.numeric(unlist(strsplit(
                                anno$exonEnds[gene_id], split = ","
                            )))

                            # interleaved vector of low and high exon borders
                            gene_ranges = c(rbind(exon_starts, exon_ends))
                            variantBegin.interval = findInterval(dupdel_variant_positions[i],
                                                                 gene_ranges)
                            variantEnd.interval = findInterval(dupdel_variant_ends[i], gene_ranges)

                            in_exon = ((variantBegin.interval %% 2 == 1) |
                                           (variantEnd.interval %% 2 == 1) |
                                           (
                                               variantEnd.interval - variantBegin.interval
                                           ) >= 1
                            )
                            uniprot_row = uniprot[which(uniprot$Gene.names...primary.. == anno[uniqueCommon[ig], ]$name2),]

                            if (dim(uniprot_row)[1] == 0) {
                                uniprot_row = rep(NA, 11)
                            }

                            row = c(row_left, gene, in_exon, uniprot_row)
                        }
                    } else {
                        row = c(row_left, rep(NA, 13))
                    }

                    #instead of rbind which messes row names for empty frames
                    sv_local_table[nrow(sv_local_table) + 1, ] = row
                } # end position loop
            }

            ###############################################
            # Part 9 - Annotate point structural variants
            ###############################################

            if (length(point_variant_positions) > 0) {
                # find genes containing analysed positions
                intervals_start = findInterval(point_variant_positions, anno$txStart[ordering_start])

                intervals_end = findInterval(point_variant_positions, anno$txEnd[ordering_end])

                # analyse begining/ending positions of duplications/deletions
                for (i in 1:length(point_variant_positions)) {
                    row_left = list(
                        chromosome,
                        start(vcf_sv[matched_point])[i],
                        rownames(vcf_sv[matched_point])[i],
                        "N",
                        unlist(rowRanges(vcf_sv)$ALT)[matched_point][i],
                        unlist(rowRanges(vcf_sv)$QUAL)[matched_point][i],
                        info(vcf_sv)$SVTYPE[matched_point][i],
                        geno(vcf_sv)$GT[matched_point][i],
                        info(vcf_sv)$END[matched_point][i]
                    )

                    # consider variant begin positions
                    starting_before = ordering_start[1:intervals_start[i]]
                    ending_after = ordering_end[(intervals_end[i] + 1):length(ordering_end)]
                    common = intersect(starting_before, ending_after)

                    # proceed only when position intersects with some gene
                    if (length(common) > 0) {
                        uniqueMask = !duplicated(anno[common,]$name2)

                        uniqueCommon = common[uniqueMask]

                        # for all genes containing the position
                        for (ig in 1:length(uniqueCommon)) {
                            # extract gene name
                            gene_id = uniqueCommon[ig]

                            gene = anno$name2[gene_id]

                            exon_starts = as.numeric(unlist(
                                strsplit(anno$exonStarts[gene_id], split = ",")
                            ))
                            exon_ends = as.numeric(unlist(strsplit(
                                anno$exonEnds[gene_id], split = ","
                            )))

                            # interleaved vector of low and high exon borders
                            gene_ranges = c(rbind(exon_starts, exon_ends))
                            interval = findInterval(point_variant_positions[i],
                                                    gene_ranges)

                            exon = (interval %% 2 == 1)
                            uniprot_row = uniprot[which(uniprot$Gene.names...primary.. == anno[uniqueCommon[ig], ]$name2),]

                            if (dim(uniprot_row)[1] == 0) {
                                uniprot_row = rep(NA, 11)
                            }

                            row = c(row_left, gene, exon, uniprot_row)
                        }
                    } else {
                        row = c(row_left, rep(NA, 13))
                    }

                    #instead of rbind which messes row names for empty frames
                    sv_local_table[nrow(sv_local_table) + 1, ] = row

                } # end position loop
            }

            sv_local_table = sv_local_table[order(as.numeric(sv_local_table[, 2])), ]
            sv_local_table = sv_local_table[sv_local_table$QUAL > 0, ]

            variation_sv = rep(0, length(sv_local_table$start_position))
            variation_sv[sv_local_table$GT == "0/1"] = 0.5
            variation_sv[sv_local_table$GT == "1/1"] = 1

            sv_table = rbind(sv_table, sv_local_table)
        }

        ###############################################
        # Part 10 - Plots
        ###############################################

        # prepare png
        png(
            filename = paste(sampleName, "_chr", chromosome,
                             ".png", sep = ""),
            bg = "white",
            width = pngWidth,
            height = pngHeight,
            units = "px",
            pointsize = 24
        )

        # plot all variants
        plot(
            positions,
            variation,
            cex = 0.2,
            pch = 16,
            col = "blue",
            cex.lab = 1,
            xlab = paste0("Chromosome_", chromosome),
            ylab = "Number of alternative reads / depth",
            ylim = c(0, 1),
            xlim = c(0, max(positions))
        )

        if (chromosome != "Y") {
            MA <- movingAverage(variation, MA_Window)
            lines(positions, MA, col = "red", lwd = 2)
        }
        abline(h = 0.75,
               lwd = 4,
               col = "red")
        abline(v = centromeres$Start[chromosomeIndex],
               lwd = 4,
               col = "orange")
        abline(v = centromeres$Stop[chromosomeIndex],
               lwd = 4,
               col = "orange")

        # plot nonsynonymous/frameshift/nonesene coding variants, simple, rare in dbsnp and exac
        points(
            final_positions,
            final_variations,
            cex = 3,
            pch = 16,
            col = "green"
        )

        if (!missing(sv_sample)) {
            # plot structural variants that passed the filter
            points(
                sv_local_table$start_position,
                variation_sv,
                cex = 3,
                pch = 16,
                col = "orange"
            )
        }

        dev.off()

    } # end chromosome loop

    final_table[final_table == ""] = NA
    final_table = final_table[-1,] # remove dummy first row
    final_table$Gene.names...primary.. <-
        NULL # remove duplicated column with gene name

    # write to output
    write.table(
        final_table,
        paste0("RareVariants_", sampleName, ".txt"),
        sep = "\t",
        row.names = FALSE,
        quote = FALSE
    )

    # complex_table[complex_table == ""] = NA
    # complex_table = complex_table[-1,] # remove dummy first row
    #
    #
    # # write to output
    # write.table(
    #     complex_table,
    #     paste0("ComplexVariants_", sampleName, ".txt"),
    #     sep = "\t",
    #     row.names = FALSE,
    #     quote = FALSE
    # )


    if (!missing(sv_sample)) {
        sv_table[sv_table == ""] = NA
        sv_table = sv_table[-1,] # remove dummy first row
        sv_table$Gene.names...primary.. <-
            NULL # remove duplicated column with gene name

        write.table(
            sv_table,
            paste0("StructuralVariants_", sampleName, ".txt"),
            sep = "\t",
            row.names = FALSE,
            quote = FALSE
        )
    }
}

