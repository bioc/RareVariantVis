callHomozygous <- function(
    sample,
    chromosomes,
    caller = "speedseq",
    MA_Window = 2000,
    HMZ_length = 100000,
    min_n_HMZ = 20) {

    sampleName = unlist(strsplit(basename(sample), '[.]'))[1]

    # Annotation files
    centromeres_file    = system.file("extdata", "CentromeresHg19.txt", package =
                                          "RareVariantVis")

    # load annotations
    centromeres = read.delim(centromeres_file,
                             header = TRUE,
                             stringsAsFactors = FALSE)

    dir.create(paste0(getwd(),"/HMZ_Calling_Output_", sampleName))

    # create empty table for output
    output_table = data.frame(
        sample = NA,
        chromosome = "",
        start_positions = NA,
        end_positions = NA,
        length = NA,
        var_in_reg = NA,
        perc_HMZ = NA,
        stringsAsFactors = FALSE
    )


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
        # Part 2 - Compute and plot
        ###############################################

        setwd(paste0(getwd(),"/HMZ_Calling_Output_", sampleName))

        if (chromosome != "Y") {
            MA <- movingAverage(variation, MA_Window)

            ### HMZ calling

            g075 = MA > 0.75
            g075_xor = NULL

            for (j in 1:length(g075) - 1) {
                g075_xor[j] = xor(g075[j], g075[j + 1])

            }

            starts_hmz = g075[2:length(g075)] & g075_xor
            starts_hmz_up = rep(FALSE, length(starts_hmz) + 1)
            starts_hmz_up[1] = MA[1] > 0.75
            starts_hmz_up[2:length(starts_hmz_up)] = starts_hmz

            stops_hmz = !g075[2:length(g075)] & g075_xor
            stops_hmz_up = rep(FALSE, length(starts_hmz) + 1)
            stops_hmz_up[1:length(stops_hmz)] = stops_hmz
            stops_hmz_up[length(stops_hmz_up)] = MA[(length(MA))] > 0.75

            sel_HMZ = (positions[stops_hmz_up] - positions[starts_hmz_up]) > HMZ_length &
                which(stops_hmz_up) - which(starts_hmz_up) > min_n_HMZ

            ### plotting

            png(
                filename = paste(sampleName, "_chr", chromosome,
                                 ".png", sep = ""),
                bg = "white",
                width = 1600,
                height = 1200,
                units = "px",
                pointsize = 24
            )

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

            lines(positions, MA, col = "red", lwd = 2)

            abline(h = 0.75, lwd = 4, col = "red")
            abline(v = centromeres$Start[chromosomeIndex],
                   lwd = 4,
                   col = "orange")
            abline(v = centromeres$Stop[chromosomeIndex],
                   lwd = 4,
                   col = "orange")

            abline(v = positions[starts_hmz_up][sel_HMZ],
                   lwd = 1,
                   col = "green")
            abline(v = positions[stops_hmz_up][sel_HMZ],
                   lwd = 1,
                   col = "red")

            dev.off()

        } else {
            png(
                filename = paste(sampleName, "_chr", chromosome,
                                 ".png", sep = ""),
                bg = "white",
                width = 1600,
                height = 1200,
                units = "px",
                pointsize = 24
            )

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

            abline(h = 0.75, lwd = 4, col = "red")
            abline(v = centromeres$Start[chromosomeIndex],
                   lwd = 4,
                   col = "orange")
            abline(v = centromeres$Stop[chromosomeIndex],
                   lwd = 4,
                   col = "orange")

            dev.off()
        }

        ###############################################
        # Part 3 - Report
        ###############################################

        if (length(positions[starts_hmz_up][sel_HMZ]) > 0) {
            perc_HMZ = NULL

            if (chromosome != "Y") {
                for (p in 1:length(positions[starts_hmz_up][sel_HMZ])) {
                    variation_curr_hmz = variation[which(starts_hmz_up)[sel_HMZ][p]:which(stops_hmz_up)[sel_HMZ][p]]
                    perc_HMZ[p] = sum(variation_curr_hmz > 0.75) / length(variation_curr_hmz)
                }
            } else {
                perc_HMZ = NA
            }

            current_table = cbind(
                rep(sampleName, length(positions[starts_hmz_up][sel_HMZ])),
                rep(chromosome, length(positions[starts_hmz_up][sel_HMZ])),
                positions[starts_hmz_up][sel_HMZ],
                positions[stops_hmz_up][sel_HMZ],
                positions[stops_hmz_up][sel_HMZ] - positions[starts_hmz_up][sel_HMZ],
                (which(stops_hmz_up) - which(starts_hmz_up))[sel_HMZ],
                round(perc_HMZ, 4) * 100
            )
            colnames(current_table) = colnames(output_table)


            output_table = rbind(output_table, current_table)
        }

    }

    final_table = output_table[!is.na(output_table$start_positions), ]

    write.table(
        final_table,
        paste0("Regions_HMZ_", sampleName, ".txt"),
        sep = "\t",
        row.names = FALSE,
        quote = FALSE
    )

}
