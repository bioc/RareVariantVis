chromosomeVis <- function(file, sampleName, chromosome, centromeres, pngWidth = 1600, pngHeight = 1200, plot = FALSE, frFilter = 0.01, dpFilter = 10, vcf = FALSE, posFilter = 0) {

    # input reading
    cat("Input file reading...\n")
    if (vcf == FALSE){
        # txt file
        if (is.character(file)){
            Table <- read.table(file, header = TRUE, sep="\t", quote="\"", fill = TRUE)
            FirstTable <- Table[order(Table$Start.position),]
        }

        # build in file
        else{

            FirstTable = file
        }
    }

    else{
        # vcf
        vcf <- readVcf(file, genome="hg19")

        firstADs = sapply( geno(vcf)[["AD"]], function(m) m[1] )
        secondADs = sapply( geno(vcf)[["AD"]], function(n) n[2] )
        combinedADs = paste(firstADs, secondADs, sep = ",")

        Table <- as.data.frame(cbind(start(rowRanges(vcf)),
                            as.character(info(vcf)[["Variant.type"]]),
                            unlist(info(vcf)[["SNP.Frequency"]]),
                            as.character(info(vcf)[["Gene.name"]]),
                            as.character(info(vcf)[["Gene.component"]]),
                            unlist(info(vcf)[["DP"]]),
                            combinedADs))

        colnames(Table) = c("Start.position", "Variant.type", "SNP.Frequency", "Gene.name", "Gene.component", "DP", "AD")

        FirstTable = data.frame(lapply(Table, as.character), stringsAsFactors=FALSE)
        FirstTable$Start.position = as.numeric(FirstTable$Start.position)
        FirstTable$SNP.Frequency = as.numeric(FirstTable$SNP.Frequency)
        FirstTable$DP = as.numeric(FirstTable$DP)
        FirstTable$Variant.type = as.character(FirstTable$Variant.type)

    }

    # testing if input file is OK
    columnTest = any(is.na(match(c("SNP.Frequency", "DP", "AD", "Start.position", "Gene.name", "Gene.component", "Variant.type"), colnames(FirstTable))))
    if(columnTest == FALSE){
        # if there is more than 2 variants in the same position, only one is kept
        FirstTable = FirstTable[!duplicated(FirstTable$Start.position),]
        cat("Input file reading completed.\n")
    }

    else{
        stop("Your input file/data frame must include Start.position, AD, DP, SNP.Frequency, Gene.name, Gene.component, Variant.type columns")
    }

    # Variables
    DBSNPfrequency = as.numeric(as.character(FirstTable$SNP.Frequency))
    DPs = as.numeric(as.character(FirstTable$DP))
    ADs = as.character(FirstTable$AD)
    ADDs = strsplit(ADs, ",")
    ADS <-t(sapply(ADDs, '[', 1:max(sapply(ADDs, length))))
    variation = as.numeric(ADS[,2]) / (as.numeric(ADS[,1]) + as.numeric(ADS[,2]))
    # PhyloP = FirstTable$phyloP

    # Filtering
    F1 = which(DBSNPfrequency < frFilter)
    F2 = which(FirstTable$Gene.component == "EXON_REGION" |
            FirstTable$Gene.component == "SA_SITE_CANONICAL" |
            FirstTable$Gene.component == "SD_SITE_CANONICAL" |
            FirstTable$Gene.component == "UTR")
    F3 = which(FirstTable$Variant.type == "Substitution - nonsynonymous" |
            FirstTable$Variant.type == "Substitution - nonsense" |
            FirstTable$Variant.type == "Complex" |
            FirstTable$Variant.type == "Deletion - frameshift" |
            FirstTable$Variant.type == "Insertion - frameshift")
    F4 = which(DPs > dpFilter)
    F5 = which(FirstTable != posFilter)
    Filters = intersect(intersect(intersect(intersect(F1,F2),F3),F4),F5)

    # moving average calculation
    MA <- movingAverage(variation, 2000)

    # plotting
    if(plot == FALSE){
    png(filename = paste(sampleName, "_chr", chromosome, ".png", sep = ""),
            bg = "white", width = pngWidth, height = pngHeight, units = "px", pointsize = 24)
    }

    plot(FirstTable$Start.position, variation, cex=0.2, pch=16, col="blue",
            cex.lab = 0.7, xlab = paste("Chromosome_", chromosome, sep=""),
            ylab = "Number of alternative reads / depth",
            main = paste("Sample_", sampleName, sep = ""), ylim = c(0,1),
            xlim = c(0, max(FirstTable$Start.position)))
    points(FirstTable$Start.position[Filters], variation[Filters], cex=2, col = "green", pch = 16)
    lines(FirstTable$Start.position, MA, col="red", lwd=2)
    abline(h=0.75, lwd=4, col="red")
    abline(v=centromeres$Start[chromosome], lwd=4, col="orange")
    abline(v=centromeres$Stop[chromosome], lwd=4, col="orange")

    # coverage plot



    if (plot == FALSE){
        dev.off()
    }

    # rare variant reporting

    write.table(FirstTable[Filters,], paste("RareVariants_",
            sampleName, "_chromosome_", chromosome, ".txt", sep="") ,
            sep='\t', row.names = 1:dim(FirstTable[Filters,])[1])

    cat("Analysis finished.\n")
    cat(paste("Your output files are in folder:\n"), getwd(), "\n", sep = "")
}
