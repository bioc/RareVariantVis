multipleVis <- function(inputFiles, outputFile, sampleNames, chromosome) {

    for (i in 1:length(inputFiles)) {
        rareVariantVis(inputFiles[i], outputFile, sampleNames[i], c(chromosome), (i > 1))
    }
}
