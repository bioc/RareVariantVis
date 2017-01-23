rareVariantVis <- function(input, outputFile, sample, chromosomes = c(as.character(1:22), "X", "Y"), append = FALSE) {

    centromeres_file    = system.file("extdata", "CentromeresHg19.txt", package = "DataRareVariantVis")
    centromeres = read.delim(centromeres_file,
                             header = TRUE,
                             stringsAsFactors = FALSE)

    cat("Reading input:", input, "...")

    if (is.character(input)){
        entireTable <- read.table(input, header = TRUE, sep="\t", quote="\"", fill = TRUE, stringsAsFactors = FALSE)
    } else {
        entireTable = input
    }

    scatter = list()
    scatterId = 1

    cat("done!\n")
    # iterate over all chromosomes
    for (cid in 1:length(chromosomes)) {
        # get subtable corresponding to current chromosome
        chromosome = chromosomes[cid]

        if (chromosome == "X") {
            chromosomeIndex = 23
        } else if (chromosome == "Y") {
            chromosomeIndex = 24
        } else {
            chromosomeIndex = as.integer(chromosome)
        }

        table = entireTable[entireTable$chromosome == chromosome,]

        if (nrow(table) > 0) {
            variationRound = round(table$final_variations, 2)
            data = as.data.frame(cbind(table$final_positions, table$final_variations))
            data$Gene.name = paste("Gene: ", as.character(table$gene_name),
                                   "<br>Position: ", table$final_positions,
                                   "<br>Alt/Depth: ", variationRound, sep="")

            colnames(data) = c("Start.position", paste("Sample: ", sample, sep = ""),
                               paste("Sample: ", sample, ".html.tooltip", sep = ""))

            scatter[[scatterId]] = gvisScatterChart(data,
                        options=list(
                        title=paste("Sample: ", sample, ", chromosome: ", chromosome, sep = ""),
                        explorer="{actions: ['dragToZoom',
                        'rightClickToReset'],
                        maxZoomIn:0.05}",
                        #chartArea="{width:'85%',height:'80%'}",
                        tooltip="{isHtml:'True'}",
                        series="[{color: 'green', targetAxisIndex: 0},
                        {targetAxisIndex:1}]",
                        crosshair="{trigger:'both'}",
                        legend="none", lineWidth=0, pointSize=8,
                        vAxis="{title:'Number of alternative reads / depth',
                        viewWindow:{min:0, max:1}}",
                        hAxis=paste("{title:'Position on chromosome', viewWindow:{min:0, max:", centromeres$WinSize[chromosomeIndex] ,"}}", sep = ""),
                        #hAxis="{title:'Position on chromosome'}",
                        width=1500, height=200))

            localAppend = append
            if (scatterId > 1) {
                localAppend = TRUE
            }
            cat(scatter[[scatterId]]$html$chart , file=paste(outputFile, sep = ""), append=localAppend)

            scatterId = scatterId + 1;

        }
    }

    cat("Visualisation saved to:", outputFile, '\n')

}
