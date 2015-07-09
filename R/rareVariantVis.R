rareVariantVis <- function(file, sample, chromosome, centromeres) {

    cat("Input file reading...\n")

    if (is.character(file)){
        Table <- read.table(file, header = TRUE, sep="\t", quote="\"", fill = TRUE)
        FirstTable <- Table[order(Table$Start.position),]
    }

    else{
        FirstTable = file
    }

    cat("Input file reading completed.\n")

    ADs = as.character(FirstTable$AD)
    ADDs = strsplit(ADs, ",")
    ADS <-t(sapply(ADDs, '[', 1:max(sapply(ADDs, length))))
    variation = as.numeric(ADS[,2]) / (as.numeric(ADS[,1]) + as.numeric(ADS[,2]))
    variationRound = round(variation, 2)

    Data = as.data.frame(cbind(FirstTable$Start.position, variation))
    Data$Gene.name = paste("Gene: ", as.character(FirstTable$Gene.name),
                           "<br>Position: ", FirstTable$Start.position,
                           "<br>Alt/Depth: ", variationRound, sep="")

    colnames(Data) = c("Start.position", paste("Sample_", sample, sep = ""),
        paste("Sample_", sample, ".html.tooltip", sep = ""))

    Scatter <- gvisScatterChart(Data,
                options=list(
                title=paste("Sample_", sample, "_chr_", chromosome, sep = ""),
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
                hAxis=paste("{title:'Position on chromosome', viewWindow:{min:0, max:", centromeres$WinSize[chromosome] ,"}}", sep = ""),
                #hAxis="{title:'Position on chromosome'}",
                width=1500, height=200))

    cat(Scatter$html$chart, file=paste(sample, "_chr", chromosome, ".html", sep = ""))

    cat("Analysis finished.\n")
    cat(paste("Your output files are in folder:\n"), getwd(), "\n", sep = "")

}
