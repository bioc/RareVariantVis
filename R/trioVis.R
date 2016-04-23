trioVis <- function(fileMother, fileIndex, fileFather, sampleMother, sampleIndex, sampleFather, chromosome, centromeres) {

    cat("Input file reading...\n")

    # Mother
    if (is.character(fileMother)){
        Table1 <- read.table(fileMother, header = TRUE, sep="\t", quote="\"", fill = TRUE)
        FirstTable1 <- Table1[order(Table1$Start.position),]
    }

    else{
        FirstTable1 = fileMother
    }

    # Index
    if (is.character(fileIndex)){
        Table2 <- read.table(fileIndex, header = TRUE, sep="\t", quote="\"", fill = TRUE)
        FirstTable2 <- Table2[order(Table2$Start.position),]
    }

    else{
        FirstTable2 = fileIndex
    }

    # Father
    if (is.character(fileFather)){
        Table3 <- read.table(fileFather, header = TRUE, sep="\t", quote="\"", fill = TRUE)
        FirstTable3 <- Table3[order(Table3$Start.position),]
    }

    else{
        FirstTable3 = fileFather
    }

    cat("Input file reading completed.\n")

    # Prepare data
    ADs1 = as.character(FirstTable1$AD)
    ADDs1 = strsplit(ADs1, ",")
    ADS1 <-t(sapply(ADDs1, '[', 1:max(sapply(ADDs1, length))))
    ADs2 = as.character(FirstTable2$AD)
    ADDs2 = strsplit(ADs2, ",")
    ADS2 <-t(sapply(ADDs2, '[', 1:max(sapply(ADDs2, length))))
    ADs3 = as.character(FirstTable3$AD)
    ADDs3 = strsplit(ADs3, ",")
    ADS3 <-t(sapply(ADDs3, '[', 1:max(sapply(ADDs3, length))))

    variation1 = as.numeric(ADS1[,2]) / (as.numeric(ADS1[,1]) + as.numeric(ADS1[,2]))
    variation2 = as.numeric(ADS2[,2]) / (as.numeric(ADS2[,1]) + as.numeric(ADS2[,2]))
    variation3 = as.numeric(ADS3[,2]) / (as.numeric(ADS3[,1]) + as.numeric(ADS3[,2]))

    polyMother = FirstTable1$Start.position
    polyIndex  = FirstTable2$Start.position
    polyFather = FirstTable3$Start.position

    motherPolyCount = 0
    indexPolyCount  = 0
    fatherPolyCount = 0

    for (i in 1:length(polyMother)){
        motherPolyCount[i] = sum(polyMother > polyMother[i] - 100000 & polyMother < polyMother[i] + 100000)
    }

    for (i in 1:length(polyIndex)){
        indexPolyCount[i] = sum(polyIndex > polyIndex[i] - 100000 & polyIndex < polyIndex[i] + 100000)
    }

    for (i in 1:length(polyFather)){
        fatherPolyCount[i] = sum(polyFather > polyFather[i] - 100000 & polyFather < polyFather[i] + 100000)
    }

    motherPolyVariant = motherPolyCount > 3
    indexPolyVariant  = indexPolyCount  > 3
    fatherPolyVariant = fatherPolyCount > 3

    C1 = cbind(FirstTable1$Start.position, variation1, motherPolyVariant)
    colnames(C1) = c("Start.position", "Mother", "MotherPoly")

    C2 = cbind(FirstTable2$Start.position, variation2, indexPolyVariant)
    colnames(C2) = c("Start.position", "Index", "IndexPoly")

    C3 = cbind(FirstTable3$Start.position, variation3, fatherPolyVariant)
    colnames(C3) = c("Start.position", "Father", "FatherPoly")

    M1 = merge(C1, C2, by = "Start.position", all = TRUE)
    M2 = merge(M1, C3, by = "Start.position", all = TRUE)

    # Inheritance

    Data1 = M2[,c(1,2)]
    Data2 = M2[,c(1,4)]
    Data3 = M2[,c(1,6)]

    Vect1 = as.character(!is.na(M2$Mother))
    Vect2 = as.character(!is.na(M2$Index))
    Vect3 = as.character(!is.na(M2$Father))

    both = !is.na(M2$Index) & !is.na(M2$Mother) & !is.na(M2$Father)
    moth = !is.na(M2$Index) & !is.na(M2$Mother) & is.na(M2$Father)
    fath = !is.na(M2$Index) & is.na(M2$Mother) & !is.na(M2$Father)
    dnov = !is.na(M2$Index) & is.na(M2$Mother) & is.na(M2$Father)

    inhe = both
    inhe[inhe == "TRUE"] = "Inherited from both parents"
    inhe[moth == "TRUE"] = "Inherited from mother"
    inhe[fath == "TRUE"] = "Inherited from father"
    inhe[dnov == "TRUE"] = "Possible De Novo"

    gene1 = Data1$Mother
    gene1[!is.na(gene1)] = as.character(FirstTable1$Gene.name)
    gene2 = Data2$Index
    gene2[!is.na(gene2)] = as.character(FirstTable2$Gene.name)
    gene3 = Data3$Father
    gene3[!is.na(gene3)] = as.character(FirstTable3$Gene.name)


    tooltip1 = paste("Gene: ", gene1,
                    "<br>Position: ", Data1$Start.position,
                    "<br>Alt/Depth: ", round(Data1$Mother, 2),
                    "<br>Problematic region: ", as.logical(M2$MotherPoly),
                    sep="")

    tooltip2 = paste("Gene: ", gene2,
                    "<br>Position: ", Data2$Start.position,
                    "<br>Alt/Depth: ", round(Data2$Index, 2),
                    "<br>Problematic region: ", as.logical(M2$IndexPoly),
                    "<br>", inhe,
                    sep="")

    tooltip3 = paste("Gene: ", gene3,
                    "<br>Position: ", Data3$Start.position,
                    "<br>Alt/Depth: ", round(Data3$Father, 2),
                    "<br>Problematic region: ", as.logical(M2$FatherPoly),
                    sep="")


    Data1$Mother.scope = dnov
    Data3$Father.scope = dnov

    Data1$pop.html.tooltip = tooltip1
    Data2$pop.html.tooltip = tooltip2
    Data3$pop.html.tooltip = tooltip3

    Data2$Index.certainty = !dnov
    Data2$Index.emphasis = dnov
    Data2$Index.scope = dnov


    # Plotting

    Scatter2 <- gvisScatterChart(Data1,
                    options=list(
                    title=paste("SampleMother_", sampleMother, "_chr_", chromosome, sep = ""),
                    explorer="{actions: ['dragToZoom',
                    'rightClickToReset'],
                    maxZoomIn:0.05}",
                    #chartArea="{width:'85%',height:'80%'}",
                    tooltip="{isHtml:'True'}",
                    series="[{color: 'red', targetAxisIndex: 0}, {targetAxisIndex:1}]",
                    crosshair="{trigger:'both'}",
                    legend="none", lineWidth=0, pointSize=8,
                    vAxis="{title:'Number of alternative reads / depth',
                    viewWindow:{min:0, max:1}}",
                    hAxis=paste("{title:'Position on chromosome', viewWindow:{min:0, max:", centromeres$WinSize[chromosome] ,"}}", sep = ""),
                    width=1500, height=200))


    Scatter3 <- gvisScatterChart(Data2,
                    options=list(
                    title=paste("SampleIndex_", sampleIndex, "_chr_", chromosome, sep = ""),
                    explorer="{actions: ['dragToZoom',
                    'rightClickToReset'],
                    maxZoomIn:0.05}",
                    #chartArea="{width:'85%',height:'80%'}",
                    tooltip="{isHtml:'True'}",
                    series="[{color: 'red', targetAxisIndex: 0},{targetAxisIndex:1}]",
                    crosshair="{trigger:'both'}",
                    legend="none", lineWidth=0, pointSize=8,
                    vAxis="{title:'Number of alternative reads / depth',
                    viewWindow:{min:0, max:1}}",
                    hAxis=paste("{title:'Position on chromosome', viewWindow:{min:0, max:", centromeres$WinSize[chromosome] ,"}}", sep = ""),
                    width=1500, height=200))

    Scatter4 <- gvisScatterChart(Data3,
                    options=list(
                    title=paste("SampleFather_", sampleFather, "_chr_", chromosome, sep = ""),
                    explorer="{actions: ['dragToZoom',
                    'rightClickToReset'],
                    maxZoomIn:0.05}",
                    #chartArea="{width:'85%',height:'80%'}",
                    tooltip="{isHtml:'True'}",
                    series="[{color: 'red', targetAxisIndex: 0},
                    {targetAxisIndex:1}]",
                    crosshair="{trigger:'both'}",
                    legend="none", lineWidth=0, pointSize=8,
                    vAxis="{title:'Number of alternative reads / depth',
                    viewWindow:{min:0, max:1}}",
                    hAxis=paste("{title:'Position on chromosome', viewWindow:{min:0, max:", centromeres$WinSize[chromosome] ,"}}", sep = ""),
                    width=1500, height=200))


    ScatterA = gvisMerge(Scatter2, Scatter3, horizontal = FALSE,
                         tableOptions = "border=\"0\"")

    ScatterB = gvisMerge(ScatterA, Scatter4, horizontal = FALSE,
                         tableOptions = "border=\"0\"")

    cat(ScatterB$html$chart, file=paste(sampleMother, "_", sampleIndex, "_", sampleFather, "_", "TrioChr", chromosome, ".html", sep = ""))

    geneCount = table(FirstTable2$Gene.name)
    geneCountMatches = match(as.character(FirstTable2$Gene.name), rownames(geneCount))

    output = cbind(FirstTable2, inhe[which(inhe != "FALSE")], indexPolyVariant, geneCount[geneCountMatches])
    colnames(output)[c(dim(output)[2] - 2, dim(output)[2] - 1, dim(output)[2])] = c("Inheritance", "In polymorphic or problematic region", "Gene Count")

    write.table(output, paste("RareVariantsTrio_",
                              sampleIndex, "_chromosome_", chromosome, ".txt", sep=""), sep = "\t")

    cat("Analysis finished.\n")
    cat(paste("Your output files are in folder:\n"), getwd(), "\n", sep = "")

}
