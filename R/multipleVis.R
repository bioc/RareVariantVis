multipleVis <- function(filesList, filesOrder=1:length(filesList), sampleNames = filesList, chromosome, centromeres, inHouse=TRUE, inHouseSampNr=5){



 if(is.character(filesList[1])){
   sampleNames=sampleNames[filesOrder]
   filesList=filesList[filesOrder]
   filesData=sapply(filesList, read.table)
 }
 else
 {
   #example
   ln=length(unique(names(filesList)))
   d1=as.matrix(filesList[1:ln])
   d2=as.matrix(filesList[(ln+1):(ln*2)])
   d3=as.matrix(filesList[(ln*2+1):(ln*3)])
   d4=as.matrix(filesList[(ln*3+1):(ln*4)])
   d5=as.matrix(filesList[(ln*4+1):(ln*5)])
   filesData=cbind(d1,d2,d3,d4,d5)
   colnames(filesData) <- sampleNames
   filesList=colnames(filesData)
 }
 HeadNames=rownames(filesData)

 cat("Input file reading...\n")

 Scatter=c()
AD.all=NULL;
Start.position.all=NULL;
Gene.name.all=NULL;
Variant.all=NULL;
Reference.all=NULL;

 for (i in 1:length(HeadNames))
  {
       assign(paste(HeadNames[i], ".all", sep=""),filesData[HeadNames[i],])
  }


 # Prepare data
 for(j in 1:length(filesList))
   {

 ADs1 = as.character(AD.all[[j]])
 ADDs1 = strsplit(ADs1, ",")
 ADS1 <-t(sapply(ADDs1, '[', 1:max(sapply(ADDs1, length))))
 variation1 = as.numeric(ADS1[,2]) / (as.numeric(ADS1[,1]) + as.numeric(ADS1[,2]))

 C1 = data.frame(as.numeric(Start.position.all[[j]]), variation1)
 colnames(C1) = c("Start.position", sampleNames[j])

 tooltip = paste("Gene: ",as.character(Gene.name.all[[j]]),
                  "<br>Position: ", Start.position.all[[j]],
                  "<br>Alt/Depth: ", round(variation1,2),
                   sep="")

 C1$pop.html.tooltip = tooltip


 Scatter[[j]] <- gvisScatterChart(C1,
                                options=list(
                                title=paste(sampleNames[j],"_chr_",chromosome, sep = ""),
                                explorer="{actions: ['dragToZoom',
                                'rightClickToReset'],
                                maxZoomIn:0.05}",
                                #chartArea="{width:'85%',height:'80%'}",
                                tooltip="{isHtml:'True'}",
                                series="[{targetAxisIndex: 0}, {targetAxisIndex:1}]",
                                crosshair="{trigger:'both'}",
                                legend="none", lineWidth=0, pointSize=8,
                                vAxis="{title:'Number of alternative reads / depth',
                                viewWindow:{min:0, max:1}}",
                                hAxis=paste("{title:'Position on chromosome', viewWindow:{min:0, max:", centromeres$WinSize[chromosome] ,"}}", sep = ""),
                                width=1500, height=200))


  }

Scatter1=Scatter[[1]]

 for(k in 2:length(filesList))
  {
    Scatter1 = gvisMerge(Scatter1, Scatter[[k]], horizontal = FALSE,tableOptions = "border=\"0\"")
  }
today <- Sys.time()
if(inHouse==TRUE){
#table
  snp=unname(unlist(Start.position.all))
  var=unname(unlist(Variant.all))
  ref=unname(unlist(Reference.all))
  dfVar<-data.frame(snp,var)
  tabSnp=table(dfVar)
  varSum=rowSums(unname(tabSnp))
  thrSumVar=which(varSum>=inHouseSampNr)
  snpref=paste(ref[thrSumVar],"/",var[thrSumVar],sep="")
  inHouseFrame=data.frame(snp[thrSumVar],snpref,varSum[thrSumVar])
  names(inHouseFrame)=c("position","snp","repetition")
  write.table(unique(inHouseFrame[order(inHouseFrame$position),]),paste(format(today, format="%d%m%Y"),"_", "inHouseTable",".txt"), sep="\t",row.names=FALSE)
}


cat(Scatter1$html$chart, file=paste(format(today, format="%d%m%Y"),"_", "multipleVis_", chromosome, ".html", sep = ""))
cat("Analysis finished.\n")
cat(paste("Your output files are in folder:\n"), getwd(), "\n", sep = "")
}
