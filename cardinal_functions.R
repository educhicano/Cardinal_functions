##use common peaklist to do reduceDimension over multiple different datasets in one folder in parallel
##normalization by TIC is default

ReduceDimensionParallel <- function(peaklist, mywd_=getwd(),no_clus=4){
  
  imzmllist <- gsub(".imzML","",list.files(path=mywd,pattern="*.imzML"))
  imzmllist <- as.list(c(imzmllist))
  ppeaklist <- peaklist
  
  require(snow)
  
  cl <- makeCluster(no_clus)
  clusterEvalQ(cl, library(Cardinal)) 
  clusterExport(cl, "imzmllist") 
  clusterExport(cl, "ppeaklist")
  
  peakdata <- parLapply(cl, imzmllist, function(x){
    rawd <-readImzML(x,folder=mywd)
    rawd <-normalize(rawd, method=c("tic"))
    peakd <- reduceDimension(rawd, ref = ppeaklist, type = "height")
    rm(rawd)
    gc()
    return(peakd)
  })
  
  stopCluster(cl)
  return(peakdata)}

##subset class data from a spatialShrunkenCentroids analysis or other clustering analysis on original data
##useful for removing MALDI noise clusters from dataset, or other well-defined background
##cardinaldata: original dataset where clustering was performed
##classvector example: sscgdata[[1]]$classes
##indicestokeep: vector of class number to keep, for instance, to keep classes 1-5, c(1:5)

reducebysegment<- function(cardinaldata, classvector, indicestokeep){
  cardinaldata$class <- classvector
  cardinaldata <- cardinaldata[, cardinaldata$class %in% indicestokeep]
  return(cardinaldata)
}

##multi-plot segment loading data from spatialShrunkenCentroids analysis
##sscgdata: your spatial shrunken analysis data, nsegments=number of segments in your analysis
##pmode: "tstatistics" "centers"
##playout: vector of layout, for instance
##modelno is the model used if multiple sscgs were performed
##plot_segments(sscgdata,modelno=2,playout=c(3,3),pmode="tstatistics")
##plots tstats of second model in a 3x3 grid
##grid size must be larger than number of segments

plot_segments <- function(sscgdata,modelno=1,playout,pmode){
  nsegments <- nlevels(sscgdata@resultData[[modelno]]$classes)
  stopifnot(playout[1]*playout[2]  >= nsegments)
  plot(sscgdata, model=modelno,mode=pmode, layout=playout, superpose=FALSE, column=1,col=rainbow(nsegments)[1])
  for(i in 1:(nsegments-1)){
    plot(sscgdata,model=modelno,mode=pmode,superpose=FALSE,column=i+1,col=rainbow(nsegments)[1+i])
  }
}

##sscgdata: your spatial shrunken centroids analysis data, no_markers: the number of top markers you wish to display, must be < than number of picked peaks
##modelno is the model used if multiple sscgs were performed
##no_markers: how many top markers, top ten: no_markers=10

sscg_topSegMarkers <- function(sscgdata,modelno=1,no_markers){
  nsegments <- nlevels(sscgdata@resultData[[modelno]]$classes)
  df <- data.frame(sscgdata@resultData[[modelno]]$tstatistics)
  colnames(df) <- paste("segment_",c(1:nsegments),sep="")
  dfs <- list()
  for(i in 1:nsegments){
    dfs[[i]] <- df[order(df[,i],decreasing = TRUE),]
    dfs[[i]] <- dfs[[i]][1:no_markers,]
  }
  names(dfs) <- paste("segment_",c(1:nsegments),sep="")
  for(i in 1:nsegments){
    print(dfs[[i]][i])
  }
  return(dfs)
}

##import spotlist from Bruker FlexImaging XML, changing Cardinal SampleNames to those in the XML
##can take multiple regions, but not overlapping regions
##Prints "TRUE" if the reduction sampling worked

sampleROIs<-function(cardinaldata, XMLdata="myxml.xml"){
  
  require(XML)
  
  XMLspots<-xmlToList(XMLdata)
  
  spots <-lapply(XMLspots[which(names(XMLspots) == "Class")],function(x){
    spots <- list()
    for(i in 1:(length(x == "Element") -1)){
      spots[[i]] <- x[[i]][["Spot"]]
    }
    spots <-unlist(spots)
    spots <-substr(spots,6,100)
    return(spots)})
  names(spots) <- NULL
  sampnames <- sapply(XMLspots[which(names(XMLspots) == "Class")],function(x)x$.attrs[["Name"]])
  
  
  cardinaldatas <- lapply(spots,function(x){
    cardinaldataf <- cardinaldata[,paste("X",coord(cardinaldata)[,1],"Y",coord(cardinaldata)[,2],sep="") %in% x]
    return(cardinaldataf)
  })
  
  for(i in 1:length(cardinaldatas)){
    levels(cardinaldatas[[i]]$sample) <- sampnames[i]
    protocolData(cardinaldatas[[i]]) <- AnnotatedDataFrame(data=data.frame(row.names=sampleNames(cardinaldatas[[i]])))
    cardinaldatas[[i]] <- regeneratePositions(cardinaldatas[[i]])
    print(validObject(cardinaldatas[[i]]))
  }
  
  cardinaldata<-do.call(combine,cardinaldatas)
  return(cardinaldata)}
