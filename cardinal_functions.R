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
##reduces cardinaldata to only the ROIs, for instance, if there are multiple ROIs in the XML
##function will only take data in ROIs, all non selected data is excluded
##can take multiple regions, but not overlapping regions
##Prints "TRUE" if the reduction sampling worked

reduce_by_XML<-function(cardinaldata, XMLdata="myxml.xml"){
  
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

##This will add a factor to the pixeldata of a cardinal ImageSet with the regions from a Bruker FlexImaging XML spotlist
##There is no reduction performed
##region names will be taken from name given in flexImaging to the XML file, defaults for first ROI is "ROI 01"
##if x and y coordinates are > 3 digits, change "width" argument in formatC() function

XML_ROIs<-function(cardinaldata, spotlist="myxml.xml"){
  
  require(XML)
  
  XMLspots<-xmlToList(spotlist)
  
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

  cardinaldata$regions <- NA
  
  for(i in seq_along(spots)){
  cardinaldata$regions[paste0("X",formatC(coord(cardinaldata)[,1], width = 3, format = "d", flag = "0"),"Y",formatC(coord(cardinaldata)[,2], width = 3, format = "d", flag = "0")) %in% spots[[i]]] <- sampnames[i]
  }
  
  cardinaldata$regions <- as.factor(cardinaldata$regions)
  
  return(cardinaldata)}


##get a table of mean, sd, cv for all mass bins in a cardinal dataset by a subset, meant for picked peak data
##grouping is how the data is to be subsetted and is the name of the vector in the cardinal pixelData
##to go simply by same name grouping="sample"

ROI_analysis <- function(cardinaldata,grouping="regions"){
  
  cardinaldata_df <- data.frame(pData(cardinaldata)[grouping],t(iData(cardinaldata)))
  colnames(cardinaldata_df) <- c(paste(grouping),fData(cardinaldata)$mz)
  cardinaldata_df<-na.omit(cardinaldata_df)[1:length(cardinaldata_df)]
  
  melted <- melt(cardinaldata_df,variable.name = "mass")
  analysisoutput<-ddply(melted, c(paste0(grouping), "mass"), summarise,
                        mean = mean(value), sd = sd(value),
                        cv = sd(value)/mean(value)
  )
  return(analysisoutput) 
}
