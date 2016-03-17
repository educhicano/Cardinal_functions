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


##returns a factor in the pData of a cardinal imageSet titled "regions"
##with the regions of a Bruker FlexImaging .xml

XML_ROIs<-function(cardinaldata, spotlist="myxml.xml",width=3){
  
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
