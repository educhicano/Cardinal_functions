##outputs the .txt and .mir file for resultSets of Cardinal data
##if the pixelData contains multiple levels in $sample, each sample will have its own output
##useful for visualization, overlaying in cardinal results in FlexImaging

##model is the analysis performed, if multiple are performed, they can be conveniently accessed by 
##row.names(modelData(resultsdata))[model]
##and the number used for the model desired
##if you want to output multiple models/analyses, run the function in a for loop

##mode is the type of output:

###### SSCG:
## probabilities : probabilities in sscg analysis, also works for predicted
## classes : hard classes in sscg analysis, also works for predicted
## scores : scores in sscg analysis, also works for predicted

###### PLS:
## pls_scores : pls scores 
## pls_fitted : fitting of pls scores, also works for predicted

###### OPLS:
## opls_scores : opls scores 
## opls_fitted : fitting of opls scores, also works for predicted

###### PCA:
## pca : pca scores


flexOutputter <- function(resultsdata, model=1,mode="probabilities"){
  
  require(XML)
  
  input_model <- row.names(modelData(resultsdata))[model]
    
  
    if(mode == "probabilities"){
    
    rdata <-  round(resultData(resultsdata)[input_model][[input_model]]$probabilities,3)
    outputtype <- "probs"
    }

    if(mode == "scores"){
  
      rdata <-  round(resultData(resultsdata)[input_model][[input_model]]$scores,3)
      outputtype <- "scores"
    }


    if(mode == "classes"){
      
    rdata <- data.frame(resultData(resultsdata)[input_model][[input_model]]$classes)
    
    for(i in 1:resultData(resultsdata)[input_model][[input_model]]$k){
    rdata[i] <- as.integer(resultData(resultsdata)[input_model][[input_model]]$classes)
    rdata[i][rdata[i] != i] <- 0
  }
    
  outputtype <- "classes"
  }
  
  if(mode == "pca"){
  
  rdata <-  round(resultData(resultsdata)[input_model][[input_model]]$scores,3)
  
  outputtype <- "PCA_comp"
  
  }
  
  if(mode == "pls_scores"){
    rdata <-resultData(resultsdata)[input_model][[input_model]]$scores
    outputtype <- "PLS_scores"
  }
  
  if(mode == "pls_fitted"){
    rdata <-resultData(resultsdata)[input_model][[input_model]]$fitted
    outputtype <- "PLS_fitted"
  }
  
  if(mode == "opls_oscores"){
    rdata <-resultData(resultsdata)[input_model][[input_model]]$Oscores
    outputtype <- "OPLS_scores"
  }
  
  if(mode == "opls_fitted"){
    rdata <-resultData(resultsdata)[input_model][[input_model]]$fitted
    outputtype <- "OPLS_fitted"
  }
  
  
  
  if(mode == "probabilities" | mode == "scores" | mode == "classes"){
    if ((length(resultData(resultsdata)[input_model][[input_model]]$classes)  < length(resultData(resultsdata)[input_model][[input_model]]$y)) == TRUE){
      
      modelname<-paste0("sscg_predicted_")
    
      columnnames <- paste0(outputtype,"_k",1:resultData(resultsdata)[input_model][[input_model]]$k)
    
    }else{
    
    modelname<-paste0("sscg_r",resultData(resultsdata)[input_model][[input_model]]$r,
                      "k",resultData(resultsdata)[input_model][[input_model]]$k,
                      "s",resultData(resultsdata)[input_model][[input_model]]$s, "_")
    
    columnnames <- paste0(outputtype,"_k",1:resultData(resultsdata)[input_model][[input_model]]$k)
  }}
  
  if(mode == "pca"){
    
    modelname <- paste0("PCA_ncomp_",resultData(resultsdata)[input_model][[input_model]]$ncomp,"_")
    
    columnnames <- paste0(outputtype,"_",1:resultData(resultsdata)[input_model][[input_model]]$ncomp)
    
  }
  
  if(mode == "pls_scores"){
    
    modelname <- paste0("PLS_ncomp_",resultData(resultsdata)[input_model][[input_model]]$ncomp,"_")
    
    columnnames <- paste0(outputtype,"_",1:resultData(resultsdata)[input_model][[input_model]]$ncomp)
    
  }
  
  if(mode == "pls_fitted"){
    
    modelname <- paste0("PLS_y_fit_")
    
    columnnames <- paste0(outputtype,"_",1:nlevels(resultData(resultsdata)[input_model][[input_model]]$y))
    
  }
  
  if(mode == "opls_oscores"){
    
    modelname <- paste0("OPLS_ncomp_",resultData(resultsdata)[input_model][[input_model]]$ncomp,"_")
    
    columnnames <- paste0(outputtype,"_",1:resultData(resultsdata)[input_model][[input_model]]$ncomp)
    
  }
  
  if(mode == "opls_fitted"){
    
    modelname <- paste0("OPLS_y_fit_")
    
    columnnames <- paste0(outputtype,"_",1:nlevels(resultData(resultsdata)[input_model][[input_model]]$y))
    
  }
  
  if(mode == "pca"){
    
    modelname <- paste0("PCA_ncomp_",resultData(resultsdata)[input_model][[input_model]]$ncomp,"_")
    
    columnnames <- paste0(outputtype,"_",1:resultData(resultsdata)[input_model][[input_model]]$ncomp)
    
  }
  

  
  df_output <- data.frame(pData(resultsdata)$sample, 
                          paste0("0_R00X",pData(resultsdata)$x,"Y",pData(resultsdata)$y),
                          rdata)
  
  colnames(df_output) <- c("Sample","Pixel",
                           columnnames)
  
  splitdf <- split(df_output, df_output$Sample)
  
  lapply(names(splitdf), 
         function (x) write.table(splitdf[[x]], 
                                  file=paste0(modelname,
                                              x,"_",mode,".txt"), 
                                  quote=FALSE, col.names=TRUE, row.names=FALSE, sep=" "))
    
  XMLs <- list()
  
  for(j in seq_along(sampleNames(resultsdata))){
    XMLs[[j]] = xmlHashTree()
    ImagingResults = addNode(xmlNode("ImagingResults"), character(), XMLs[[j]])
    for (i in 1:length(columnnames)) {
      addNode(xmlNode("Result",attrs= c(Type= "TextImport", 
                                        Name=columnnames[i], 
                                        Color=NA, 
                                        Show="0", 
                                        MinIntensity="0", 
                                        IntensityThreshold="100", 
                                        AbsIntens="0", 
                                        LogScale="0", 
                                        InputFile=paste0(modelname,
                                                         sampleNames(resultsdata)[j],"_",mode,".txt"), 
                                        SpotCol="1", 
                                        DataCol=i+1)),
              ImagingResults, XMLs[[j]])
    }
  }
  
  
  for(i in seq_along(sampleNames(resultsdata))) {
    sink(paste0(modelname,
                sampleNames(resultsdata)[i],"_",mode,".mir")) 
    print(XMLs[[i]]) 
    sink()
  }
}
