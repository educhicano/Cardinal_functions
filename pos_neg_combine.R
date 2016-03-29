require(Cardinal)
require(plyr)

##script for combining positive and negative ionization mode IMS datasets acquired using x-y grid offset
##datasets must have the exact same x and y coordinates
##sample name will be the same as in the positive dataset
##the mass of range of negative data  1000 m/z higher(i.e. negative m/z 700 is 1700)
##if positive mode data extends into mass range of negative mode + 1000, increase mz_offset

mz_offset <- 1000

##ppd = positive peak data
##npd = negative peak data
##pnpd = pos and neg peak data

ppd_df <- arrange(cbind(pData(ppd),t(iData(ppd))),x,y)
npd_df <- arrange(cbind(pData(npd),t(iData(npd))),x,y)
pn_df <- cbind(ppd_df[1:ncol(ppd_df)],npd_df[4:ncol(npd_df)])

##make new MSIageSet:

pnpd <- MSImageSet(spectra=t(pn_df[4:ncol(pn_df)]), coord=pn_df[c("x","y","sample")], mz=c(mz(ppd),mz(npd)+mz_offset),
                   processingData = processingData(ppd),
                   protocolData = protocolData(ppd),
                   experimentData = experimentData(ppd))

##get rid of processing data.frames:

rm(ppd_df,npd_df,pn_df)
