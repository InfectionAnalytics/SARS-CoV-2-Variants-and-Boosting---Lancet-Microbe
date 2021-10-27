findXVals<-function(xAxisGroup='Serum',dataset=WTdata, nameStart='',BarWidthTotal=.65){
  if (xAxisGroup=='Serum'){
    colourGroup='Laboratory'
  } else if (xAxisGroup=='Laboratory'){
    colourGroup='Serum'
  }
  dataset$Laboratory <-factor(dataset$Laboratory)
  
  alreadyPresent<-distinct(dataset[colnames(dataset)%in%c(colourGroup,xAxisGroup,'Variant')])
  # Number of colour groups
  ncolours<-alreadyPresent %>% group_by(.dots=c(xAxisGroup,'Variant')) %>% summarise(n=n())
  widthInfp<-ncolours%>% group_by(Variant)%>% summarise(largestGroup=max(n))
  width<-min(BarWidthTotal/(widthInfp$largestGroup-1))
  
  groupInfo<-alreadyPresent %>% group_by(.dots=c('Variant',xAxisGroup)) %>% summarise_at(colourGroup,list(colourGrouping=identity,rank=rank))
  colnames(groupInfo)[3]=colourGroup
  groupInfo<-groupInfo[order(groupInfo$Variant,unlist(groupInfo[,colnames(groupInfo)==xAxisGroup]),groupInfo$rank),]
  
  groupInfo$offset<-(groupInfo$rank-1)*width
  medians<-groupInfo %>% group_by(.dots=c('Variant',xAxisGroup)) %>% summarise(median=median(offset))
  groupInfo=merge(groupInfo,medians,by=c('Variant',xAxisGroup),all.x=T)
  groupInfo$totalOffset=groupInfo$offset+unclass(unlist(groupInfo[,colnames(groupInfo)==xAxisGroup]))-groupInfo$median
  xoffsets<-subset(groupInfo,select=c('Variant',xAxisGroup,colourGroup,'totalOffset'))
  colnames(xoffsets)[4]<-paste0(nameStart,colnames(xoffsets)[4],xAxisGroup)
  xoffsets
}
