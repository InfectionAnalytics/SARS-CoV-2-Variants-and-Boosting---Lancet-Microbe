plotValueByFactor<-function(dataset=WTdata, value='WTneutL', by='Laboratory',xby=NULL, valueName=NULL, titleVal=NULL, censoredCol=NULL){
  if (is.null(xby)){  
    if (by=='Laboratory'){
      xby='Serum'
    } else if (by=='Serum'){
      xby='Laboratory'
    }
  } else if (xby=='Laboratory'){
    by='Serum'
  }
  if (is.null(valueName)){
    valueName=value
  }
  
  cuts = 4^c(1:3)
  foldCuts = sort(c(cuts,1,1/cuts),decreasing=F)
  log10cuts =log10(foldCuts)
  foldLabels=as.character(rev(foldCuts))
  #foldLabels[foldCuts<1]=paste0('1/',as.character(1/foldCuts[foldCuts<1]))
  
  dataset=dataset[!is.na(dataset[,xby]),]
  
  plotName=paste0(value,'_ColouredBy',by)
  
  dataset$groupMeansY=dataset[,paste('muL',value,substr(xby,1,3),'Var',sep='_')]
  dataset$allMeanY=dataset[,paste('muL',value,'Var',sep='_')]
  dataset$yValue = dataset[,value]
  if (is.null(censoredCol)){
    dataset$censored=F
  } else{
    dataset$censored=dataset[,censoredCol]
}
  plotFileNameCurrent=paste0(plotFolder,plotName,datestr,".png")
  
  dataset$colourGroup = dataset[,by]
  eval(parse(text=paste0('colourValues=',tolower(by),'Colours')))
  dataset$groupMeansX=unclass(dataset[,xby])
  dataset$allMeanX=dataset[,xby]
  dataset$totalOffset=dataset[,paste0('changeFromWTLcfConvtotalOffset',xby)]
  xshift=0
  leftcol=col2
  if (xby=='Serum'){
    leftcol=col1
    leftcol='lightblue1'
  }
  dataset$xValue=dataset[,xby]
  if (xby=='Serum'){
    dataset$xValue=factor(dataset$xValue, serumCommonNames)
    dataset$colourGroup=mapLaboratory(dataset$colourGroup)
  } else if (xby=='Laboratory'){
    #dataset$xValue=mapLaboratory(dataset$colourGroup)
    dataset$colourGroup=factor(dataset$colourGroup, serumCommonNames, serumTechnicalNames)
  }
  
  dataset<-dataset[,c('Variant','Serum','Laboratory','xValue','yValue','colourGroup','groupMeansY','allMeanY','groupMeansX','allMeanX','totalOffset','censored')]
  variantLevels=paste0(paste(WHONames,pangoNames,sep=' ('),')')
  dataset$Variant <-mapVariant(dataset$Variant)
  dataset$Variant<-factor(dataset$Variant,levels=variantLevels)
  
  # replace the NA values so that something appears
  #dataset$xValue[is.na(dataset$xValue)]<-dataset$allMeanX[is.na(dataset$xValue)]
  #dataset$yValue[is.na(dataset$yValue)]<-dataset$allMeanY[is.na(dataset$yValue)]
  dataset$xValue[is.na(dataset$xValue)]<-0
  dataset$yValue[is.na(dataset$yValue)]<-0
  
  cfConvalescentBoxPlot<-(ggplot(dataset,aes(x=xValue,y=yValue,colour=(colourGroup)))
                            +geom_boxplot(position=position_dodge2(preserve='single',padding=.3),outlier.shape=NA, mapping=aes(width=2),inherit.aes = TRUE)
                            +facet_wrap(~Variant,ncol=1,scales='fixed')
  )
  yGap=.5
  midptY=yGap*round(median(dataset$yValue,na.rm=T)/yGap)
  maxPos=yGap*ceiling(max(dataset$yValue,na.rm=T)/yGap)
  maxNeg=yGap*floor(min(dataset$yValue,na.rm=T)/yGap)
  if (abs(maxPos)>=abs(maxNeg)){
    max=maxPos
  } else{
    diff=midptY-maxNeg
    max=midptY+diff
  }
  h2=2*(max-midptY)+.1
  
  alphaVal=.2
  # Now add in the backgrounds
  convcol=leftcol
  secondcol=col3
  if(endsWith(value,'Pfizer')){
    convcol=col2
    secondcol=leftcol
  }
  cfConvalescentBoxPlot=(cfConvalescentBoxPlot
                         +geom_tile(mapping=aes(x=1,y=midptY,width=1,height=h2),fill=convcol,alpha=alphaVal, colour=convcol)
                         +geom_tile(mapping=aes(x=3,y=midptY,width=1,height=h2),fill=col2,alpha=alphaVal,colour=col2)
                         +geom_tile(mapping=aes(x=5,y=midptY,width=1,height=h2),fill=col2,alpha=alphaVal,colour=col2)
                         +geom_tile(mapping=aes(x=2,y=midptY,width=1,height=h2),fill=secondcol,alpha=alphaVal,colour=secondcol)
                         +geom_tile(mapping=aes(x=4,y=midptY,width=1,height=h2),fill=col3,alpha=alphaVal,colour=col3)
                         +geom_tile(mapping=aes(x=6,y=midptY,width=1,height=h2),fill=col3,alpha=alphaVal,colour=col3)
  )
  if (length(unique(dataset[,xby]))>6){
    cfConvalescentBoxPlot=cfConvalescentBoxPlot+geom_tile(mapping=aes(x=7,y=midptY,width=1,height=h2),fill=col2,alpha=alphaVal,colour=col2)
  }
  if(length(unique(dataset[,xby]))>7){
    cfConvalescentBoxPlot=cfConvalescentBoxPlot+geom_tile(mapping=aes(x=8,y=midptY,width=1,height=h2),fill=col3,alpha=alphaVal,colour=col3)
  }
  if(length(unique(dataset[,xby]))>8) {
    cfConvalescentBoxPlot=cfConvalescentBoxPlot+geom_tile(mapping=aes(x=9,y=midptY,width=1,height=h2),fill=col2,alpha=alphaVal,colour=col2)
  }
  if(length(unique(dataset[,xby]))>9) {
    cfConvalescentBoxPlot=cfConvalescentBoxPlot+geom_tile(mapping=aes(x=10,y=midptY,width=1,height=h2),fill=col3,alpha=alphaVal,colour=col3)
  }
  if(length(unique(dataset[,xby]))>10) {
    cfConvalescentBoxPlot=cfConvalescentBoxPlot+geom_tile(mapping=aes(x=11,y=midptY,width=1,height=h2),fill=col2,alpha=alphaVal,colour=col2)
  }
  xlims=c(min(unclass(dataset$xValue)-xshift),max(unclass(dataset$xValue)-xshift))
  
  if (is.null(titleVal)){
    titleStr =  TeX(paste0(valueName,' coloured by ',by,''))
  } else {
    titleStr =  TeX(titleVal)
  }
  valueStr = TeX(valueName)
  
  format_XLabels<-function(){
    ifelse(xby=='Serum',format_SerumDisplay(),format_LaboratoryDisplay()) 
  }
  
  cfConvalescentBoxPlot = (cfConvalescentBoxPlot+geom_boxplot(position=position_dodge2(preserve='single',padding=.3),outlier.shape=NA, mapping=aes(width=2),inherit.aes = TRUE)
                                      +geom_point(mapping=aes(x=totalOffset-xshift, alpha=shapeAlpha,shape=censored),size=shapeSize, position=position_jitter())
                                      +geom_hline(width=1,yintercept=0,colour='black',size=zeroWidth,show.legend = T)
                                      +scale_colour_manual(values=colourValues, name=paste0(by,' (Panel A)'))
                                      +scale_alpha_identity(guide=F)
                                      +scale_shape_manual(values=shapes,guide=F)
                                      +geom_segment(mapping=aes(x=unclass(dataset$xValue)-xshift-.25,xend=unclass(dataset$xValue)-xshift+.25,y=groupMeansY,yend=groupMeansY),colour='blue',size=blueWidth,linetype='solid')
                                      +geom_hline(mapping=aes(x=allMeanX,yintercept=allMeanY),colour='black',size=avgWidth,linetype='dotdash')
                                      +theme_anglexaxis 
                                      #+geom_text(aes(y=log10(1/64), label = paste0('Change=',round(10^groupMeansY,2),'fold')))
                                      +scale_y_continuous(expand=c(0,0), breaks=log10cuts,labels=foldLabels,lim=c(midptY-h2/2,midptY+h2/2))
                                      +scale_x_discrete(expand=c(0,0),na.translate=F, labels=format_XLabels())
                                      #+annotate('text',label='Loss of Neut in Serum > Loss in Convalescent',x=3,y=-1.3,size=1.5,colour='red')
                                      #+annotate('text',label='Loss of Neut in Serum < Loss in Convalescent',x=4,y=.8,size=1,colour='darkgreen')
                                      #+labs(y=TeX(r'(log_{10}(Fold Drop) Serum:Convalescent)'),title=paste0(plotNameStart,' by ',by,' :: Fold Drop compared to Conv',x=''))
                                      +labs(y=valueStr,title=titleStr,x=xby)
  )
  
  plotWidth=4.25
  plotHeight=7

  list(name=plotFileNameCurrent,plot=cfConvalescentBoxPlot)
}
