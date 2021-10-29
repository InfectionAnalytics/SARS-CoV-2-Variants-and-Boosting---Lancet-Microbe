makeDropPlot <-function(data=WTdata,xVar='WTneutL',yVar='variantneutL', cf='Conv', byLab=F, useWTdenom=F,shiftdata=NULL, rotated=F){
  #plotData = select(data,c('Serum','Laboratory','Variant','PaperRef'))
  ending='_Ser_Var'
  summaryVars=c('Serum','Variant')
  plotNameEnding = ''
  
  pointShape=1
  if (byLab) {
    ending = '_Lab_Ser_Var'
    summaryVars=c('Laboratory',summaryVars)
    plotNameEnding = 'WithLab'
    pointShape=21
  }
  if(is.null(shiftdata)){
    shiftdata=data
  }
  
  plotDataX=getCensoredMeans(data,summaryVars,paste0(xVar,'ChangeFrom',cf),censoredBelow = data$censoredWT)
  plotData = select(plotDataX,c(summaryVars, 'nstudy',paste0('muL_',xVar,'ChangeFrom',cf,ending),paste0('seL_',xVar,'ChangeFrom',cf,ending)))
  colnames(plotData)[length(summaryVars)+1+c(1,2)]=c('xmu','xse')
  if(byLab){
    plotData$fillVar=plotData$Laboratory
  } else {
    plotData$fillVar = plotData$Serum
  }
  
  
  if (!useWTdenom){
    yvar= paste0(yVar,'ChangeFrom',cf)
    ycens = data$censoredWT
  } else {
    data$variantneutLChangeFromConvWT = data$variantneutL-data$muL_WTneutL_Convalescent_Lab_Var
    yvar = 'variantneutLChangeFromConvWT'
    ycens = data$censoredMutant
  }
  plotDataY=getCensoredMeans(data,summaryVars,yvar,censoredBelow = ycens)
  plotDataY = select(plotDataY,c(summaryVars, 'nstudy',paste0('muL_',yvar,ending),paste0('seL_',yvar,ending)))
  colnames(plotDataY)[length(summaryVars)+1+c(1,2)]=c('ymu','yse')
  plotData = merge(plotData,plotDataY, by=summaryVars)
  plotData$nstudyMatch = plotData$nstudy.x==plotData$nstudy.y
  plotData = rename(plotData[,!(colnames(plotData)=='nstudy.y')], nstudy=nstudy.x)
  
  #plotData$xmu<-data[,which(colnames(data)==paste0('muL_',xVar,'ChangeFrom',cf,ending))]
  #plotData$ymu<-data[,which(colnames(data)==paste0('muL_',yVar,'ChangeFrom',cf,ending))]
  #plotData$xse<-data[,which(colnames(data)==paste0('seL_',xVar,'ChangeFrom',cf,ending))]
  #plotData$yse<-data[,which(colnames(data)==paste0('seL_',yVar,'ChangeFrom',cf,ending))]
  
  #plotData=mergeMeans(plotData,meanVariantChange,summaryVars)
  
  if(useWTdenom){
  #  plotData$ymu<-plotData[,paste0('muL_variantneutLChangeFromConvWT',ending)]
  #  plotData$yse<-plotData[,paste0('seL_variantneutLChangeFromConvWT',ending)]
    plotNameEnding = paste0(plotNameEnding,'WTdenom')
  }
  scaleBreaks=c(1/32,1/16,.125,.25,.5,1,2,4,8)
  minLine=(min(min(plotData$ymu-plotData$yse,na.rm=T),min(plotData$xmu-plotData$xse,na.rm=T)))
  maxLine=(max(max(plotData$ymu+plotData$yse,na.rm=T),max(plotData$xmu+plotData$xse,na.rm=T)))
  
  minLine2=(min(min(plotData$ymu-plotData$yse-plotData$shift,na.rm=T),min(plotData$xmu-plotData$xse,na.rm=T)))
  maxLine2=(max(max(plotData$ymu+plotData$yse-plotData$shift,na.rm=T),max(plotData$xmu+plotData$xse,na.rm=T)))
  
  
  # Add in the number of studies contributing
  #plotData<-merge(plotData, plotData %>% group_by(Serum,Variant) %>% summarise(nstudy=length(unique(PaperRef))),by=c('Serum','Variant'))
  # Convert to fold changes
  
  plotData[,c('xmu','xse','ymu','yse')] = 10^plotData[,c('xmu','xse','ymu','yse')]
  minLine=10^minLine
  maxLine=10^maxLine
  minLine2=10^minLine2
  maxLine2=10^maxLine2
  minLine2=.5
  maxLine2=4
  #plotData=distinct(plotData[,c("Serum","Variant","Laboratory","PaperRef","fillVar","xmu","ymu","xse","yse","nstudy")])
  
  #Add in a shift for the line
  if(useWTdenom){
    shifts = getCensoredMeans(shiftdata, 'Variant', 'changeFromWTL',censoredBelow = shiftdata$censoredMutant, censoredAbove= shiftdata$censoredWT)
    shifts = select(shifts,c('Variant','muL_changeFromWTL_Var','lower95L','upper95L'))
    shifts<-rename(shifts,shift=muL_changeFromWTL_Var)
    shifts[,c(2:ncol(shifts))]=10^shifts[,c(2:ncol(shifts))]
    plotData<-merge(plotData,shifts, by='Variant')
  } else {
    plotData$shift=1
  }
  
  lineX = c(seq(.06,.4,.01),seq(.4,4,.1))
  lineData = data.frame(lineX=lineX,Variant=shifts[1,'Variant'], shift=shifts[1,'shift'], L=shifts[1,'lower95L'],U=shifts[1,'upper95L'])
  for (i in c(2:nrow(shifts))){
    tempdf = data.frame(lineX=lineX,Variant=shifts[i,'Variant'], shift=shifts[i,'shift'], L=shifts[i,'lower95L'],U=shifts[i,'upper95L'])
    lineData = rbind(lineData, tempdf)
  }
    
  lw = .2
  line1col='blue'
  line2col='red'
  fillCol='red'
  alphaVal=.3
  
  titleLab = paste('Fold change in',yVar,'against',xVar,'and variant strains compared to',cf,'neut titre')
  yLab= expression('Neutralisation titre to variant (/convalescent to ancestral in same study)')
  xLab=expression('Neutralisation titre to ancestral virus: (/ convalescent to ancestral in same study)')
  if(rotated){xLab=expression('Neutralisation to Ancestral')}#: \n (/ Conv to Ancestral)')}
  #xlab(expression(atop("A long string of text for the purpose", paste("of illustrating my point" [reported]))))
  if(useWTdenom){
    yLab= expression('Neutralisation titre to variant: \n (/ convalescent to ancestral in same study)')
    if(rotated){
      yLab= expression('Neutralisation titre to variant: (/ convalescent to ancestral in same study)')
    }
  }
  titleLab=NULL
  variantLevels=paste0(paste(WHONames,pangoNames,sep=' ('),')')
  plotData$Variant2 <-mapVariant(plotData$Variant)
  plotData$Variant2<-factor(plotData$Variant2,levels=variantLevels)
  lineData$Variant2<-mapVariant(lineData$Variant)
  lineData$Variant2<-factor(lineData$Variant2,levels=variantLevels)
  dropPlot=(ggplot(plotData,aes(x=xmu,y=ymu,colour=Serum, fill=fillVar))
  +labs(title = titleLab, x=xLab, y=yLab)
  +theme(axis.title.x = element_text(vjust=1, hjust = .5))
  )
  
  dropPlot=(dropPlot
            #+geom_point(aes(fill=fillVar))
            +scale_colour_manual(values=serumColours, labels=format_SerumDisplay(), name= 'Serum (Panel B)')
            +geom_errorbar(aes(ymin=ymu/yse,ymax=ymu*yse))
            +geom_errorbarh(aes(xmin=xmu/xse,xmax=xmu*xse))
            +geom_line(data=lineData, aes(x=lineX,y=lineX, fill=NA),colour=line1col,size=lw, linetype='dashed')
            +geom_line(data=lineData, aes(x=lineX,y=shift*lineX, fill=NA),colour=line2col,size=lw, linetype='dashed')
            +geom_ribbon(data=lineData, aes(x=lineX,y=lineX,ymin=lineX*L,ymax=lineX*U, fill=fillCol),fill=fillCol,colour=NA,alpha=alphaVal, show.legend = F)
            +scale_x_continuous(breaks=scaleBreaks, trans='log10', labels = format_fold(),limits=c(max(minLine,min(scaleBreaks)),min(maxLine,max(scaleBreaks))))
            +scale_y_continuous(breaks=scaleBreaks, trans='log10', labels = format_fold(),limits=c(max(minLine,min(scaleBreaks)),min(maxLine,max(scaleBreaks))))
  )
  if (byLab){
    dropPlot = dropPlot 
  } else {
    dropPlot = dropPlot +scale_fill_manual(values=serumColours, guide=F, name = 'Serum (Panel B)')
  }
  if (!byLab){
    dropPlot = dropPlot+geom_text(aes(label=paste0('(',nstudy,')'), guide=F),guide=F,nudge_y=.2,nudge_x=-.05, size=2,fontface='bold' , show.legend=F)
  }
  dropPlot=dropPlot+geom_point(shape=pointShape)
  
  # Add in the vertical arrows
  arrowShift = 1
  if (byLab){
    arrowShift=2
  }
  if(useWTdenom){
    dropPlot = (dropPlot
                +geom_segment(aes(x = 1/arrowShift, xend = 1/arrowShift, y = .9/arrowShift, yend = 1.2*shift/arrowShift),colour = "purple",size=.3, arrow = arrow(length = unit(0.03, "npc"))) 
                +geom_text(aes(x = .25/arrowShift,y=1/arrowShift,label='Average drop in IC50'),colour='purple',size=1.5)
                +geom_text(aes(x = .25/arrowShift,y=.7/arrowShift,label=paste0('for ',Variant,' (',round(1/shift,1),')')),colour='purple',size=1.5))
  }
  if (rotated){ 
    facet_rows = length(unique(dropPlot$data$Variant2))
    facet_cols = 1
  } else {
    facet_rows = 1
    facet_cols = length(unique(dropPlot$data$Variant2))
  }
  dropPlot = dropPlot+facet_wrap(~Variant2, nrow = facet_rows, ncol=facet_cols)
  
  if(byLab){
    dropPlot=(dropPlot+scale_fill_manual(name='Laboratory',values=laboratoryColours,breaks=levels(plotData$fillVar), labels=format_LaboratoryDisplay())
                    +theme(axis.text.x = element_text(angle=45,hjust=1,vjust=.5), axis.title.x=element_text(vjust=1)) )
  }
  
  if (rotated){
    plotWidth=3
    plotHeight=7
  } else{
    plotWidth=8
    plotHeight=2
  }
  plotName=paste0('dropPlot_',yVar,'vs',xVar,'cf',cf,plotNameEnding,'_',datestr)
  dropPlotName=paste0(plotFolder,plotName,".pdf")
  list(name=dropPlotName,plot=dropPlot)
  
}