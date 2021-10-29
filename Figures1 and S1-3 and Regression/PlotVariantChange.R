# Add in ancestral
censoredChangeMeans_Var=getCensoredMeans(WTdata,c('Variant'),'changeFromWTL',censoredBelow = WTdata$censoredMutant, censoredAbove= WTdata$censoredWT)
plotData=censoredChangeMeans_Var
plotData1=plotData[1,]
plotData1$Variant='WT'
plotData1[,-c(1)]=0
plotData=rbind(plotData1,plotData)

scaleBreakLabels=c(1/16,1/8,1/4,1/2,1)
scaleBreaks=log10(scaleBreakLabels)

cuts = 2^c(1:4)
foldCuts = sort(c(cuts,1,1/cuts),decreasing=T)
log10cuts =log10(foldCuts)
foldLabels=as.character(foldCuts)

foldLabels[foldCuts<1]=paste0('',as.character(1/foldCuts[foldCuts<1]))


labelSize=2
variantChangePlot=(ggplot(data=plotData, aes(x=factor(Variant,levels = pangoNames),y=muL_changeFromWTL_Var))
  +geom_text(aes(label=paste0('(',round(10^-muL_changeFromWTL_Var,1),')')),nudge_x=-.42, show.legend = F,alpha=1,size=labelSize+1)
  +geom_rect(xmin=0,xmax=.95,ymin=-.1,ymax=.1,fill='white')
  +geom_point(aes(alpha=n),shape=16)
  +geom_errorbar(aes(ymin=lower95L,ymax=upper95L,width=.5,alpha=n),size=1)
  +geom_text(aes(label=paste0('samples=',n))     ,nudge_y=.1, show.legend = F,alpha=1,size=labelSize)
  +geom_text(aes(label=paste0('studies=',nstudy)),nudge_y=.15, show.legend = F,alpha=1,size=labelSize)
  +labs(title='Fold Drop in neutralisation IC50 compared to Ancestral Virus',x='Variant',y='Fold Drop in IC50 compared to Ancestral Virus')
  +scale_y_continuous(breaks=log10cuts,labels=foldLabels, limits=log10(c(1/16,1)))
  +scale_x_discrete(labels=format_VariantDisplay(), breaks=pangoNames)
  +scale_alpha_continuous(range=c(.4,1),guide=F)
  +scale_size_continuous(guide=F)
  +geom_hline(mapping=aes(yintercept=0),colour='black',size=.5,linetype='solid')
  +theme(legend.text = element_text(size=4),
       legend.spacing.y=unit(1,"mm"),
       legend.key.size = unit(2,"mm"),
       legend.key.width = unit(2,"mm"),
       legend.key.height = unit(2,"mm"),
       plot.title = element_text(size=8,face='bold',hjust=.5))
  +theme_anglexaxis
)

plotWidth=4
plotHeight=4
ggsave('SuppFigure_3.pdf',width=plotWidth, height=plotHeight)

  