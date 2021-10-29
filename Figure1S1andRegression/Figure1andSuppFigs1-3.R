source('./loadGlobalVariantParameters.R')
source('./loadVariantDataFrame.R')
source('./getCensoredMeans.R')
source('./censoredRegression.R')
source('./regressionSummaryHelperFunctions.R')
source('./findXVals.R')
source('./dropPlot.R')
source('./plotValueByFactor.R')
source('./runRegressionsForAllVariants_VariantNeutL.R')

WTdata =loadVariantDataFrame() 
#### Box Plots ####
boxPlotVals = c('changeFromWTL',
                     'changeFromWTLcfConv',
                     'changeFromWTLcfPfizer'
)

boxPlotNames=c(
  'Reduction in neutralization titre (compared to Ancestral)',
  'Difference in reduction in neutralization titre between vaccine and convalescent serum (in same assay)',
  'Difference in reduction in neutralization titre between vaccine / convalescent serum and BNT162b2 (in same assay)'
)

fig1a=plotValueByFactor(dataset=WTdata, value=boxPlotVals[1], by='Laboratory',valueName=boxPlotNames[1], titleVal='', censoredCol='censoredMutant')
fig1b=makeDropPlot(data=WTdata,xVar='WTneutL',yVar='variantneutL', cf='Conv', byLab=F,useWTdenom = T, rotated=T)

fig1<-plot_grid(
  fig1a$plot+ theme(legend.position="none", axis.title.x=element_text(hjust=.5,vjust=1)),
  fig1b$plot+theme(legend.position="nonw", axis.title.x = element_text(vjust = 0, hjust=.5)),
  align = 'h',
  labels = c("A", "B"),
  hjust = -1,
  nrow = 1,
  rel_widths = c(1.3, 1)
)

legend1 <- get_legend(fig1a$plot + theme(legend.box.margin = margin(0, 0, 0, 0)))
legend2 <- get_legend(fig1b$plot + theme(legend.box.margin = margin(0, 0, 0, 0)))
legends<-plot_grid(NULL,legend1,legend2,NULL,ncol = 1,align='h',axis='t')
fig1_full<-plot_grid(fig1,legends,
                   align = 'vh',
                   rel_widths = c(4, 1),
                   ncol=2
)
ggsave('./Figure_1.pdf',fig1_full,width=7, height=10) 


supp_plots=list()
for (i in c(1:2)){
  plotValue = boxPlotVals[i]
  plotValName = boxPlotNames[i]
  supp_plots[[i]]=plotValueByFactor(dataset=WTdata, value=plotValue, by='Laboratory',valueName=plotValName, titleVal='', censoredCol='censoredMutant')
}  

# Now do a plot coloured by Serum - this is for supp materials
supp_plots[[3]]=plotValueByFactor(dataset=WTdata, value='changeFromWTL', by='Serum',valueName='Reduction in neutralization titre (compared to Ancestral)', titleVal='', censoredCol='censoredMutant')

tmp=supp_plots[[2]]
supp_plots[[2]]=supp_plots[[3]]
supp_plots[[3]]=tmp

supFig1<-plot_grid(supp_plots[[1]]$plot+ theme(legend.position="none", axis.title.x = element_text(vjust=0)),
               supp_plots[[2]]$plot+ theme(legend.position="none", axis.title.x = element_text(vjust=0)),
               supp_plots[[3]]$plot+ theme(legend.position="none", axis.title.x = element_text(vjust=0)),
               align = 'vh',
               #rel_widths = c(1, 1,1),
               labels = c("A", "B","C"),
               hjust = -1,
               nrow=1
)
#ggsave(paste0(manuscriptPlotFolder,'SuppFigure_1a.png'),fs1,width=12, height=10) 

legend1 <- get_legend(supp_plots[[1]]$plot +guides(col=guide_legend(title='Laboratory \n(Panels A and C)',title.hjust=.5))+ theme(legend.box.margin = margin(0, 0, 0, 0)))
legend2 <- get_legend(supp_plots[[2]]$plot +guides(col=guide_legend(title='Serum (Panel B)',title.hjust=.5))+ theme(legend.box.margin = margin(0, 0, 0, 0)))
legends<-plot_grid(NULL,legend1,legend2,NULL,ncol = 1,align='h',axis='t')

supFig1_full<-plot_grid(supFig1,legends,
                    align = 'v',
                    rel_widths = c(5, .8),
                    rel_heights = c(3, 1),
                    ncols=2
)
ggsave('./SuppFigure_1.pdf',supFig1_full,width=13, height=10)  


# Make Supplementary Figure 2
figS2=makeDropPlot(data=WTdata,xVar='WTneutL',yVar='variantneutL', cf='Conv', byLab=T,useWTdenom = T, rotated=T)
ggsave('./SuppFigure_2.pdf',figS2$plot, width=12, height=10) 

# Plot Supplementary Figure 3
source('./PlotVariantChange.R')


