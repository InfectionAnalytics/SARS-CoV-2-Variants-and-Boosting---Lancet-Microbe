library(reshape2)
library(ggplot2)
library(dplyr)
library(stringr)
library(lubridate)
library(zoo)
library(cowplot)
library(ggpubr)
library(latex2exp)
library(doBy)
library(RColorBrewer)
library(ggnewscale)
library(stats)

codeFolder = './'
plotFolder = './'
source('./dropPlot.R')
source('./plotValueByFactor.R')
source('./getCensoredMeans.R')
source('./censoredRegression.R')
source('./runRegressionsForAllVariants_VariantNeutL.R')

pangoNames=c('WT','B.1.1.7','B.1.351','P.1','B.1.617.2','B.1.617.1')
WHONames=c('Ancestral','Alpha','Beta','Gamma','Delta','Kappa')
serumCommonNames=c('Convalescent',	'Pfizer',	'AstraZeneca',	'Moderna',	'Novavax',	'Covaxin')
serumTechnicalNames=c('Convalescent',	'BNT162b2',	'ChAdOx1 nCoV-19',	'mRNA-1273',	'NVX-CoV2373',	'BBV152')

# Other global params
# Columns we can drop
excelOrigin<-"1899-12-30"
options(warn=-1)
options(dplyr.summarise.inform=F)

# For plotting
# Load the global plotting parameters
datestr=stringr::str_replace_all(Sys.Date(), c("-" = ""))
# Plotting parameters
defaultWidth=5
defaultHeight=6
colourPalette='Paired'
colourPalette='Set3'
laboratoryPallette = 'Set1'
serumPallette = 'Set2'
labColN=brewer.pal.info[laboratoryPallette,]$maxcolors
labMaxNo=11
serumColN=brewer.pal.info[laboratoryPallette,]$maxcolors
laboratoryColours <- colorRampPalette(brewer.pal(labColN,laboratoryPallette))(labMaxNo)
serumColours<-colorRampPalette(brewer.pal(serumColN,serumPallette))(8)

# This position function makes the width of columns the same
constWidthPosition=position_dodge2(preserve = "single", padding = .1)
doPlots=T
maxPlots=20
sizeFrame=data.frame(width=rep(defaultWidth,maxPlots),height=rep(defaultHeight,maxPlots))


col1='blue'
col2='grey97'
col3='grey94'
ypos=-.7
h=2*(1-ypos)+.1
blueWidth=.3
avgWidth=.2
shapeAlpha=.5
zeroWidth=.5
shapeSize=.4
shapes=c(16,4)

# Themes
VItheme<-theme_bw()+theme(legend.title = element_text(size=6, face='bold'), #change legend title font size
                          legend.text = element_text(size=6), #change legend text font size)
                          #plot.title = element_text(size=9, hjust=.5, face='bold'), #change title text font size)
                          #axis.title = element_text(size=8), #change title text font size)
                          axis.text = element_text(size=7), #change axis text font size)
                          axis.title = element_text(size=8), #change axis text font size)
                          legend.spacing.y=unit(2,"mm"),
                          legend.key.size = unit(2,"mm"),
                          legend.key.width = unit(4,"mm"),
                          legend.key.height = unit(3,"mm"),
                          #legend.box="vertical",
                          strip.background = element_rect(fill='azure'),
                          strip.text = element_text(size=5,face='bold',),
                          plot.title = element_text(size=9,face='bold',hjust=.5)
                          
                          
)
theme_anglexaxis<-theme(axis.text.x = element_text(angle = 45, hjust = 1))
theme_set(VItheme)

mapVariant<-function(thisPangoName='WT', includePango=T){
  mappedName=WHONames[match(thisPangoName,pangoNames)]
  if (includePango){
    mappedName=paste0(mappedName,' (',thisPangoName,')')
  }
  mappedName
}

mapSerum<-function(serumCommonName='WT', includeCommon=F){
  mappedName=serumTechnicalNames[match(serumCommonName,serumCommonNames)]
  if (includeCommon){
    mappedName=paste0(mappedName,' (',serumCommonName,')')
  }
  mappedName
}

mapLaboratory<-function(labName){
  paste('Lab',strtrim(labName,2))
}
format_VariantDisplay <- function() {
  function(x) mapVariant(x) 
}
format_SerumDisplay <- function() {
  function(x) mapSerum(x) 
}
format_LaboratoryDisplay <- function() {
  function(x) mapLaboratory(x) 
}

format_fraction <- function() {
  function(x)  ifelse(x>=1,as.character(round(x,0)),paste0('1/',as.character(round(1/x,0))))
}

format_fold <- function() {
  function(x)  ifelse(x>=1,as.character(round(x,0)),substr(as.character(round(x,3)),2,100))
}
