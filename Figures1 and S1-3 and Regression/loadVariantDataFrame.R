baseFolder<-'./'
codeFolder<-'./'

#Load Data
loadVariantDataFrame <- function() {
WTdata = read.csv('./variantNeutDataForImport.csv')
WTdata$Serum<-factor(WTdata$Serum,levels=serumCommonNames,labels=serumCommonNames,ordered=F)
WTdata$Laboratory<-factor(WTdata$Laboratory) 

# This is only for plotting
xoffsetsForVaccine=findXVals('Serum',WTdata,'drop')
xoffsetsForLaboratory=findXVals('Laboratory',WTdata,'drop')

# This uses the censored means of the convalescent data
censoredChangeMeans_LabSerumVar=getCensoredMeans(WTdata,c('Laboratory','Serum','Variant'),'changeFromWTL',censoredBelow = WTdata$censoredMutant)
censoredChangeMeans_SerumVar=getCensoredMeans(WTdata,c('Serum','Variant'),'changeFromWTL',censoredBelow = WTdata$censoredMutant, censoredAbove= WTdata$censoredWT)
censoredChangeMeans_LabVar=getCensoredMeans(WTdata,c('Laboratory','Variant'),'changeFromWTL',censoredBelow = WTdata$censoredMutant, censoredAbove= WTdata$censoredWT)
censoredChangeMeans_Var=getCensoredMeans(WTdata,c('Variant'),'changeFromWTL',censoredBelow = WTdata$censoredMutant, censoredAbove= WTdata$censoredWT)

censoredChangeB117Means_LabSerumVar=getCensoredMeans(WTdata,c('Laboratory','Serum','Variant'),'changeFromB117L',censoredBelow = WTdata$censoredMutant)
censoredChangeB117Means_SerumVar=getCensoredMeans(WTdata,c('Serum','Variant'),'changeFromB117L',censoredBelow = WTdata$censoredMutant, censoredAbove= WTdata$censoredB117)
censoredChangeB117Means_LabVar=getCensoredMeans(WTdata,c('Laboratory','Variant'),'changeFromB117L',censoredBelow = WTdata$censoredMutant, censoredAbove= WTdata$censoredB117)
censoredChangeB117Means_Var=getCensoredMeans(WTdata,c('Variant'),'changeFromB117L',censoredBelow = WTdata$censoredMutant, censoredAbove= WTdata$censoredB117)

censoredWTneutLMeans_LabSerumVar=getCensoredMeans(WTdata,c('Laboratory','Serum','Variant'),'WTneutL',censoredBelow = WTdata$censoredWT)
censoredWTneutLMeans_SerumVar=getCensoredMeans(WTdata,c('Serum','Variant'),'WTneutL',censoredBelow = WTdata$censoredWT)
censoredWTneutLMeans_LabVar=getCensoredMeans(WTdata,c('Laboratory','Variant'),'WTneutL',censoredBelow = WTdata$censoredWT)
censoredWTneutLMeans_Var=getCensoredMeans(WTdata,c('Variant'),'WTneutL',censoredBelow = WTdata$censoredWT)

censoredB117neutLMeans_LabSerumVar=getCensoredMeans(WTdata,c('Laboratory','Serum','Variant'),'B117neutL',censoredBelow = WTdata$censoredB117)
censoredB117neutLMeans_SerumVar=getCensoredMeans(WTdata,c('Serum','Variant'),'B117neutL',censoredBelow = WTdata$censoredB117)
censoredB117neutLMeans_LabVar=getCensoredMeans(WTdata,c('Laboratory','Variant'),'B117neutL',censoredBelow = WTdata$censoredB117)
censoredB117neutLMeans_Var=getCensoredMeans(WTdata,c('Variant'),'B117neutL',censoredBelow = WTdata$censoredB117)

censoredVariantneutLMeans_LabSerumVar=getCensoredMeans(WTdata,c('Laboratory','Serum','Variant'),'variantneutL',censoredBelow = WTdata$censoredMutant)
censoredVariantneutLMeans_SerumVar=getCensoredMeans(WTdata,c('Serum','Variant'),'variantneutL',censoredBelow = WTdata$censoredMutant)
censoredVariantneutLMeans_LabVar=getCensoredMeans(WTdata,c('Laboratory','Variant'),'variantneutL',censoredBelow = WTdata$censoredMutant)
censoredVariantneutLMeans_Var=getCensoredMeans(WTdata,c('Variant'),'variantneutL',censoredBelow = WTdata$censoredMutant)

#Obtain the mean log10fold change for mutant:WT by Serum, and Variant
#meansByVaccineVariant<-WTdata %>% group_by(Serum,Variant) %>% summarise (meanVaccineVariantlog10Change = mean(changeFromWTL,na.rm=T))
#Obtain the mean log10fold change for mutant:WT by Laboratory, and Variant
#meansByLaboratoryVariant<-WTdata %>% group_by(Laboratory,Variant) %>% summarise (meanLaboratoryVariantlog10Change = mean(changeFromWTL,na.rm=T))
#Obtain the mean log10fold change for mutant:WT by Variant
#meansByVariant<-WTdata %>% group_by(Variant) %>% summarise (meanVariantlog10Change = mean(changeFromWTL,na.rm=T))

mergeMeans=function(data,mergeData,byCols, includesig=T){
  useCols=startsWith(colnames(mergeData),'muL')|(colnames(mergeData) %in% byCols)
  if(includesig){
    useCols=useCols|startsWith(colnames(mergeData),'sigL')|startsWith(colnames(mergeData),'seL')
  }
  mergeData=mergeData[,useCols]
  data=merge(data,mergeData,by=byCols,all.x=T,all.y=T)
}
WTdata=mergeMeans(WTdata,censoredChangeMeans_Var,c('Variant'))
WTdata=mergeMeans(WTdata,censoredChangeMeans_SerumVar,c('Serum','Variant'))
WTdata=mergeMeans(WTdata,censoredChangeMeans_LabVar,c('Laboratory','Variant'))
WTdata=mergeMeans(WTdata,censoredChangeMeans_LabSerumVar,c('Laboratory','Serum','Variant'))

WTdata=mergeMeans(WTdata,censoredChangeB117Means_Var,c('Variant'))
WTdata=mergeMeans(WTdata,censoredChangeB117Means_SerumVar,c('Serum','Variant'))
WTdata=mergeMeans(WTdata,censoredChangeB117Means_LabVar,c('Laboratory','Variant'))
WTdata=mergeMeans(WTdata,censoredChangeB117Means_LabSerumVar,c('Laboratory','Serum','Variant'))
                          
WTdata=mergeMeans(WTdata,censoredWTneutLMeans_Var,c('Variant'))
WTdata=mergeMeans(WTdata,censoredWTneutLMeans_SerumVar,c('Serum','Variant'))
WTdata=mergeMeans(WTdata,censoredWTneutLMeans_LabVar,c('Laboratory','Variant'))

WTdata=mergeMeans(WTdata,censoredB117neutLMeans_Var,c('Variant'))
WTdata=mergeMeans(WTdata,censoredB117neutLMeans_SerumVar,c('Serum','Variant'))
WTdata=mergeMeans(WTdata,censoredB117neutLMeans_LabVar,c('Laboratory','Variant'))

WTdata=mergeMeans(WTdata,censoredVariantneutLMeans_Var,c('Variant'))
WTdata=mergeMeans(WTdata,censoredVariantneutLMeans_SerumVar,c('Serum','Variant'))
WTdata=mergeMeans(WTdata,censoredVariantneutLMeans_LabVar,c('Laboratory','Variant'))

WTdata=merge(WTdata,xoffsetsForVaccine, by=c('Serum','Variant','Laboratory'),all.x=T,all.y=T)
WTdata=merge(WTdata,xoffsetsForLaboratory, by=c('Laboratory','Variant','Serum'),all.x=T,all.y=T)

#### Add data relative to convalescent means ####
# Extract only the relevant convalescent Data
#convData=distinct(WTdata[WTdata$Serum=='Convalescent',])
convalescentWTneutLMeans=censoredWTneutLMeans_LabSerumVar[censoredWTneutLMeans_LabSerumVar$Serum=='Convalescent',]
colnames(convalescentWTneutLMeans)[str_detect(colnames(convalescentWTneutLMeans),'Lab_Ser_Var')]=str_replace(colnames(convalescentWTneutLMeans)[str_detect(colnames(convalescentWTneutLMeans),'Lab_Ser_Var')],'Lab_Ser_Var','Convalescent_Lab_Var')
convalescentWTneutLMeans=select(convalescentWTneutLMeans,c('Laboratory','Variant','muL_WTneutL_Convalescent_Lab_Var'))

convalescentB117neutLMeans=censoredB117neutLMeans_LabSerumVar[censoredB117neutLMeans_LabSerumVar$Serum=='Convalescent',]
colnames(convalescentB117neutLMeans)[str_detect(colnames(convalescentB117neutLMeans),'Lab_Ser_Var')]=str_replace(colnames(convalescentB117neutLMeans)[str_detect(colnames(convalescentB117neutLMeans),'Lab_Ser_Var')],'Lab_Ser_Var','Convalescent_Lab_Var')
convalescentB117neutLMeans=select(convalescentB117neutLMeans,c('Laboratory','Variant','muL_B117neutL_Convalescent_Lab_Var'))

convalescentVariantneutLMeans=censoredVariantneutLMeans_LabSerumVar[censoredVariantneutLMeans_LabSerumVar$Serum=='Convalescent',]
colnames(convalescentVariantneutLMeans)[str_detect(colnames(convalescentVariantneutLMeans),'Lab_Ser_Var')]=str_replace(colnames(convalescentVariantneutLMeans)[str_detect(colnames(convalescentVariantneutLMeans),'Lab_Ser_Var')],'Lab_Ser_Var','Convalescent_Lab_Var')
convalescentVariantneutLMeans=select(convalescentVariantneutLMeans,c('Laboratory','Variant','muL_variantneutL_Convalescent_Lab_Var'))

convalescentChangeMeans=censoredChangeMeans_LabSerumVar[censoredChangeMeans_LabSerumVar$Serum=='Convalescent',]
colnames(convalescentChangeMeans)[str_detect(colnames(convalescentChangeMeans),'Lab_Ser_Var')]=str_replace(colnames(convalescentChangeMeans)[str_detect(colnames(convalescentChangeMeans),'Lab_Ser_Var')],'Lab_Ser_Var','Convalescent_Lab_Var')

convalescentChangeB117Means=censoredChangeB117Means_LabSerumVar[censoredChangeB117Means_LabSerumVar$Serum=='Convalescent',]
colnames(convalescentChangeB117Means)[str_detect(colnames(convalescentChangeB117Means),'Lab_Ser_Var')]=str_replace(colnames(convalescentChangeB117Means)[str_detect(colnames(convalescentChangeB117Means),'Lab_Ser_Var')],'Lab_Ser_Var','Convalescent_Lab_Var')

# Merge convalescent data against  WT
WTdata=mergeMeans(WTdata,convalescentWTneutLMeans, by=c('Laboratory','Variant'),includesig=F)
WTdata$WTneutLChangeFromConv = WTdata$WTneutL-WTdata$muL_WTneutL_Convalescent_Lab_Var

# Merge convalescent data against  B117
WTdata=mergeMeans(WTdata,convalescentB117neutLMeans, by=c('Laboratory','Variant'),includesig=F)
WTdata$B117neutLChangeFromConv = WTdata$B117neutL-WTdata$muL_B117neutL_Convalescent_Lab_Var

#  Merge convalescent data against  variant
WTdata=mergeMeans(WTdata,convalescentVariantneutLMeans, by=c('Laboratory','Variant'),includesig=F)
WTdata$variantneutLChangeFromConv = WTdata$variantneutL-WTdata$muL_variantneutL_Convalescent_Lab_Var

# Merge change from WT to variant
WTdata=mergeMeans(WTdata,convalescentChangeMeans, by=c('Laboratory','Variant'),includesig=F)
WTdata$changeFromWTLcfConv = WTdata$changeFromWTL-WTdata$muL_changeFromWTL_Convalescent_Lab_Var

# Merge change from B117 to variant
WTdata=mergeMeans(WTdata,convalescentChangeB117Means, by=c('Laboratory','Variant'),includesig=F)
WTdata$changeFromB117LcfConv = WTdata$changeFromB117L-WTdata$muL_changeFromB117L_Convalescent_Lab_Var

# Now calculate mean drops compared to convalescent for WT
censoredWTneutLChangeFromConvMeans_LabSerumVar=getCensoredMeans(WTdata,c('Laboratory','Serum','Variant'),'WTneutLChangeFromConv',censoredBelow = WTdata$censoredWT)
censoredWTneutLChangeFromConvMeans_SerumVar=getCensoredMeans(WTdata,c('Serum','Variant'),'WTneutLChangeFromConv',censoredBelow = WTdata$censoredWT)
censoredWTneutLChangeFromConvMeans_LabVar=getCensoredMeans(WTdata,c('Laboratory','Variant'),'WTneutLChangeFromConv',censoredBelow = WTdata$censoredWT)
censoredWTneutLChangeFromConvMeans_Var=getCensoredMeans(WTdata,c('Variant'),'WTneutLChangeFromConv',censoredBelow = WTdata$censoredWT)

# Now calculate mean drops compared to convalescent for variant
censoredvariantneutLChangeFromConvMeans_LabSerumVar=getCensoredMeans(WTdata,c('Laboratory','Serum','Variant'),'variantneutLChangeFromConv',censoredBelow = WTdata$censoredMutant)
censoredvariantneutLChangeFromConvMeans_SerumVar=getCensoredMeans(WTdata,c('Serum','Variant'),'variantneutLChangeFromConv',censoredBelow = WTdata$censoredMutant)
censoredvariantneutLChangeFromConvMeans_LabVar=getCensoredMeans(WTdata,c('Laboratory','Variant'),'variantneutLChangeFromConv',censoredBelow = WTdata$censoredMutant)
censoredvariantneutLChangeFromConvMeans_Var=getCensoredMeans(WTdata,c('Variant'),'variantneutLChangeFromConv',censoredBelow = WTdata$censoredMutant)

# Now calculate mean drops compared to convalescent for change from WT to variant
censoredChangeFromWTLcfConvMeans_LabSerumVar=getCensoredMeans(WTdata,c('Laboratory','Serum','Variant'),'changeFromWTLcfConv',censoredBelow = WTdata$censoredMutant, censoredAbove= WTdata$censoredWT)
censoredChangeFromWTLcfConvMeans_SerumVar=getCensoredMeans(WTdata,c('Serum','Variant'),'changeFromWTLcfConv',censoredBelow = WTdata$censoredMutant, censoredAbove= WTdata$censoredWT)
censoredChangeFromWTLcfConvMeans_LabVar=getCensoredMeans(WTdata,c('Laboratory','Variant'),'changeFromWTLcfConv',censoredBelow = WTdata$censoredMutant, censoredAbove= WTdata$censoredWT)
censoredChangeFromWTLcfConvMeans_Var=getCensoredMeans(WTdata,c('Variant'),'changeFromWTLcfConv',censoredBelow = WTdata$censoredMutant, censoredAbove= WTdata$censoredWT)

# Now calculate mean drops compared to convalescent for change from B117 to variant
#censoredChangeFromB117LcfConvMeans_LabSerumVar=getCensoredMeans(WTdata,c('Laboratory','Serum','Variant'),'changeFromB117LcfConv',censoredBelow = WTdata$censoredMutant, censoredAbove= WTdata$censoredB117)
#censoredChangeFromB117LcfConvMeans_SerumVar=getCensoredMeans(WTdata,c('Serum','Variant'),'changeFromB117LcfConv',censoredBelow = WTdata$censoredMutant, censoredAbove= WTdata$censoredB117)
#censoredChangeFromB117LcfConvMeans_LabVar=getCensoredMeans(WTdata,c('Laboratory','Variant'),'changeFromB117LcfConv',censoredBelow = WTdata$censoredMutant, censoredAbove= WTdata$censoredB117)
#censoredChangeFromB117LcfConvMeans_Var=getCensoredMeans(WTdata,c('Variant'),'changeFromB117LcfConv',censoredBelow = WTdata$censoredMutant, censoredAbove= WTdata$censoredB117)

#And now merge all the above
WTdata=mergeMeans(WTdata,censoredWTneutLChangeFromConvMeans_Var, by=c('Variant'),includesig=F)
WTdata=mergeMeans(WTdata,censoredWTneutLChangeFromConvMeans_SerumVar,c('Serum','Variant'))
WTdata=mergeMeans(WTdata,censoredWTneutLChangeFromConvMeans_LabVar,c('Laboratory','Variant'))
WTdata=mergeMeans(WTdata,censoredWTneutLChangeFromConvMeans_LabSerumVar,c('Laboratory','Serum','Variant'))

WTdata=mergeMeans(WTdata,censoredvariantneutLChangeFromConvMeans_Var, by=c('Variant'),includesig=F)
WTdata=mergeMeans(WTdata,censoredvariantneutLChangeFromConvMeans_SerumVar,c('Serum','Variant'))
WTdata=mergeMeans(WTdata,censoredvariantneutLChangeFromConvMeans_LabVar,c('Laboratory','Variant'))
WTdata=mergeMeans(WTdata,censoredvariantneutLChangeFromConvMeans_LabSerumVar,c('Laboratory','Serum','Variant'))

WTdata=mergeMeans(WTdata,censoredChangeFromWTLcfConvMeans_Var, by=c('Variant'),includesig=F)
WTdata=mergeMeans(WTdata,censoredChangeFromWTLcfConvMeans_SerumVar,c('Serum','Variant'))
WTdata=mergeMeans(WTdata,censoredChangeFromWTLcfConvMeans_LabVar,c('Laboratory','Variant'))

#WTdata=mergeMeans(WTdata,censoredChangeFromB117LcfConvMeans_Var, by=c('Variant'),includesig=F)
#WTdata=mergeMeans(WTdata,censoredChangeFromB117LcfConvMeans_SerumVar,c('Serum','Variant'))
#WTdata=mergeMeans(WTdata,censoredChangeFromB117LcfConvMeans_LabVar,c('Laboratory','Variant'))

#offsets relative to convalescent by serum
cfConvOffsetsForVaccine=findXVals('Serum',WTdata,'changeFromWTLcfConv')
WTdata=merge(WTdata,cfConvOffsetsForVaccine, by=c('Serum','Variant','Laboratory'),all.x=T,all.y=T)

#offsets relative to convalescent by lab
cfConvOffsetsForLaboratory=findXVals('Laboratory',WTdata,'changeFromWTLcfConv',.5)
WTdata=merge(WTdata,cfConvOffsetsForLaboratory, by=c('Laboratory','Variant','Serum'),all.x=T,all.y=T)

#### Add Data relative to Pfizer means ####

# Extract only the relevant Pfizer Data
pfizerWTneutLMeans=censoredWTneutLMeans_LabSerumVar[censoredWTneutLMeans_LabSerumVar$Serum=='Pfizer',]
colnames(pfizerWTneutLMeans)[str_detect(colnames(pfizerWTneutLMeans),'Lab_Ser_Var')]=str_replace(colnames(pfizerWTneutLMeans)[str_detect(colnames(pfizerWTneutLMeans),'Lab_Ser_Var')],'Lab_Ser_Var','Pfizer_Lab_Var')
pfizerWTneutLMeans=select(pfizerWTneutLMeans,c('Laboratory','Variant','muL_WTneutL_Pfizer_Lab_Var'))

pfizerVariantneutLMeans=censoredVariantneutLMeans_LabSerumVar[censoredVariantneutLMeans_LabSerumVar$Serum=='Pfizer',]
colnames(pfizerVariantneutLMeans)[str_detect(colnames(pfizerVariantneutLMeans),'Lab_Ser_Var')]=str_replace(colnames(pfizerVariantneutLMeans)[str_detect(colnames(pfizerVariantneutLMeans),'Lab_Ser_Var')],'Lab_Ser_Var','Pfizer_Lab_Var')
pfizerVariantneutLMeans=select(pfizerVariantneutLMeans,c('Laboratory','Variant','muL_variantneutL_Pfizer_Lab_Var'))

pfizerChangeMeans=censoredChangeMeans_LabSerumVar[censoredChangeMeans_LabSerumVar$Serum=='Pfizer',]
colnames(pfizerChangeMeans)[str_detect(colnames(pfizerChangeMeans),'Lab_Ser_Var')]=str_replace(colnames(pfizerChangeMeans)[str_detect(colnames(pfizerChangeMeans),'Lab_Ser_Var')],'Lab_Ser_Var','Pfizer_Lab_Var')

# Merge pfizer data against  WT
WTdata=mergeMeans(WTdata,pfizerWTneutLMeans, by=c('Laboratory','Variant'),includesig=F)
WTdata$WTneutLChangeFromPfizer = WTdata$WTneutL-WTdata$muL_WTneutL_Pfizer_Lab_Var

#  Merge pfizer data against  variant
WTdata=mergeMeans(WTdata,pfizerVariantneutLMeans, by=c('Laboratory','Variant'),includesig=F)
WTdata$variantneutLChangeFromPfizer = WTdata$variantneutL-WTdata$muL_variantneutL_Pfizer_Lab_Var

# Merge change from WT to variant
WTdata=mergeMeans(WTdata,pfizerChangeMeans, by=c('Laboratory','Variant'),includesig=F)
WTdata$changeFromWTLcfPfizer = WTdata$changeFromWTL-WTdata$muL_changeFromWTL_Pfizer_Lab_Var

# Now calculate mean drops compared to pfizer for WT
censoredWTneutLChangeFromPfizerMeans_LabSerumVar=getCensoredMeans(WTdata,c('Laboratory','Serum','Variant'),'WTneutLChangeFromPfizer',censoredBelow = WTdata$censoredWT)
censoredWTneutLChangeFromPfizerMeans_SerumVar=getCensoredMeans(WTdata,c('Serum','Variant'),'WTneutLChangeFromPfizer',censoredBelow = WTdata$censoredWT)
censoredWTneutLChangeFromPfizerMeans_LabVar=getCensoredMeans(WTdata,c('Laboratory','Variant'),'WTneutLChangeFromPfizer',censoredBelow = WTdata$censoredWT)
censoredWTneutLChangeFromPfizerMeans_Var=getCensoredMeans(WTdata,c('Variant'),'WTneutLChangeFromPfizer',censoredBelow = WTdata$censoredWT)

# Now calculate mean drops compared to pfizer for variant
censoredvariantneutLChangeFromPfizerMeans_LabSerumVar=getCensoredMeans(WTdata,c('Laboratory','Serum','Variant'),'variantneutLChangeFromPfizer',censoredBelow = WTdata$censoredMutant)
censoredvariantneutLChangeFromPfizerMeans_SerumVar=getCensoredMeans(WTdata,c('Serum','Variant'),'variantneutLChangeFromPfizer',censoredBelow = WTdata$censoredMutant)
censoredvariantneutLChangeFromPfizerMeans_LabVar=getCensoredMeans(WTdata,c('Laboratory','Variant'),'variantneutLChangeFromPfizer',censoredBelow = WTdata$censoredMutant)
censoredvariantneutLChangeFromPfizerMeans_Var=getCensoredMeans(WTdata,c('Variant'),'variantneutLChangeFromPfizer',censoredBelow = WTdata$censoredMutant)

# Now calculate mean drops compared to pfizer for chnage from WT to variant
censoredChangeFromWTLcfPfizerMeans_LabSerumVar=getCensoredMeans(WTdata,c('Laboratory','Serum','Variant'),'changeFromWTLcfPfizer',censoredBelow = WTdata$censoredMutant, censoredAbove= WTdata$censoredWT)
censoredChangeFromWTLcfPfizerMeans_SerumVar=getCensoredMeans(WTdata,c('Serum','Variant'),'changeFromWTLcfPfizer',censoredBelow = WTdata$censoredMutant, censoredAbove= WTdata$censoredWT)
censoredChangeFromWTLcfPfizerMeans_LabVar=getCensoredMeans(WTdata,c('Laboratory','Variant'),'changeFromWTLcfPfizer',censoredBelow = WTdata$censoredMutant, censoredAbove= WTdata$censoredWT)
censoredChangeFromWTLcfPfizerMeans_Var=getCensoredMeans(WTdata,c('Variant'),'changeFromWTLcfPfizer',censoredBelow = WTdata$censoredMutant, censoredAbove= WTdata$censoredWT)

#And now merge all the above
WTdata=mergeMeans(WTdata,censoredWTneutLChangeFromPfizerMeans_Var, by=c('Variant'),includesig=F)
WTdata=mergeMeans(WTdata,censoredWTneutLChangeFromPfizerMeans_SerumVar,c('Serum','Variant'))
WTdata=mergeMeans(WTdata,censoredWTneutLChangeFromPfizerMeans_LabVar,c('Laboratory','Variant'))
WTdata=mergeMeans(WTdata,censoredWTneutLChangeFromPfizerMeans_LabSerumVar,c('Laboratory','Serum','Variant'))

WTdata=mergeMeans(WTdata,censoredvariantneutLChangeFromPfizerMeans_Var, by=c('Variant'),includesig=F)
WTdata=mergeMeans(WTdata,censoredvariantneutLChangeFromPfizerMeans_SerumVar,c('Serum','Variant'))
WTdata=mergeMeans(WTdata,censoredvariantneutLChangeFromPfizerMeans_LabVar,c('Laboratory','Variant'))
WTdata=mergeMeans(WTdata,censoredvariantneutLChangeFromPfizerMeans_LabSerumVar,c('Laboratory','Serum','Variant'))

WTdata=mergeMeans(WTdata,censoredChangeFromWTLcfPfizerMeans_Var, by=c('Variant'),includesig=F)
WTdata=mergeMeans(WTdata,censoredChangeFromWTLcfPfizerMeans_SerumVar,c('Serum','Variant'))
WTdata=mergeMeans(WTdata,censoredChangeFromWTLcfPfizerMeans_LabVar,c('Laboratory','Variant'))

#offsets relative to pfizer by serum
cfPfizerOffsetsForVaccine=findXVals('Serum',WTdata,'changeFromWTLcfPfizer')
WTdata=merge(WTdata,cfPfizerOffsetsForVaccine, by=c('Serum','Variant','Laboratory'),all.x=T,all.y=T)

#offsets relative to pfizer by lab
cfPfizerOffsetsForLaboratory=findXVals('Laboratory',WTdata,'changeFromWTLcfPfizer',.5)
WTdata=merge(WTdata,cfPfizerOffsetsForLaboratory, by=c('Laboratory','Variant','Serum'),all.x=T,all.y=T)
}