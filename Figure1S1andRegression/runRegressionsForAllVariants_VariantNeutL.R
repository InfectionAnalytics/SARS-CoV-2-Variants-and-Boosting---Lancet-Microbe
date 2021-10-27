runRegressionForAllVariants_VariantNeutL<-function(data, regressions=NULL, doRegressions=T, useModels=NULL,yVar='variantneutL', usePseudovirus=T, useConvalescent=F, xVarsList = NULL, baseVars=NULL, offsetVars=NULL, useIntercept=T){
variants = c('B.1.1.7','B.1.351','P.1','B.1.617.2','B.1.617.1')
#yVar='changeFromWTL'
#yVar='variantneutL'

if (is.null(baseVars)){
  baseVars = c('WTneutL','muL_changeFromWTL_Var')#,'WTLcat')
  baseVars = c('muL_changeFromWTL_Var')#,'WTLcat')
  #baseVars=c()
}
baseVarArray = c('WTneutL','muL_changeFromWTL_Var')#,'WTLcat')
baseVarNameArray = c('WTneut','avgVarDrop')#,'WTLcat')
for (i in c(1:length(baseVars))){
  baseVarNames = c()
  if (length(baseVars)>0){
    for(j in c(1:length(baseVars))){
      basevar=baseVars[j]
      baseVarNames[j]=baseVarNameArray[baseVarArray==basevar]
    }
  }
}
if(is.null(xVarsList)){
  xVarsList = list(c(),c('Serum'),c('Laboratory'),c('muL_WTneutL_Convalescent_Lab'),
                 c('Serum','muL_WTneutL_Convalescent_Lab'),c('Laboratory','muL_WTneutL_Convalescent_Lab'),
                 c('Serum','Laboratory'),
                 c('Serum','Laboratory','muL_WTneutL_Convalescent_Lab'))
  xVarsList = list(c(),c('Serum'),c('Laboratory'),c('WTLcat'),
                 c('Serum','WTLcat'),c('Laboratory','WTLcat'),
                 c('Serum','Laboratory'),
                 c('Serum','Laboratory','WTLcat'))
  xVarsList = list(c(),c('muL_WTneutL_Convalescent_Lab'),
                   c('Serum'),
                   c('Serum','muL_WTneutL_Convalescent_Lab'))
}


modelNames=c()
variables=c('Variant','Serum','Laboratory','WTLcat','WTneutL','muL_changeFromWTL_Var', 'Group','muL_WTneutL_Convalescent_Lab', 'muL_changeFromWTL_Convalescent_Lab_Var')
variableNames=c('Variant','Serum','Lab','WTCat','WTneut','avgVarDrop', 'Group','avgConvWTneut','avgConvChangeFromWT')
if (is.null(offsetVars)){
  offsetVars = c('WTneutL')
}
for (i in c(1:length(xVarsList))){
  istr=ifelse(i<10,paste0('0',as.character(i)),as.character(i))
  xVars = xVarsList[[i]]
  xvarnames = c()
   if (length(xVars)>0){
    for(j in c(1:length(xVars))){
      xvar=xVars[j]
      xvarnames[j]=variableNames[variables==xvar]
    }
  }
  xvarnames=c(baseVarNames,xvarnames)
  name=paste(xvarnames,collapse='+')
  modelNames[i]=paste(istr,name)
}

# Check which models we want to use - this isn't currently used for anythign except including the base variables (that are in all models)
# but we can change teh models that we try by changinf useModels
if (is.null(useModels)){
  useModels = c(1:length(modelNames))
  xVarsListNew=list()
  i=1
  for (modelNo in useModels){
    if (!is.na(baseVars)){
      xVarsListNew[[i]]=c(baseVars,xVarsList[[modelNo]])
    } else {
      xVarsListNew[[i]]=c(xVarsList[[modelNo]])
    }
    i=i+1
  }
  xVarsList=xVarsListNew
  modelNames=modelNames[useModels]
}

baseModelsList=list()
baseModelsNames=c()
# Add in the three basic models
# if(useIntercept){
#   baseModelsList=list(c())
#   baseModelsNames=c('000 Intercept')
# } 

for(i in c(1:length(baseVars))){
  if(!is.na(baseVars[i])){
    baseModelsList=c(baseModelsList,c(baseVars[i]))
    baseModelsNames=c(baseModelsNames,paste0('00',i,baseVarNames[i]))
  }
}

xVarsList=c(xVarsList,baseModelsList)
modelNames=c(modelNames,baseModelsNames)

#modelNames = c('01 Serum','02 Lab','Serum+Lab',
#               '06 muConvDrop','Serum+muConvDrop','Lab+muConvDrop',
#               '03 Variant',
#               '07 Serum+Variant','08 Lab+Variant','0x9 Serum+Lab+Variant','10 Serum+muConvDrop+Variant','11 Variant+muConvDrop')
k=1

if (is.null(regressions)){
  doRegressions=T
}

regressionData = getDataForRegression(data,useConvalescent, writeData=F)

if (!usePseudovirus) {
  regressionData=subset(regressionData, Isolate != 'Pseudovirus')
  regressionData=rename(regressionData,muL_changeFromWTLAllIsolates_Var=muL_changeFromWTL_Var,muL_changeFromWTL_Var=muL_changeFromWTLNoPseudo_Var, muL_WTneutL_ConvalescentAllIsolates_Lab=muL_WTneutL_Convalescent_Lab, muL_WTneutL_Convalescent_Lab=muL_WTneutL_ConvalescentNoPseudo_Lab)
}
requiredCols=unique(c(unlist(xVarsList),offsetVars))
#Only use data where all variables we might use are present
regressionData<-regressionData[rowSums(is.na(regressionData[,requiredCols]))==0,]
regressionData$Variant = factor(regressionData$Variant)
resultSummary=makeEmptyRegressionSummaryDataFrame(regressionData, NULL,xVarsList)

if (doRegressions){
  if (exists('regressions')){regressionsBAK =regressions}
  regressions=list()

  for(j in c(1:length(xVarsList))){
      # Check if we want to include convalescent
      # if (sum(str_count(xVarsList[[j]],'Convalescent'))>=1){
      #   regressionData = WTdata[WTdata$Serum!='Convalescent',]
      # }
    
      reg=runVariantRegression(dataset=regressionData, xVars=xVarsList[[j]],yVar=yVar, offsetVars = offsetVars, intercept=useIntercept, mutants=variants)
      
      regressions[[k]]=reg$censored
      if (is.na(reg$censored)){
        regressions[[k]]=reg$linear
        regressions[[k]]$isCensored=F
      } else{
        regressions[[k]]$isCensored=T
      }
      regressions[[k]]$formula = reg$formula
      regressions[[k]]$mutant = reg$mutants
      regressions[[k]]$xVars = xVarsList[[j]]
      regressions[[k]]$variants = variants
      regressions[[k]]$modelName = modelNames[j]
      
      k=k+1
      #saveRDS(regressions,file=paste0(outputFolder,'regressionsSave.RDS'))
  }
} else {
  regressionsNew = list()
  i=1
  for (k in c(1:length(regressions))){
    if (sum(useModels==k)>0){
      regressionsNew[i] =regressions[k]
      i=i+1
    }
  }
  regressions=regressionsNew
}
resultSummary = makeResultSummaryFromRegressionList(regressions)
if (usePseudovirus){pseudoString='WithPseudo'} else {pseudoString='NoPseudo'}
if (useConvalescent){convString='WithConv'} else {convString='NoConv'}
infoString = paste0(pseudoString,convString)
#write.csv(resultSummary,paste0(outputFolder,yVar,'allVariantRegressionSummary_',infoString,'_pVals',datestr,".csv"))

makeAICpValPlots(resultSummary, infoString)
regressions
}

makeAICpValPlots<-function(resultSummary, infoString=''){
  # Check if we have entered a dataframe (resultSummary) or a list (regressions)
  if (class(resultSummary)=='list'){
    resultSummary = makeResultSummaryFromRegressionList(resultSummary)
  }
  AICplot=makeAICplot(resultSummary)
  pValPlot=makePvalPlot(resultSummary)
  
  plotWidth=10
  plotHeight=4
  print('Made Plots')
  #ggsave(paste0(plotFolder,yVar,'allVariantModels_',infoString,'_AICs',datestr,".png"),AICplot,width=plotWidth, height=plotHeight)
  #ggsave(paste0(plotFolder,yVar,'allVariantModels_',infoString,'_pVals',datestr,".png"),pValPlot,width=plotWidth, height=plotHeight)
  #write.csv(resultSummary,paste0(outputFolder,yVar,'allVariantRegressionSummary',pseudoString,'_pVals',datestr,".csv"))
}

makeAICplot<-function(resultSummary){
  alphaRange=1*resultSummary$isCensored
  alphaRange[alphaRange==0]=.4
  alphaRange=sort(unique(alphaRange))
  if (length(alphaRange)==1){
    alphaRange=c(alphaRange,alphaRange)
  }
  AICplot=(ggplot(resultSummary,aes(x=modelName,y=AIC-minAIC,fill=modelName, alpha=isCensored))
           +geom_col()
           +geom_text(aes(y=5+(AIC-minAIC),label=round(AIC,2)))
           +scale_alpha_discrete(range=alphaRange)
           +geom_point(aes(y=0*minCensAIC), show.legend = F)
           +geom_point(aes(y=minCensAIC), show.legend = F)
           #+facet_wrap(~mapVariant(mutant),nrow=1, scales='free_y')
           +theme_anglexaxis
  )
}

makePvalPlot <- function(resultSummary){
  resultSummaryPvals<-makepValsResultSummary(resultSummary)
  pValPlot=(ggplot(resultSummaryPvals,aes(x=modelName,y=coef,alpha=pValGp,fill=modelName))
   +geom_tile(width=1,height=1)
   #+facet_wrap(~mapVariant(mutant),nrow=1)
   +theme_anglexaxis
   +scale_alpha_manual(values=c(1,.6,.25,.1,0),breaks=levels(resultSummaryPvals$pValGp) )#scale_alpha_binned(breaks=c(-10,-5,-3,-2.3,-2.1,0),range=c(1,0))
   +theme(legend.text=element_text(size=5))
  )
  xmin=.5
  xmax=length(unique(resultSummaryPvals$modelName))+.5
  pValPlot = (pValPlot
              +geom_hline(mapping=aes(xmin=xmin,xmax=xmax,yintercept=which(startsWith(as.character(unique(resultSummaryPvals$coef)),'Serum'))[1]-.5),color='black')
              +geom_hline(mapping=aes(xmin=xmin,xmax=xmax,yintercept=which(startsWith(as.character(unique(resultSummaryPvals$coef)),'Laboratory'))[1]-.5),color='black')
              +geom_hline(mapping=aes(xmin=xmin,xmax=xmax,yintercept=which(startsWith(as.character(unique(resultSummaryPvals$coef)),'muL'))[1]-.5),color='black')
              +geom_hline(mapping=aes(xmin=xmin,xmax=xmax,yintercept=which(startsWith(as.character(unique(resultSummaryPvals$coef)),'Variant'))[1]-.5),color='black')
              +geom_hline(mapping=aes(xmin=xmin,xmax=xmax,yintercept=which(startsWith(as.character(unique(resultSummaryPvals$coef)),'WT'))[1]-.5),color='black')
              +geom_hline(mapping=aes(xmin=xmin,xmax=xmax,yintercept=which(startsWith(as.character(unique(resultSummaryPvals$coef)),'WTLcat'))[1]-.5),color='black')
              #+geom_rect(mapping=aes(xmin=xmin,xmax=xmax,ymin=17.55,ymax=18.45),fill='grey',alpha=.5)
              +geom_tile(width=1,height=.9)
  )
}
  
makepValsResultSummary <- function(resultSummary){
  resultSummaryPvals = melt(resultSummary,measure.vars = which(endsWith(colnames(resultSummary),'pVal')))
  resultSummaryPvals = rename(resultSummaryPvals, pVal = value, coef=variable)
  resultSummaryPvals[is.na(resultSummaryPvals)]=1
  # Now we need to set the pValue Group
  pValBreaks=c(0,.0001,.001,.01,.05,1)
  pValGpLabels=paste0('<',pValBreaks[2:length(pValBreaks)])
  pValGpLabels[length(pValGpLabels)]='NS'
  resultSummaryPvals$pValGp = cut(resultSummaryPvals$pVal,pValBreaks,labels=pValGpLabels)
  #resultSummaryPvals$pValGp=factor([resultSummaryPvals$logpVal< -10]=-10
  #resultSummaryPvals$logpVal[is.na(resultSummaryPvals$logpVal)]=0
  resultSummaryPvals$wasUsed=is.na(resultSummaryPvals$pVal)
  resultSummaryPvals
}
compareModels<-function(resultSummary,m1,m2,printResults=T){
  if (class(resultSummary)=='list'){
    resultSummary=makeResultSummaryFromRegressionList(resultSummary)
  }
  m=c(m1,m2);
  chisqVal=2*(resultSummary$neglogLik[m[1]]-resultSummary$neglogLik[m[2]]);
  dfVal=resultSummary$df[m[1]]-resultSummary$df[m[2]];
  pVal=pchisq(chisqVal, df=dfVal, lower.tail=FALSE)
  if(printResults){
    print(paste('Comparing Models',resultSummary$modelName[m[1]],'and',resultSummary$modelName[m[2]]))
    print(resultSummary$modelFormula[m[1]])
    print(resultSummary$modelFormula[m[2]])
    print(paste('Chi-Sq Val =',chisqVal,'df =',dfVal))
    print(paste('pVal =',pVal))
    print(paste('Rsq extra = ',round((1-resultSummary$neglogLik[m[2]]/resultSummary$neglogLik[m[1]])*100,2),'%'))
  }
  pVal
}

getDataForRegression<-function(data, keepConvalescent=F, writeData=T){
  # Select the regression data that has all the required variables
  #regressionData = WTdata
  #regressionData = WTdata[WTdata$Serum!='Convalescent',]
  if (!keepConvalescent){
    regressionData = subset(data, Serum!='Convalescent')
    convString = 'NoConv'
  } else {
    regressionData=data
    convString='withConv'
  }
  
  # Get mean WT convalescent by lab
  censoredWTneutLConvMeans_Lab=getCensoredMeans(subset(data,Serum=='Convalescent'),c('Laboratory'),'WTneutL',censoredBelow = subset(data,Serum=='Convalescent')$censoredWT)
  censoredWTneutLConvMeans_Lab<-select(censoredWTneutLConvMeans_Lab,Laboratory=Laboratory,muL_WTneutL_Convalescent_Lab=muL_WTneutL_Lab)
  
  # Get Q1 and Q4 for WT by lab and serum
  censoredWTneutLMeans_Lab_Ser=getCensoredMeans(data,c('Laboratory','Serum'),'WTneutL',censoredBelow = data$censoredWT)
  censoredWTneutLMeans_Lab_Ser<-select(censoredWTneutLMeans_Lab_Ser,Laboratory,Serum,Q1cut_WTneutL_Lab_Ser,Q4cut_WTneutL_Lab_Ser)
  
  censoredWTneutLConvMeansNoPseudo_Lab=getCensoredMeans(subset(data,Serum=='Convalescent'& Isolate!='Pseudovirus'),c('Laboratory'),'WTneutL',censoredBelow = subset(data,Serum=='Convalescent'& Isolate!='Pseudovirus')$censoredWT)
  censoredWTneutLConvMeansNoPseudo_Lab<-select(censoredWTneutLConvMeansNoPseudo_Lab,Laboratory=Laboratory,muL_WTneutL_ConvalescentNoPseudo_Lab=muL_WTneutL_Lab)
  
  # Get mean drop per variant
  censoredChangeMeansNoPseudo_Var=getCensoredMeans(subset(data,Isolate!='Pseudovirus'),c('Variant'),'changeFromWTL',censoredBelow = subset(data,Isolate!='Pseudovirus')$censoredMutant, censoredAbove= subset(data,Isolate!='Pseudovirus')$censoredWT)
  censoredChangeMeansNoPseudo_Var<-select(censoredChangeMeansNoPseudo_Var,Variant=Variant,muL_changeFromWTLNoPseudo_Var=muL_changeFromWTL_Var)
  
  # Merge lab means with data
  regressionData<-merge(regressionData,censoredWTneutLConvMeans_Lab, by='Laboratory',all.x=F, all.y=F)
  regressionData<-merge(regressionData,censoredWTneutLMeans_Lab_Ser, by=c('Laboratory','Serum'),all.x=T, all.y=F)
  regressionData<-merge(regressionData,censoredWTneutLConvMeansNoPseudo_Lab, by='Laboratory',all.x=T,all.y=F)
  regressionData<-merge(regressionData,censoredChangeMeansNoPseudo_Var, by='Variant',all.x=T, all.y=F)
  regressionData$WTLcat='mid'
  regressionData$WTLcat[regressionData$WTneutL>=regressionData$Q1cut_WTneutL_Lab_Ser]='high'
  regressionData$WTLcat[regressionData$WTneutL<=regressionData$Q4cut_WTneutL_Lab_Ser]='low'
  regressionData$WTLcat=factor(regressionData$WTLcat,levels=c('mid','high','low'))
  requiredVariables = c('Variant','Laboratory','Serum','Group','Isolate','LOD','censoredWT','censoredMutant','WTneutL','variantneutL','changeFromWTL','muL_WTneutL_Convalescent_Lab','Q1cut_WTneutL_Lab_Ser','Q4cut_WTneutL_Lab_Ser','muL_WTneutL_ConvalescentNoPseudo_Lab','muL_changeFromWTL_Var','muL_changeFromWTLNoPseudo_Var','muL_changeFromWTL_Lab_Var','WTneutLChangeFromConv','PaperRef','variantneutLChangeFromConv','muL_WTneutL_Convalescent_Lab_Var','muL_changeFromWTL_Convalescent_Lab_Var','WTLcat')
  requiredVariables = requiredVariables[requiredVariables %in% colnames(regressionData)]
  regressionData = select(regressionData,requiredVariables)
  if (writeData) {
    write.csv(regressionData, paste0(outputFolder,'newRegressionData_',convString,'.csv'),row.names = F)
  }
  regressionData
}

makeNewAICplot<-function(regsummary,str=''){
  if(class(regsummary)!='data.frame'){
    regsummary=regsummary(variantNeutTitreDropWTOffsetRegNoConv)
  }
  regsummary=regsummary %>% mutate(modelName=str_replace_all(modelName,'NA','offsets'))
  AICplot=makeAICplot(regsummary)
  aicdiffs=regsummary$AIC-regsummary$minAIC
  
  for (i in seq(1,11,2)){
    AICplot=AICplot+geom_bracket(data=regsummary,xmin=i,xmax=i+1,y.position = aicdiffs[i]+.1*max(aicdiffs),label=round(compareModels(regsummary,i,i+1,F),3),type='text')
  }
  offsetstr=str_replace_all(str_split(regsummary$modelFormula[1],'offset\\(')[[1]],' \\)',' ')
  offsetstr=trimws(paste0(offsetstr[2:length(offsetstr)],collapse=''))
  AICplot=AICplot+labs(title=paste(str_split(regsummary$modelFormula[1],' ')[[1]][1],str,'Offsets:',offsetstr))+ylim(c(0,1.2*max(aicdiffs)))
}
