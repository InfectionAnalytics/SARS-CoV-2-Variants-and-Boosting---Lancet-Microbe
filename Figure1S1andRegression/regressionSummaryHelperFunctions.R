makeResultSummaryFromRegressionList<-function(regressions, resultSummary=NULL){
  if(is.null(resultSummary)){
    xVarsList=list()
    varsList=list()
    for(k in c(1:length(regressions))){
      # Extract the xVariables from the formula
      xVarsList[[k]]=str_trim(strsplit(strsplit(regressions[[k]]$formula,'~')[[1]][2],'\\+')[[1]])
      varsList[[k]]=rownames(regressions[[k]]$coeff)
    }
    resultSummary<-makeEmptyRegressionSummaryDataFrame(data=NULL,variants=NULL, xVarsList,summaryColList=varsList)
    for(k in c(1:length(regressions))){
      if (is.null(regressions[[k]]$modelName)){
        regressions[[k]]$modelName='NoName'
      } 
      theseVariants=regressions[[k]]$mutant
      if (length(theseVariants)>1){
        theseVariants=paste0(theseVariants,collapse=',')
      }
      resultSummary<-insertRegressionResultsIntoSummaryDataFrame(resultSummary,k,regResult=regressions[[k]],variables=regressions[[k]]$xVars,modelName=regressions[[k]]$modelName, variants=theseVariants)
    }
  }
  if (length(unique(resultSummary$mutant))<=1){
    resultSummary$minAIC=as.numeric(resultSummary %>% summarise(minAIC = min(AIC, na.rm=T)))
  } else {
    resultSummary = merge(resultSummary,resultSummary %>% group_by(mutant) %>% summarise(minAIC = min(AIC, na.rm=T)),by='mutant')
  }
  resultSummary$isMinCensAIC = resultSummary$AIC==min(subset(resultSummary,isCensored)$AIC)
  resultSummary$minCensAIC=NA
  resultSummary$minCensAIC[resultSummary$isMinCensAIC==T]=min(subset(resultSummary,isCensored)$AIC)-resultSummary$minAIC
  resultSummary
  
}

makeEmptyRegressionSummaryDataFrame <- function(data=regressionData,variants=NULL, modelXVars, summaryColList=NULL){
  xVars = unique(unlist(modelXVars))
  if (!is.null(summaryColList)){
    summaryColNames = unique(unlist(summaryColList))
    summaryColNames = summaryColNames[!startsWith(summaryColNames,'intercept')]
  }
  
  nVariants = max(1,length(variants))
  resultSummary = data.frame(rownums=c(1:(nVariants*length(modelXVars))))
  resultSummary$AIC = NA
  resultSummary$neglogLik = NA
  resultSummary$df = NA
  resultSummary$isCensored = NA
  
  resultSummary$modelName = NA
  resultSummary$modelFormula = NA
  resultSummary$mutant = NA
  
  resultSummary$intercept_coeff=NA
  resultSummary$intercept_pVal=NA
  
  #If the names of the summary columns haven't been entered, get them from the dataset
  if(is.null(summaryColList)){
    for(xVarNo in c(1:length(xVars))){
      xvar = xVars[xVarNo]
      if(is.factor(data[,xvar])){
        # Add in coeff and pVal columns for each of the factor values
        for (val in as.character(unique(data[,xvar]))){
          resultSummary$newCoeff = NA
          resultSummary$newpVal = NA
          colnames(resultSummary)[colnames(resultSummary)=='newCoeff']=paste0(xvar,val,'coeff')
          colnames(resultSummary)[colnames(resultSummary)=='newpVal']=paste0(xvar,val,'pVal')
        }  
      } else {
        # Add in coeff and pVal columns for this variable
        resultSummary$newCoeff = NA
        resultSummary$newpVal = NA
        colnames(resultSummary)[colnames(resultSummary)=='newCoeff']=paste0(xvar,'coeff')
        colnames(resultSummary)[colnames(resultSummary)=='newpVal']=paste0(xvar,'pVal')
      }
    }
  } else {
    #Otherwise just add them from the summaryCols
    for (i in c(1:length(summaryColNames))){
      # Add in coeff and pVal columns for this variable
      resultSummary$newCoeff = NA
      resultSummary$newpVal = NA
      colnames(resultSummary)[colnames(resultSummary)=='newCoeff']=paste0(summaryColNames[i],'coeff')
      colnames(resultSummary)[colnames(resultSummary)=='newpVal']=paste0(summaryColNames[i],'pVal')  
    }
  }
  resultSummary
}

insertRegressionResultsIntoSummaryDataFrame<-function(resultSummary,k=1,regResult=regressions[[k]],variables=xVarsList[[j]],modelName='noName', variants='ALL'){
  if (!is.na(regResult)){
    resultSummary$AIC[k]=regResult$AIC
    resultSummary$neglogLik[k]=regResult$neglogLik
    resultSummary$df[k]=regResult$df
    resultSummary$modelName[k]=modelName
    resultSummary$modelFormula[k]=regResult$formula
    resultSummary$mutant[k]=variants
    resultSummary$isCensored[k]=regResult$isCensored
    
    # extract p values
    pInd=4
    coefInd = 1
    # This assumes we always have an intercept - ok for now
    resultSummary$intercept_pVal[k] = regResult$coeff[1,pInd]
    resultSummary$intercept_coeff[k] = regResult$coeff[1,coefInd]
    
    xNames = rownames(regResult$coeff)[!startsWith(rownames(regResult$coeff),'intercept')]
    for (xname in xNames){
      ind = which(rownames(regResult$coeff)==xname)
      if(length(ind)>0) {
        resultSummary[k,paste0(xname,'coeff')] = regResult$coeff[ind,coefInd]
        resultSummary[k,paste0(xname,'pVal')] = regResult$coeff[ind,pInd]
      }
    }
  }
  resultSummary
}
#runRegressionForAllVariantsTogether(regressions,F)
#regressions=runRegressionForAllVariantsTogether()
