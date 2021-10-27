runVariantRegression = function(dataset=WTdata,xVars=c('Serum','muL_changeFromWTL_Convalescent_Lab_Var'),yVar='changeFromWTL',offsetVars=NULL,intercept=F,mutants=NULL){
  #source(paste0(codeFolder,'NLLcensoredRegression.R'))

  #useCols=c("Laboratory","subsequentVariant","Serum",   "PaperRef","LOD",       "WTneutL",     "variantneutL","changeFromWTL","LaboratoryBAK","censoredMutant","censoredWT",      
  #        "meanlog10ConvalescentDrop","meanlog10ConvalescentDropCen","meanLog10ConvalescentWT","vacToConv_log10FoldChangeDropCen")
  
  if (is.null(mutants)){
    mutants = unique(dataset$Variant)
  }
  regressionData=dataset[!is.na(dataset$LOD)&(dataset$Variant %in% mutants)&dataset$Serum!='Covaxin',]
  regressionData$Variant=factor(regressionData$Variant)
  if (length(xVars)==0){
    modelFormula = paste(yVar,'~1')
  } else {
    i=1
    xVars=sort(xVars)
    modelFormula = paste(yVar,'~',xVars[i])
  
    while (length(xVars)>i){
      i=i+1
      modelFormula = paste(modelFormula,'+',xVars[i])
    }
    if (intercept==F){
      modelFormula = paste(modelFormula,'+ 0')
    }
  }
    j=0
    while (length(offsetVars)>j){
      j=j+1
      #Not sure why we would only do this if its int he xvars
      #if (offsetVars[j] %in% xVars){
      if(T){
        modelFormula = paste(modelFormula,'+offset(',offsetVars[j],')')
      }
    }
  
  model=lm(as.formula(modelFormula), data=regressionData)
  regData = data.frame(y=regressionData[,yVar])

  useIntercept=intercept
  isFactors=rep(F,length(xVars))
  if (length(xVars)>0){
    for( i in c(1:length(xVars))){
      xvar=xVars[i]
      if (xvar=='0'){  #NOTE: This would cause a problem if one of names of the RHS regression values was "0" and intercept=F
        useIntercept=F
      } else if (is.factor(regressionData[,xvar])){
        for (val in as.character(unique(regressionData[,xvar]))){
          regData$newVal = regressionData[,xvar]==val
          colnames(regData)[colnames(regData)=='newVal']=paste0(xvar,val)
          isFactors[i]=T
        }
      } else {
        regData$newVal=regressionData[,xvar]
        colnames(regData)[colnames(regData)=='newVal']=xvar
      }
    }
  }
  # Add intercept column if necessary
  interceptName='intercept'
  if (useIntercept){
    regData$intercept = 1
    # Remove one of the factors if there is one.
    if(sum(isFactors)>0){
      interceptName='intercept_'
      for ( i in c(1:sum(isFactors))){
        factor=xVars[isFactors][i]
        refName=paste0(factor,as.character(sort(unique(regressionData[,factor]))))[1]
        regData=select(regData, -refName)
        interceptName=paste0(interceptName,refName)
      }
      colnames(regData)[colnames(regData)=='intercept']=interceptName
    }
  }
  j=0
  #add offset variables
  if(length(xVars)>0) {
    while (length(offsetVars)>j){
      j=j+1
      if (offsetVars[j] %in% xVars){
        regData$offset =regressionData[,offsetVars[j]]
        colnames(regData)[colnames(regData)=='offset']=paste0('offset_',offsetVars[j])
      }
    }
  }

  #get the censoted values
  censLow=regressionData[,'censoredMutant']
  censHigh=regressionData[,'censoredWT']

  # Get the starting parameters
  initParams=model$coefficients
  names(initParams)[names(initParams)=='(Intercept)']=interceptName
  initParams=initParams[order(names(initParams))]
  #names(initParams)[startsWith(names(initParams),xVars[isFactors][1])]=colnames(regData)[startsWith(colnames(regData),xVars[isFactors][1])]
  initParams[length(initParams)+1]=sigma(model)
  names(initParams)[length(initParams)]='sig'
  
  # Get the censoring values
  #aUncens<-NLLcensoredRegression(initParams,regData, censLow, censHigh)
  #aCens<-NLLcensoredRegression(initParams,regData, censLow, censHigh)
  
  lowerB=rep(-5,length(initParams))
  lowerB[length(lowerB)]=1e-5
  upperB=rep(5,length(initParams))
  if(sum(is.na(initParams))>0){
    # Dont regress if NA values
    resUncens=NA
    resCens=NA
  } else {
    resUncens=runCensoredRegression(initParams,'NLLcensoredRegression',regData,lowerB,upperB)
    resCens=runCensoredRegression(initParams,'NLLcensoredRegression',regData,lowerB, upperB,censLo=censLow, censHi=censHigh)
  }
  reslm=makeNewModelStruct(model)
  result=list(linear=reslm,uncensored=resUncens, censored=resCens, formula=modelFormula, mutants=mutants)
}

NLLcensoredRegression = function(params, regData,censLow=NULL,censHigh=NULL) {
  # Values predicted by the model
  ydata=regData$y
  
  xdataInds = !(colnames(regData) %in% c('y'))&!startsWith(colnames(regData),'offset')
  xdata=regData[,xdataInds]
  if (is.data.frame(xdata)) {
    xdata=xdata[,order(colnames(regData)[xdataInds])]
  } else {
    xdata=data.frame(xdata)
    colnames(xdata)<-colnames(regData)[xdataInds]
  }
  
  offsetData= regData[,startsWith(colnames(regData),'offset')]
  
  coeffs=head(params,-1)
  coeffs=as.data.frame(matrix(rep(coeffs,nrow(xdata)),ncol=length(coeffs),byrow=T))
  estimates = rowSums(xdata*coeffs)
  #Add offset data, depending on if there is one offset or more
  noffsets=sum(startsWith(colnames(regData),'offset'))
  if (noffsets==1){
    estimates=estimates+offsetData
  } else if (noffsets>=2){
    estimates=estimates+rowSums(offsetData)
  }
  sig=tail(params,1)
  
  # determine censoring
  if (is.null(censLow)){
    censLow=rep(FALSE,length(ydata))
  }
  if (is.null(censHigh)){
    censHigh=rep(FALSE,length(ydata))
  }
  uncensored = !(censLow|censHigh)
  
  # Negative log-likelihood 
  NLLinRange=-sum(dnorm(ydata[uncensored], mean = estimates[uncensored], sd = sig,log=T), na.rm=T)
  NLLbelowLOD=-sum(pnorm(ydata[censLow], mean = estimates[censLow], sd = sig,log=T), na.rm=T)
  NLLaboveLOD=-sum(pnorm(ydata[censHigh], mean = estimates[censHigh], lower.tail = F, sd = sig,log=T), na.rm=T)
  NLL=NLLinRange+NLLbelowLOD+NLLaboveLOD
  if(is.nan(NLL)){
    NLL=length(data$y)
  }
  #print(c(params,NLL))
  NLL
}

runCensoredRegression = function(initParams,func,data,lowerB=-Inf,upperB=Inf,censLo=NULL,censHi=NULL,fitmethod="L-BFGS-B",hess=T){
  # make sure bounds are the right length
  if (length(lowerB)==1){
    lowerB=rep(lowerB,length(initParams))
  }
  if (length(upperB)==1){
    upperB=rep(upperB,length(upperB))
  }
  regressionString = paste0('regressionResult = optim(par = initParams, fn = ',func,', regData=data,censLow=censLo,censHigh=censHi, lower=lowerB, upper=upperB, method=fitmethod, hessian=hess)')
  eval(parse(text=regressionString))
  ests=head(regressionResult$par,-1)
  stderrs=head(sqrt(diag(solve(regressionResult$hessian))),-1)
  tvals=ests/stderrs
  coeffs=data.frame(Estimate=head(regressionResult$par,-1),StdErr=stderrs,tValue=tvals)
  df = nrow(data)-length(regressionResult$par)
  pvals=2*pt(-abs(tvals), df=df)
  
  regressionCoeff=data.frame(Estimate=ests,StdErr=stderrs,tValue=tvals,pValue=pvals, significance=isSignificant(pvals))
  regressionNeglogLik=regressionResult$value
  regResult=list(coeff=regressionCoeff, neglogLik=regressionNeglogLik, AIC=2*(length(regressionResult$par)+regressionNeglogLik),df=df, sigma=regressionResult$par[length(regressionResult$par)], offset=colnames(data)[startsWith(colnames(data),'offset')])
}

isSignificant = function(pvals){
  isSig=rep(' ',length(pvals))
  isSig[pvals<.1]='.'
  isSig[pvals<.05]='*'
  isSig[pvals<.01]='**'
  isSig[pvals<.001]='***'
  isSig
}

makeNewModelStruct=function(model){
  vars = as.character(attr(summary(model)$terms,'variables'))
  offsets=vars[startsWith(vars,'offset')]
  
  struct=list(coeff=as.data.frame(summary(model)$coeff),neglogLik=-logLik(model), AIC=2*(length(model$coefficients)+1-logLik(model)),df=summary(model)$df[2], sigma=sigma(model), offset=offsets)
  struct$coeff$significance=isSignificant(summary(model)$coeff[,4])
  struct
}

