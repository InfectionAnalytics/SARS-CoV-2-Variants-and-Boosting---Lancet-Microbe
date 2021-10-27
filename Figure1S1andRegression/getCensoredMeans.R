getCensoredMeans=function(dataset=WTdata,groups=c("Laboratory","subsequentVariant"),value='changeFromWTL', censoredBelow=NULL, censoredAbove=NULL){

  ####Model Fold Change
  likelihoodfunc<-function(values,censoredBelow, censoredAbove,mu_log10,sig){
    uncensored=!(censoredBelow|censoredAbove)
    NLL=-sum(dnorm(values[uncensored],mu_log10,sig,log=T))
    NLL=NLL - sum(pnorm(values[censoredBelow],mu_log10,sig,log=T))
    NLL=NLL - sum(pnorm(values[censoredAbove],mu_log10,sig,lower.tail = F,log=T))
    NLL
  }


  ###Make Table
  if(length(groups)==1){
    censoredMeans<-data.frame(x=unique(dataset[,groups]))
    colnames(censoredMeans)=groups
  } else{
    censoredMeans=unique(dataset[,groups])
  }
  # add censored below and censored above
  # check censored has been entered
  if (is.null(censoredBelow)){
    censoredBelow = rep(F, nrow(dataset))
  } else if (length(censoredBelow)!=nrow(dataset)){
    simpleError('Error: length of censoredBelow must be the same as the number of rows in the dataset')
  }
  
  if (is.null(censoredAbove)){
    censoredAbove= rep(F, nrow(dataset))
  }  else if (length(censoredAbove)!=nrow(dataset)){
    simpleError('Error: length of censoredAbove must be the same as the number of rows in the dataset')
  }
  uncensored = !(censoredBelow|censoredAbove)
  dataset$censoredBelow = censoredBelow
  dataset$censoredAbove = censoredAbove
  dataset$uncensored = uncensored
  
  # Remove NA values
  dataset=dataset[!is.na(select(dataset,value)),]
  #Determine the number of studies in each group
  eval(expr = str2expression(paste0('censoredMeans<-dataset %>% group_by(',paste0(groups,collapse=','),') %>% summarise(nstudy=length(unique(PaperRef)))')))  
  censoredMeans<-data.frame(censoredMeans)

  
  censoredMeans$mu_guess<-NA
  censoredMeans$sig_guess<-NA
  censoredMeans$muL<-NA
  censoredMeans$sigL<-NA
  censoredMeans$exitflag<-NA
  censoredMeans$negloglikelihood<-NA
  censoredMeans$n<-NA
  censoredMeans$Q1lower=NA
  censoredMeans$Q4upper=NA

  for (i in 1:nrow(censoredMeans)){
    #Exclude WT censored data - Not anymore, now we have a censored above variable
    tempData=dataset #[dataset$censoredWT==0,]
    tempCensoredBelow=dataset$censoredBelow #[dataset$censoredWT==0]
    tempCensoredAbove=dataset$censoredAbove
    for (j in c(1:length(groups))){
      #print(paste0('groups[j]=',groups[j],' censoredMeans[i,j]=',censoredMeans[i,j]))
      useInds=tempData[,groups[j]]==censoredMeans[i,j]
      tempData=tempData[useInds,]
      tempCensoredBelow=tempCensoredBelow[which(useInds)]
      tempCensoredAbove=tempCensoredAbove[which(useInds)]
    }
    useInds=!is.na(tempData[,value])
    tempData=tempData[useInds,]
    tempCensoredBelow=tempCensoredBelow[which(useInds)]
    tempCensoredAbove=tempCensoredAbove[which(useInds)]
    tempUncensored = !(tempCensoredBelow|tempCensoredAbove)
    # Can't continue with no data
    if (nrow(tempData)!=0 & nrow(tempData[tempUncensored,])!=0 ){
      initialguess<-c(mean(tempData[tempUncensored,value]),sd(tempData[tempUncensored,value]))
      if (initialguess[2]==0 | is.na(initialguess[2])){
        censoredMeans$mu_guess[i]<-initialguess[1]
        censoredMeans$sig_guess[i]<-initialguess[2]
        censoredMeans$muL[i]<-initialguess[1]
        censoredMeans$sigL[i]<-initialguess[2]
      } else {
        tempfit<-nlm(function(x)likelihoodfunc(tempData[,value],tempCensoredBelow, tempCensoredAbove,x[1],x[2]),initialguess)
        censoredMeans$mu_guess[i]<-initialguess[1]
        censoredMeans$sig_guess[i]<-initialguess[2]
        censoredMeans$muL[i]<-tempfit$estimate[1]
        censoredMeans$sigL[i]<-tempfit$estimate[2]
        censoredMeans$exitflag[i]<-tempfit$code
        censoredMeans$negloglikelihood[i]<-tempfit$minimum
        censoredMeans$n[i]<-nrow(tempData)
      }
      censoredMeans$Q1cut[i]=quantile(tempData[,value])[4]
      censoredMeans$Q4cut[i]=quantile(tempData[,value])[2]
    }
    #print(i)
  }
  
  censoredMeans$seL=censoredMeans$sigL/sqrt(censoredMeans$n)
  censoredMeans$lower95L=censoredMeans$muL-censoredMeans$seL*1.96
  censoredMeans$upper95L=censoredMeans$muL+censoredMeans$seL*1.96
  colnames(censoredMeans)[colnames(censoredMeans)=='muL']=paste('muL',value,paste(substr(groups,1,3),collapse = '_'),sep='_')
  colnames(censoredMeans)[colnames(censoredMeans)=='sigL']=paste('sigL',value,paste(substr(groups,1,3),collapse = '_'),sep='_')
  colnames(censoredMeans)[colnames(censoredMeans)=='seL']=paste('seL',value,paste(substr(groups,1,3),collapse = '_'),sep='_')
  colnames(censoredMeans)[colnames(censoredMeans)=='Q1cut']=paste('Q1cut',value,paste(substr(groups,1,3),collapse = '_'),sep='_')
  colnames(censoredMeans)[colnames(censoredMeans)=='Q4cut']=paste('Q4cut',value,paste(substr(groups,1,3),collapse = '_'),sep='_')
  
  censoredMeans
}