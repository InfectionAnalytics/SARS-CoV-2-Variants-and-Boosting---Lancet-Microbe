library(stringr)
library(mvtnorm)
library(lemon)
library(ggplot2)
library(car)


#### Set directory to same directory as the r-script
currentdirectory<-dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(currentdirectory)

########################## Import Tables For fitting to severe cases only ####################################
SummaryTable_Efficacy_NeutRatio_SD_SEM_Severe<-read.csv(paste0(currentdirectory,"/SummaryTable_Efficacy_NeutRatio_SD_SEM_Severe.csv"), fileEncoding="UTF-8-BOM") 
censoredChangeMeans_Var<-read.csv(paste0(currentdirectory,"/censoredChangeMeans_Var_20210812.csv"),fileEncoding="UTF-8-BOM")

### The Neutralisation ratio of vaccine to convalescence using reported neut titres.
SummaryTable_Efficacy_NeutRatio_SD_SEM_Severe$NeutRatio_Reported=log10(SummaryTable_Efficacy_NeutRatio_SD_SEM_Severe$NeutMean/SummaryTable_Efficacy_NeutRatio_SD_SEM_Severe$NeutConv)


### Import Tables For fitting Models (mild infection cases)
SummaryTable_Efficacy_NeutRatio_SD_SEM<-read.csv(paste0(currentdirectory,"/SummaryTable_Efficacy_NeutRatio_SD_SEM.csv"), fileEncoding="UTF-8-BOM") 
Variant_SevereEfficacy<-read.csv(paste0(currentdirectory,"/Variant_Efficacy_Severe_20210721.csv"), fileEncoding="UTF-8-BOM") 
colnames(Variant_SevereEfficacy)[1]<-"Study"

### The Neutralisation ratio of vaccine to convalescence using reported neut titres.
SummaryTable_Efficacy_NeutRatio_SD_SEM$NeutRatio_Reported=log10(SummaryTable_Efficacy_NeutRatio_SD_SEM$NeutMean/SummaryTable_Efficacy_NeutRatio_SD_SEM$NeutConv)

SummaryTable_Efficacy_NeutRatio_SD_SEM_Severe$TechnicalName<-SummaryTable_Efficacy_NeutRatio_SD_SEM$TechnicalName[match(SummaryTable_Efficacy_NeutRatio_SD_SEM_Severe$Study,SummaryTable_Efficacy_NeutRatio_SD_SEM$Study)]


SummaryTable_Efficacy_NeutRatio_SD_SEM_Severe$TechnicalName[SummaryTable_Efficacy_NeutRatio_SD_SEM_Severe$TechnicalName==" rAd26-S+rAd5-S"]="Gam-COVID-Vac"
SummaryTable_Efficacy_NeutRatio_SD_SEM$TechnicalName[SummaryTable_Efficacy_NeutRatio_SD_SEM$TechnicalName==" rAd26-S+rAd5-S"]="Gam-COVID-Vac"

###Boosting data
BoostingTable<-read.csv(paste0(currentdirectory,"/BoostingTableSummary_20210726.csv"),fileEncoding="UTF-8-BOM")

ModernaLevelPreBoost=SummaryTable_Efficacy_NeutRatio_SD_SEM$NeutRatio_Reported[SummaryTable_Efficacy_NeutRatio_SD_SEM$Study=="Moderna"]
ModernaLevelPreBoost_Pf=SummaryTable_Efficacy_NeutRatio_SD_SEM$NeutRatio_Reported[SummaryTable_Efficacy_NeutRatio_SD_SEM$Study=="Pfizer"]
AverageLevel=0.5*(ModernaLevelPreBoost_Pf+ModernaLevelPreBoost)
BoostingTable$ReferenceNormalisation<-SummaryTable_Efficacy_NeutRatio_SD_SEM$NeutRatio_Reported[match(BoostingTable$Reference,SummaryTable_Efficacy_NeutRatio_SD_SEM$Study)]
BoostingTable$ReferenceNormalisation[BoostingTable$Reference=="Pfizer/Moderna"]=AverageLevel
BoostingTable$EstimatedLevelAfter<-log10((10^BoostingTable$ReferenceNormalisation)*BoostingTable$FoldChange)

BoostingAverage=mean(BoostingTable$EstimatedLevelAfter[BoostingTable$PreviousImmune=="Convalescent"])
BoostingRatio_SE=sd(BoostingTable$EstimatedLevelAfter[BoostingTable$PreviousImmune=="Convalescent"])/sqrt(length(BoostingTable$EstimatedLevelAfter[BoostingTable$PreviousImmune=="Convalescent"]))

####Variant Key
VariantLabel<-data.frame("Variant"=c("Ancestral","B.1.1.7","B.1.617.2","B.1.351"),
                         "Label"=c("Ancestral","Alpha (B.1.1.7)","Delta (B.1.617.2)","Beta (B.1.351)"))



##########################################################################################
##################### The Models

####Logistic Model
ProbRemainUninfected=function(logTitre,logk,C50){1/(1+exp(-exp(logk)*(logTitre-C50)))}

LogisticModel_PercentUninfected=function(mu_titre,sig_titre,logk,C50){
  NumInteration<-max(length(mu_titre),length(sig_titre),length(logk),length(C50))
  Output<-NULL
  
  if (length(C50)==1) {
    C50=rep(C50,NumInteration)
  }
  
  if (length(logk)==1) {
    logk=rep(logk,NumInteration)
  }
  
  if (length(sig_titre)==1) {
    sig_titre=rep(sig_titre,NumInteration)
  }
  
  if (length(mu_titre)==1) {
    mu_titre=rep(mu_titre,NumInteration)
  }
  
  for (i in 1:NumInteration) {
    Step=sig_titre[i]*0.001
    IntegralVector=seq(mu_titre[i]-5*sig_titre[i],mu_titre[i]+5*sig_titre[i],by=Step)
    Output[i]=sum(ProbRemainUninfected(IntegralVector,logk[i],C50[i])*dnorm(IntegralVector,mu_titre[i],sig_titre[i]))*Step
  }
  Output
}


####################################################################################################################################
#######  Logistic model for Raw Efficacy Counts (Different n50 for both severe and mild infection; shared slope) ###################

FittingLogistic_Raw_Combined_Diff_EC50_Only<-function(logRisk0_Severe,logRisk0,logk,C50_Severe,C50,N_C_Severe,N_V_Severe,Inf_C_Severe,Inf_V_Severe,MeanVector_Severe,SDVector_Severe,N_C,N_V,Inf_C,Inf_V,MeanVector,SDVector)
{
  
  Risk0=exp(logRisk0_Severe)
  
  if (length(SDVector_Severe)==1) {
    SDVector_Severe=rep(SDVector_Severe,length(Efficacy))
  }
  
  IndexNA=(is.na(N_C_Severe) | is.na(MeanVector_Severe) | is.na(SDVector_Severe))
  N_C_SEVERE=N_C_Severe[!IndexNA]
  N_V_SEVERE=N_V_Severe[!IndexNA]
  Inf_V_Severe=Inf_V_Severe[!IndexNA]
  Inf_C_Severe=Inf_C_Severe[!IndexNA]
  MeanVector_Severe=MeanVector_Severe[!IndexNA]
  SDVector_Severe=SDVector_Severe[!IndexNA]
  
  
  LL_Severe=0
  for (i in 1:length(N_C_SEVERE)) {
    
    LL_Severe=LL_Severe-log(dbinom(Inf_C_Severe[i],N_C_Severe[i],Risk0[i]))-log(dbinom(Inf_V_Severe[i],N_V_Severe[i],Risk0[i]*(1-LogisticModel_PercentUninfected(MeanVector_Severe[i],SDVector_Severe[i],logk,C50_Severe))))
  }
  
  
  Risk0=exp(logRisk0)
  
  if (length(SDVector)==1) {
    SDVector=rep(SDVector,length(Efficacy))
  }
  
  IndexNA=(is.na(N_C) | is.na(MeanVector) | is.na(SDVector))
  N_C=N_C[!IndexNA]
  N_V=N_V[!IndexNA]
  Inf_V=Inf_V[!IndexNA]
  Inf_C=Inf_C[!IndexNA]
  MeanVector=MeanVector[!IndexNA]
  SDVector=SDVector[!IndexNA]
  
  
  LL=0
  for (i in 1:length(N_C)) {
    
    LL=LL-log(dbinom(Inf_C[i],N_C[i],Risk0[i]))-log(dbinom(Inf_V[i],N_V[i],Risk0[i]*(1-LogisticModel_PercentUninfected(MeanVector[i],SDVector[i],logk,C50))))
  }
  
  
  LL_TOT = LL+LL_Severe
  
}



#########Fit model for severe (exactly repeating the Nature Med paper - not fit to new data just generating the same fit again in order to have hessian/Cov matrix for CI calculations)
#Initial values (just pick random initial values)
LogisticEstimate=c("C50_Severe"=log10(runif(1,0,0.001)),"logk"=log(runif(1,0,1)),"C50"=log10(runif(1,0,0.5)) )

ModelWithData=function(p) {
Output<-FittingLogistic_Raw_Combined_Diff_EC50_Only(p[1:nrow(SummaryTable_Efficacy_NeutRatio_SD_SEM_Severe)],
                                            p[1:nrow(SummaryTable_Efficacy_NeutRatio_SD_SEM)+nrow(SummaryTable_Efficacy_NeutRatio_SD_SEM_Severe)],
                                            p[nrow(SummaryTable_Efficacy_NeutRatio_SD_SEM)+nrow(SummaryTable_Efficacy_NeutRatio_SD_SEM_Severe)+1],
                                            p[nrow(SummaryTable_Efficacy_NeutRatio_SD_SEM)+nrow(SummaryTable_Efficacy_NeutRatio_SD_SEM_Severe)+2],
                                            p[nrow(SummaryTable_Efficacy_NeutRatio_SD_SEM)+nrow(SummaryTable_Efficacy_NeutRatio_SD_SEM_Severe)+3],
                                            SummaryTable_Efficacy_NeutRatio_SD_SEM_Severe$NumCont[!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM_Severe$NumCont)],
                                            SummaryTable_Efficacy_NeutRatio_SD_SEM_Severe$NumVac[!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM_Severe$NumCont)],
                                            SummaryTable_Efficacy_NeutRatio_SD_SEM_Severe$InfCont[!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM_Severe$NumCont)],
                                            SummaryTable_Efficacy_NeutRatio_SD_SEM_Severe$InfVac[!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM_Severe$NumCont)],
                                            SummaryTable_Efficacy_NeutRatio_SD_SEM_Severe$NeutRatio_Reported[!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM_Severe$NumCont)],
                                            SummaryTable_Efficacy_NeutRatio_SD_SEM_Severe$PooledSD[!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM_Severe$NumCont)], 
                                            SummaryTable_Efficacy_NeutRatio_SD_SEM$NumCont[!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM$NumCont)],
                                            SummaryTable_Efficacy_NeutRatio_SD_SEM$NumVac[!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM$NumCont)],
                                            SummaryTable_Efficacy_NeutRatio_SD_SEM$InfCont[!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM$NumCont)],
                                            SummaryTable_Efficacy_NeutRatio_SD_SEM$InfVac[!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM$NumCont)],
                                            SummaryTable_Efficacy_NeutRatio_SD_SEM$NeutRatio_Reported[!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM$NumCont)],
                                            SummaryTable_Efficacy_NeutRatio_SD_SEM$PooledSD[!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM$NumCont)])
Output
  }


#Minimize the negative of the log-likelihood value to fit both severe and mild cases with different n50 for each; but with a shared slope. 
#Will use reported mean of titre from the original studies.
FittedLogistic_RawEfficacy_MeanRept_SDPool_Diff_EC50_Only<-nlm(function(p){
  FittingLogistic_Raw_Combined_Diff_EC50_Only(p[1:nrow(SummaryTable_Efficacy_NeutRatio_SD_SEM_Severe)],
                                              p[1:nrow(SummaryTable_Efficacy_NeutRatio_SD_SEM)+nrow(SummaryTable_Efficacy_NeutRatio_SD_SEM_Severe)],
                                              p[nrow(SummaryTable_Efficacy_NeutRatio_SD_SEM)+nrow(SummaryTable_Efficacy_NeutRatio_SD_SEM_Severe)+1],
                                              p[nrow(SummaryTable_Efficacy_NeutRatio_SD_SEM)+nrow(SummaryTable_Efficacy_NeutRatio_SD_SEM_Severe)+2],
                                              p[nrow(SummaryTable_Efficacy_NeutRatio_SD_SEM)+nrow(SummaryTable_Efficacy_NeutRatio_SD_SEM_Severe)+3],
                                              SummaryTable_Efficacy_NeutRatio_SD_SEM_Severe$NumCont[!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM_Severe$NumCont)],
                                              SummaryTable_Efficacy_NeutRatio_SD_SEM_Severe$NumVac[!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM_Severe$NumCont)],
                                              SummaryTable_Efficacy_NeutRatio_SD_SEM_Severe$InfCont[!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM_Severe$NumCont)],
                                              SummaryTable_Efficacy_NeutRatio_SD_SEM_Severe$InfVac[!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM_Severe$NumCont)],
                                              SummaryTable_Efficacy_NeutRatio_SD_SEM_Severe$NeutRatio_Reported[!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM_Severe$NumCont)],
                                              SummaryTable_Efficacy_NeutRatio_SD_SEM_Severe$PooledSD[!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM_Severe$NumCont)], 
                                              SummaryTable_Efficacy_NeutRatio_SD_SEM$NumCont[!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM$NumCont)],
                                              SummaryTable_Efficacy_NeutRatio_SD_SEM$NumVac[!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM$NumCont)],
                                              SummaryTable_Efficacy_NeutRatio_SD_SEM$InfCont[!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM$NumCont)],
                                              SummaryTable_Efficacy_NeutRatio_SD_SEM$InfVac[!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM$NumCont)],
                                              SummaryTable_Efficacy_NeutRatio_SD_SEM$NeutRatio_Reported[!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM$NumCont)],
                                              SummaryTable_Efficacy_NeutRatio_SD_SEM$PooledSD[!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM$NumCont)])},
  c(log(SummaryTable_Efficacy_NeutRatio_SD_SEM_Severe$InfCont[!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM_Severe$NumCont)]/SummaryTable_Efficacy_NeutRatio_SD_SEM_Severe$NumCont[!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM_Severe$NumCont)]),log(SummaryTable_Efficacy_NeutRatio_SD_SEM$InfCont[!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM$NumCont)]/SummaryTable_Efficacy_NeutRatio_SD_SEM$NumCont[!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM$NumCont)]),LogisticEstimate)
  ,hessian=TRUE,iterlim=1000)

indexOfSevereParameters<-nrow(FittedLogistic_RawEfficacy_MeanRept_SDPool_Diff_EC50_Only$hessian)+c(-2,-1)
CovS<-solve(FittedLogistic_RawEfficacy_MeanRept_SDPool_Diff_EC50_Only$hessian)[indexOfSevereParameters,indexOfSevereParameters]


###Bootstrap including for predictive model
N=10000
ModelParamtemp=rmvnorm(N,mean=c(FittedLogistic_RawEfficacy_MeanRept_SDPool_Diff_EC50_Only$estimate[indexOfSevereParameters]),sigma=CovS)
SDrandom=rnorm(N,mean=SummaryTable_Efficacy_NeutRatio_SD_SEM$PooledSD[1],sd=SummaryTable_Efficacy_NeutRatio_SD_SEM$SE_PooledSD[1])
standarderrorfromdata=max(SummaryTable_Efficacy_NeutRatio_SD_SEM$SEM)

VariantList<-c("Ancestral","B.1.1.7","B.1.351","B.1.617.2")

NeutValuelog10=seq(log10(0.03),log10(80),by=0.1)
UpperIntervalsVariant_Predict_Severe<-data.frame("Variant"=rep(VariantList,each=length(NeutValuelog10)),"NeutRatio_Reported"=rep(NeutValuelog10,length(VariantList)),"Lower"=NA,"Efficacy"=NA,"Upper"=NA)
UpperIntervalsVariant_Predict_Severe$Study="Upper95Predict"
UpperIntervalsVariant_Predict_Severe$Method="RCT"
UpperIntervalsVariant_Predict_Severe$VariantLabel<-VariantLabel$Label[match(UpperIntervalsVariant_Predict_Severe$Variant,VariantLabel$Variant)]

LowerBound=0.025
UpperBound=0.975

for (i in 1:nrow(UpperIntervalsVariant_Predict_Severe)) {
  if (UpperIntervalsVariant_Predict_Severe$Variant[i]=="Ancestral") {
    MeanRandom=rnorm(N,mean=UpperIntervalsVariant_Predict_Severe$NeutRatio_Reported[i],sd=sqrt(standarderrorfromdata^2))
    tempEvaluateFunction=LogisticModel_PercentUninfected(MeanRandom,SDrandom,ModelParamtemp[,1],ModelParamtemp[,2])
    UpperIntervalsVariant_Predict_Severe$Lower[i]<-100*quantile(tempEvaluateFunction,LowerBound)
    UpperIntervalsVariant_Predict_Severe$Upper[i]<-100*quantile(tempEvaluateFunction,UpperBound)
    UpperIntervalsVariant_Predict_Severe$Efficacy[i]<-100*LogisticModel_PercentUninfected(UpperIntervalsVariant_Predict_Severe$NeutRatio_Reported[i],
                                                                                   SummaryTable_Efficacy_NeutRatio_SD_SEM$PooledSD[1],
                                                                                   tail(FittedLogistic_RawEfficacy_MeanRept_SDPool_Diff_EC50_Only$estimate,3)[1],
                                                                                   tail(FittedLogistic_RawEfficacy_MeanRept_SDPool_Diff_EC50_Only$estimate,3)[2])
  } else {
    MeanRandom=rnorm(N,mean=UpperIntervalsVariant_Predict_Severe$NeutRatio_Reported[i]+censoredChangeMeans_Var$muL_changeFromWTL_Var[censoredChangeMeans_Var$Variant==UpperIntervalsVariant_Predict_Severe$Variant[i]],sd=sqrt(standarderrorfromdata^2+censoredChangeMeans_Var$seL_changeFromWTL_Var[censoredChangeMeans_Var$Variant==UpperIntervalsVariant_Predict_Severe$Variant[i]]^2))
    tempEvaluateFunction=LogisticModel_PercentUninfected(MeanRandom,SDrandom,ModelParamtemp[,1],ModelParamtemp[,2])
    UpperIntervalsVariant_Predict_Severe$Lower[i]<-100*quantile(tempEvaluateFunction,LowerBound)
    UpperIntervalsVariant_Predict_Severe$Upper[i]<-100*(quantile(tempEvaluateFunction,UpperBound))
    UpperIntervalsVariant_Predict_Severe$Efficacy[i]<-100*LogisticModel_PercentUninfected(UpperIntervalsVariant_Predict_Severe$NeutRatio_Reported[i]+censoredChangeMeans_Var$muL_changeFromWTL_Var[censoredChangeMeans_Var$Variant==UpperIntervalsVariant_Predict_Severe$Variant[i]],
                                                                                   SummaryTable_Efficacy_NeutRatio_SD_SEM$PooledSD[1],
                                                                                   tail(FittedLogistic_RawEfficacy_MeanRept_SDPool_Diff_EC50_Only$estimate,3)[1],
                                                                                   tail(FittedLogistic_RawEfficacy_MeanRept_SDPool_Diff_EC50_Only$estimate,3)[2])
  }
}

UpperIntervalsVariant_Predict_Severe$Study="Upper95Predict"
UpperIntervalsVariant_Predict_Severe$Method="RCT"




##Create table of Efficacy Data for Severe
SummaryTable_Efficacy_NeutRatio_SD_SEM_Severe$Method<-"RCT"
SummaryTable_Efficacy_NeutRatio_SD_SEM_Severe$Variant<-"Ancestral"
EfficacyTablewithNeut<-join(Variant_SevereEfficacy,SummaryTable_Efficacy_NeutRatio_SD_SEM_Severe[,c("Study","NeutRatio_Reported","TechnicalName")],by="Study")
Efficacy_Neut_Variant_Severe<-rbind(SummaryTable_Efficacy_NeutRatio_SD_SEM_Severe[,c("Study","TechnicalName","Variant","Method","NeutRatio_Reported","Efficacy")],EfficacyTablewithNeut[,c("Study","TechnicalName","Variant","Method","NeutRatio_Reported","Efficacy")])
Efficacy_Neut_Variant_Severe$TechnicalName[Efficacy_Neut_Variant_Severe$Study=="Sputnik"]="Gam-COVID-Vac"
Efficacy_Neut_Variant_Severe$VariantLabel<-VariantLabel$Label[match(Efficacy_Neut_Variant_Severe$Variant,VariantLabel$Variant)]

######Plot
group.colors=c("black","dodgerblue","mediumpurple","darkorange2","black","turquoise","yellow","green","gray","black","purple")
group.colors.fill=c("black","black","dodgerblue","dodgerblue","white","mediumpurple","white","darkorange2","darkorange2","white","white","black","white","black",NA,"white")
xticks=c(0.0625,0.125,0.25,0.5,1,2,4,8,16,32,64)
ticklabels<-xticks

##Title names and order
UpperIntervalsVariant_Predict_Severe$VariantLabel=factor(UpperIntervalsVariant_Predict_Severe$VariantLabel,levels=c("Ancestral","Alpha (B.1.1.7)","Delta (B.1.617.2)","Beta (B.1.351)"))
Efficacy_Neut_Variant_Severe$VariantLabel=factor(Efficacy_Neut_Variant_Severe$VariantLabel,levels=c("Ancestral","Alpha (B.1.1.7)","Delta (B.1.617.2)","Beta (B.1.351)"))

Figure2b<-ggplot(data=UpperIntervalsVariant_Predict_Severe, aes(x=10^(NeutRatio_Reported),color=VariantLabel)) +
  geom_line(aes(y=Efficacy)) +
  geom_ribbon(aes(ymin=Lower, ymax=Upper,fill=VariantLabel), alpha = 0.15,col=NA,show.legend=FALSE)+
  # geom_rect(aes(xmin=10^(min(BoostingTable$EstimatedLevelAfter)),
  #               xmax=10^(max(BoostingTable$EstimatedLevelAfter)),
  #               ymax=100,
  #               ymin=0),fill="red",alpha=0.01,col=NA) +
  scale_x_log10(breaks=c(xticks),labels=ticklabels) + 
  coord_cartesian(xlim=c(0.17,8)) +
  scale_y_continuous(lim=c(0,100)) +
  theme_linedraw() +
  # geom_vline(xintercept=10^BoostingAverage,linetype=2,size=0.3,color="red") +
  # geom_vline(data=BoostingTable[BoostingTable$PreviousImmune=="Convalescent",],aes(xintercept = 10^(EstimatedLevelAfter),color=PaperNumber),linetype=1,size=0.08,color="red") +
  # geom_vline(data=BoostingTable[BoostingTable$PreviousImmune!="Convalescent",],aes(xintercept = 10^(EstimatedLevelAfter),color=PaperNumber),linetype=1,size=0.12,color="navyblue") +
  geom_point(data=Efficacy_Neut_Variant_Severe,aes(x=10^NeutRatio_Reported,y=100*Efficacy,shape=interaction(Method,TechnicalName),fill=VariantLabel,alpha=TechnicalName,size=TechnicalName)) +
  labs(x="Neutralisation of ancestral (fold of convalescent)\n[International units, IU]",
       y="Efficacy (%)",
       alpha="Vaccine/Serum") +
  scale_color_manual(values=group.colors,guide="none") +
  scale_size_manual(values=c(rep(1.5,3),4,rep(1.5,10)),guide="none") +
  scale_shape_manual(values=list(16,15,0,18,-as.hexmode("002A"),25,24,1,1,1,1),
                      guide="none")+                     
                     # guide=guide_legend(order=1,
                     #                    override.aes = list(size=c(rep(1.5,4),4,rep(1.5,2))))) +
  scale_fill_manual(values=group.colors,guide="none") +
  scale_alpha_manual(values=c(1,1,1,1,1,1,1,1,1),
                     guide=guide_legend(override.aes = list(shape=list(1,0,5,-as.hexmode("002A"),6,2),size=c(1.5,1.5,1.5,4,1.5,1.5)) ) ) +
  theme(axis.line = element_line(colour = "black"),
        # axis.text = element_text(size=2),
        panel.grid.major = element_line(colour = "gray95"),
        panel.grid.minor = element_line(colour = "gray95"),
        panel.border = element_blank(),
        panel.background = element_blank(),
        panel.spacing.y = unit(-1.5, "lines"),
        panel.spacing.x = unit(-1, "lines"),
        legend.spacing.y = unit(0, 'cm'),
        legend.key = element_rect(size = 1),
        legend.key.size = unit(0.45, "cm"),
        legend.margin = margin(t=0,b=0.7,unit="cm"),
        legend.background = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(color="black")
        # plot.margin = margin(b=1,l=1,unit="cm")
  ) +
  facet_rep_wrap(.~VariantLabel,ncol=4)



pdf("Figure2b.pdf",height=2.5,width=9)
Figure2b
dev.off()
Figure2b


####Boosting Figure
xticks=c(0.0625,0.125,0.25,0.5,1,2,4,8,16,32,64)
ticklabels<-xticks

BoostingTable$PreviousImmune<-factor(BoostingTable$PreviousImmune)
levels(BoostingTable$PreviousImmune)<-c("Convalescent","CoronaVac","mRNA-1273")

BoostingTable$BoosterTechnical<-factor(BoostingTable$Booster)
levels(BoostingTable$BoosterTechnical)<-c("CoronaVac","mRNA-1273","BNT162b2","Pfizer/Moderna")

shiftfactor=0.05

Figure3b<-ggplot(data=UpperIntervalsVariant_Predict_Severe, aes(x=10^(NeutRatio_Reported),color=VariantLabel)) +
  # geom_ribbon(aes(ymin=Lower, ymax=Upper,fill=VariantLabel), alpha = 0.15,col=NA,show.legend=FALSE)+
  annotate(geom="rect",
           xmin=10^(min(BoostingTable$EstimatedLevelAfter[BoostingTable$PreviousImmune=="Convalescent"])),
           xmax=10^(max(BoostingTable$EstimatedLevelAfter[BoostingTable$PreviousImmune=="Convalescent"])),
           ymax=100,
           ymin=0,
           fill="red",
           alpha=0.1,
           color=NA) +
  annotate(geom="rect",
           xmin=10^(min(SummaryTable_Efficacy_NeutRatio_SD_SEM$NeutRatio_Reported)),
           xmax=10^(max(SummaryTable_Efficacy_NeutRatio_SD_SEM$NeutRatio_Reported)),
           ymax=100,
           ymin=0,
           fill="black",
           alpha=0.1,
           color=NA) +
  scale_y_continuous() +
  scale_x_log10(breaks=c(xticks),labels=ticklabels) + 
  coord_cartesian(xlim=c(0.2,37),ylim=c(0,103)) +
  theme_linedraw() +
  geom_vline(data=BoostingTable[BoostingTable$PreviousImmune=="Convalescent",],aes(xintercept = 10^(EstimatedLevelAfter)),linetype=1,size=0.08,color="pink") +
  geom_vline(data=BoostingTable[BoostingTable$PreviousImmune!="Convalescent",],aes(xintercept = 10^(EstimatedLevelAfter)),linetype=1,size=0.12,color="navyblue") +
  geom_vline(data=SummaryTable_Efficacy_NeutRatio_SD_SEM,aes(xintercept = 10^(NeutRatio_Reported)),linetype=1,size=0.08,color="gray80") +
  geom_vline(xintercept=10^BoostingAverage,linetype=2,size=0.3,color="red") +
  geom_line(aes(y=Efficacy)) +
  annotate(geom="rect",
           xmin=0.1,
           xmax=37,
           ymin=101,
           ymax=110,
           fill="white",
           col=NA) +
  geom_point(data=SummaryTable_Efficacy_NeutRatio_SD_SEM,
             aes(x=(10^NeutRatio_Reported),
                 y=103,
                 shape=TechnicalName,
                 size=TechnicalName),
             color="black") +
  geom_point(data=BoostingTable,
             aes(x=(10^EstimatedLevelAfter),
                 y=103,
                 shape=PreviousImmune,
                 size=PreviousImmune),
             color="black",
             show.legend=FALSE) +
  geom_point(data=BoostingTable[BoostingTable$BoosterTechnical!="Pfizer/Moderna",],
             aes(x=(10^EstimatedLevelAfter),
                 y=106.5,
                 shape=BoosterTechnical,
                 size=BoosterTechnical),
             color="red",
             show.legend=FALSE) +
  geom_point(data=BoostingTable[BoostingTable$BoosterTechnical=="Pfizer/Moderna",],
             aes(x=(1+shiftfactor)*(10^EstimatedLevelAfter),
                 y=106.5),
             size=0.8,
             shape=22,
             color="red",
             show.legend=FALSE) +
  geom_point(data=BoostingTable[BoostingTable$BoosterTechnical=="Pfizer/Moderna",],
             aes(x=(1-shiftfactor)*(10^EstimatedLevelAfter),
                 y=106.5),
             size=0.8,
             shape=25,
             color="red",
             show.legend=FALSE) +
  labs(x="Neutralisation of ancestral (fold of convalescent)\n[International units, IU]",
       y="Efficacy (%)",
       color="Variant",
       shape="Vaccine/Serum") +
  scale_color_manual(values=group.colors) +
  scale_size_manual(values=c(rep(1.5,5),4,rep(1.5,10)),guide="none") +
  scale_shape_manual(values=list(21,22,23,3,4,-as.hexmode("002A"),25,24,8,8,8,8),
                     guide=guide_legend(order=1,
                                        override.aes = list(size=c(rep(1.5,5),4,rep(1.5,2))))) +
 theme(axis.line = element_line(colour = "black"),
        # axis.text = element_text(size=2),
        panel.grid.major = element_line(colour = "gray95"),
        panel.grid.minor = element_line(colour = "gray95"),
        panel.border = element_blank(),
        panel.background = element_blank(),
        panel.spacing.y = unit(-1.5, "lines"),
        panel.spacing.x = unit(-1, "lines"),
        legend.spacing.y = unit(0, 'cm'),
        legend.key = element_rect(size = 1),
        legend.key.size = unit(0.45, "cm"),
        legend.margin = margin(t=0,b=-0.1,unit="cm"),
        legend.background = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(color="black")
        # plot.margin = margin(b=1,l=1,unit="cm")
  )


pdf("Figure3b.pdf",height=3,width=4.8)
Figure3b
dev.off()
Figure3b








### Logistic model for Raw Efficacy Counts
FittingLogistic_Raw<-function(logRisk0,logk,C50,N_C,N_V,Inf_C,Inf_V,MeanVector,SDVector){
  
  Risk0=exp(logRisk0)
  
  if (length(SDVector)==1) {
    SDVector=rep(SDVector,length(N_C))
  }
  
  IndexNA=(is.na(N_C) | is.na(MeanVector) | is.na(SDVector))
  N_C=N_C[!IndexNA]
  N_V=N_V[!IndexNA]
  Inf_V=Inf_V[!IndexNA]
  Inf_C=Inf_C[!IndexNA]
  MeanVector=MeanVector[!IndexNA]
  SDVector=SDVector[!IndexNA]
  
  if (length(C50)==1) {
    C50=rep(C50,length(N_C))
  }
  
  if (length(logk)==1) {
    logk=rep(logk,length(N_C))
  }
  
  LL=0
  for (i in 1:length(N_C)) {
    
    LL=LL-log(dbinom(Inf_C[i],N_C[i],Risk0[i]))-log(dbinom(Inf_V[i],N_V[i],Risk0[i]*(1-LogisticModel_PercentUninfected(MeanVector[i],SDVector[i],logk[i],C50[i]))))
  }
  LL
}



LogisticEstimate=c("logk"=log(2.7),"C50"=log10(0.5))
FittedLogistic_RawEfficacy_MeanRept_SDPool<-nlm(function(p){FittingLogistic_Raw(p[1:sum(!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM$NumCont))],p[sum(!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM$NumCont))+1],p[sum(!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM$NumCont))+2],SummaryTable_Efficacy_NeutRatio_SD_SEM$NumCont[!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM$NumCont)],SummaryTable_Efficacy_NeutRatio_SD_SEM$NumVac[!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM$NumCont)],SummaryTable_Efficacy_NeutRatio_SD_SEM$InfCont[!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM$NumCont)],SummaryTable_Efficacy_NeutRatio_SD_SEM$InfVac[!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM$NumCont)],
                                                                                SummaryTable_Efficacy_NeutRatio_SD_SEM$NeutRatio_Reported[!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM$NumCont)],SummaryTable_Efficacy_NeutRatio_SD_SEM$PooledSD[!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM$NumCont)])},c(log(SummaryTable_Efficacy_NeutRatio_SD_SEM$InfCont[!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM$NumCont)]/SummaryTable_Efficacy_NeutRatio_SD_SEM$NumCont[!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM$NumCont)]),LogisticEstimate),hessian=TRUE)




#### Prediciting Efficacy
TimeCourse=seq(0,360,by=1)
halflife=108
TimeOfBooster=180

InitialTitre_func=function(mu,targeteff){(log(targeteff)-log(LogisticModel_PercentUninfected(mu,
                                                                                             SummaryTable_Efficacy_NeutRatio_SD_SEM$PooledSD[1],
                                                                                             tail(FittedLogistic_RawEfficacy_MeanRept_SDPool$estimate,2)[1],
                                                                                             tail(FittedLogistic_RawEfficacy_MeanRept_SDPool$estimate,2)[2])))^2
}

InitialEffVector=c(95,90,80,70)
StartingNeut<-NULL
for (i in 1:length(InitialEffVector)) {
  tempfit<-nlm(function(x){InitialTitre_func(x,0.01*InitialEffVector[i])},0)
  StartingNeut[i]<-tempfit$estimate
  
}

VariantsForPlot<-c("Ancestral","B.1.1.7","B.1.351","B.1.617.2")

DecayingProtection=data.frame("Variant"=rep(VariantsForPlot,each=length(TimeCourse)*length(InitialEffVector)),
                             "InitialEff"=rep(rep(InitialEffVector,each=length(TimeCourse)),length(VariantsForPlot)),
                             "InitialNeutAncest"=rep(rep(StartingNeut,each=length(TimeCourse)),length(VariantsForPlot)),
                             "DaysSince"=rep(TimeCourse,length(InitialEffVector)*length(VariantsForPlot)),
                             "Efficacy"=NA,
                             "Boost"=FALSE)

BoostingTabletemp<-DecayingProtection[DecayingProtection$InitialEff==90,]
BoostingTabletemp$Boost<-TRUE

DecayingProtection<-rbind(DecayingProtection,BoostingTabletemp)

DecayingProtection$VariantFoldChange=censoredChangeMeans_Var$muL_changeFromWTL_Var[match(DecayingProtection$Variant,censoredChangeMeans_Var$Variant)]
DecayingProtection$VariantFoldChange[DecayingProtection$Variant=="Ancestral"]=0
DecayingProtection$InitialNeut=DecayingProtection$InitialNeutAncest+DecayingProtection$VariantFoldChange
DecayingProtection$NeutDecay<-log10((10^DecayingProtection$InitialNeut)*exp(-(log(2)/halflife)*DecayingProtection$DaysSince))
DecayingProtection$NeutDecay[DecayingProtection$Boost & DecayingProtection$DaysSince>=TimeOfBooster]=log10(10^(BoostingAverage+DecayingProtection$VariantFoldChange[DecayingProtection$Boost & DecayingProtection$DaysSince>=TimeOfBooster])*exp(-(log(2)/halflife)*(DecayingProtection$DaysSince[DecayingProtection$Boost & DecayingProtection$DaysSince>=TimeOfBooster]-TimeOfBooster)))
DecayingProtection$NeutDecayAncest<-log10((10^DecayingProtection$InitialNeutAncest)*exp(-(log(2)/halflife)*DecayingProtection$DaysSince))
DecayingProtection$NeutDecayAncest[DecayingProtection$Boost & DecayingProtection$DaysSince>=TimeOfBooster]=log10((10^(BoostingAverage))*exp(-(log(2)/halflife)*(DecayingProtection$DaysSince[DecayingProtection$Boost & DecayingProtection$DaysSince>=TimeOfBooster]-TimeOfBooster)))

DecayingProtection$Efficacy<-LogisticModel_PercentUninfected(DecayingProtection$NeutDecay,
                                                            SummaryTable_Efficacy_NeutRatio_SD_SEM$PooledSD[1],
                                                            tail(FittedLogistic_RawEfficacy_MeanRept_SDPool_Diff_EC50_Only$estimate,3)[1],
                                                            tail(FittedLogistic_RawEfficacy_MeanRept_SDPool_Diff_EC50_Only$estimate,3)[2])

DecayingProtection$Lower<-NA
DecayingProtection$VariantLabel<-VariantLabel$Label[match(DecayingProtection$Variant,VariantLabel$Variant)]

for (i in 1:nrow(DecayingProtection)) {
  tempsubset=UpperIntervalsVariant_Predict_Severe[UpperIntervalsVariant_Predict_Severe$Variant==DecayingProtection$Variant[i],]
 
    SmallDiffindx=order(abs(DecayingProtection$NeutDecayAncest[i]-tempsubset$NeutRatio_Reported))[1:2]
    SmallDiffx=tempsubset$NeutRatio_Reported[SmallDiffindx]
    SmallDiffy=tempsubset$Lower[SmallDiffindx]
    tempgradient<-((SmallDiffy[2]-SmallDiffy[1])/(SmallDiffx[2]-SmallDiffx[1]))
    Interpolatey=tempgradient*(DecayingProtection$NeutDecayAncest[i]-SmallDiffx[1])+SmallDiffy[1]
    DecayingProtection$Lower[i]=Interpolatey
  
}


group.colors=c("black","dodgerblue","mediumpurple","darkorange2","black","turquoise","yellow","green","gray","black","purple")
cc1 <- scales::seq_gradient_pal("white", "black", "Lab")(seq(0,1,length.out=length(InitialEffVector)+2))
cc2 <- scales::seq_gradient_pal("white", "dodgerblue", "Lab")(seq(0,1,length.out=length(InitialEffVector)+2))
cc3 <- scales::seq_gradient_pal("white", "mediumpurple", "Lab")(seq(0,1,length.out=length(InitialEffVector)+2))
cc4 <- scales::seq_gradient_pal("white", "darkorange2", "Lab")(seq(0,1,length.out=length(InitialEffVector)+2))
cc<-c(tail(cc1,4),tail(cc2,4),tail(cc3,4),tail(cc4,4))

DecayingProtection$VariantLabel<-factor(DecayingProtection$VariantLabel,levels=c("Ancestral","Alpha (B.1.1.7)","Delta (B.1.617.2)","Beta (B.1.351)"))

Figure4b<-ggplot(data=DecayingProtection,aes(x=DaysSince,y=100*Efficacy,color=interaction(as.factor(InitialEff),VariantLabel))) +
  geom_line(aes(linetype=Boost,alpha=as.factor(InitialEff))) +
  geom_ribbon(data=DecayingProtection[DecayingProtection$Boost==FALSE,],aes(ymin=Lower,ymax=100*Efficacy,fill=interaction(as.factor(InitialEff),VariantLabel)),alpha=0.15,col=NA) +
  geom_ribbon(data=DecayingProtection[DecayingProtection$Boost==TRUE,],aes(ymin=Lower,ymax=100*Efficacy,fill=interaction(as.factor(InitialEff),VariantLabel)),alpha=0.15,col=NA) +
  geom_line(aes(linetype=Boost,alpha=as.factor(InitialEff))) +
  facet_rep_wrap(~VariantLabel,ncol=4) +
  scale_linetype(guide="none") +
  scale_color_manual(values=cc,guide="none") +
  scale_fill_manual(values=cc,guide="none") +
  scale_alpha_manual(values=c(1,1,1,1),guide=guide_legend(override.aes=list(shape=NA,linetype=1,color=tail(cc1,4)))) +
  labs(y="Efficacy",
       x="Time (d)",
       alpha="Initial efficacy\nagainst ancestral") +
  theme_linedraw() +
  theme(axis.line = element_line(colour = "black"),
        # axis.text = element_text(size=2),
        panel.grid.major = element_line(colour = "gray95"),
        panel.grid.minor = element_line(colour = "gray95"),
        panel.border = element_blank(),
        panel.background = element_blank(),
        panel.spacing.y = unit(-0.5, "lines"),
        panel.spacing.x = unit(-1, "lines"),
        # legend.spacing.y = unit(0, 'cm'),
        # legend.key = element_rect(size = 1),
        # legend.key.size = unit(0.45, "cm"),
        # legend.margin = margin(t=0,b=0.7,unit="cm"),
        legend.background = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(color="black")
  )

pdf("Figure4b.pdf",height=2.5,width=10)
Figure4b
dev.off()
Figure4b








