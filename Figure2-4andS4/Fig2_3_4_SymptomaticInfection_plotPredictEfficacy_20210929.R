#Plot efficacy figures (move per variant)
library(stringr)
library(mvtnorm)
library(lemon)
library(ggplot2)
library(magic)
library(plyr)

#### Set directory to same directory as the r-script
currentdirectory<-dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(currentdirectory)

###### Import efficacy data for Variant and WT
#Variant data
Variant_Efficacy<-read.csv(paste0(currentdirectory,"/Variant_Efficacy_20210521.csv"),fileEncoding="UTF-8-BOM")
BoostingTable<-read.csv(paste0(currentdirectory,"/BoostingTableSummary_20210726.csv"),fileEncoding="UTF-8-BOM")
censoredChangeMeans_Var<-read.csv(paste0(currentdirectory,"/censoredChangeMeans_Var_20210812.csv"),fileEncoding="UTF-8-BOM")
  
#WT data (from nature med paper)
SummaryTable_Efficacy_NeutRatio_SD_SEM<-read.csv(paste0(currentdirectory,"/SummaryTable_Efficacy_NeutRatio_SD_SEM.csv"),fileEncoding="UTF-8-BOM")
SummaryTable_Efficacy_NeutRatio_SD_SEM$NeutRatio_Reported=log10(SummaryTable_Efficacy_NeutRatio_SD_SEM$NeutMean/SummaryTable_Efficacy_NeutRatio_SD_SEM$NeutConv)
SummaryTable_Efficacy_NeutRatio_SD_SEM$RatioReported_LB=10^((SummaryTable_Efficacy_NeutRatio_SD_SEM$NeutRatio_Reported)-1.96*SummaryTable_Efficacy_NeutRatio_SD_SEM$SEM)
SummaryTable_Efficacy_NeutRatio_SD_SEM$RatioReported_UB=10^((SummaryTable_Efficacy_NeutRatio_SD_SEM$NeutRatio_Reported)+1.96*SummaryTable_Efficacy_NeutRatio_SD_SEM$SEM)
SummaryTable_Efficacy_NeutRatio_SD_SEM$TechnicalName[SummaryTable_Efficacy_NeutRatio_SD_SEM$TechnicalName==" rAd26-S+rAd5-S"]="Gam-COVID-Vac"


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


### Logistic model for fitting Raw Efficacy Counts
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



#################################################################
#######################################

####Making Table of variant data for efficacy plot
Efficacy_Neut_Variant<-join(Variant_Efficacy,censoredChangeMeans_Var[,c("Variant","muL_changeFromWTL_Var","seL_changeFromWTL_Var")],by=c("Variant"))
colnames(Efficacy_Neut_Variant)[1]="Study"
Efficacy_Neut_Variant<-join(Efficacy_Neut_Variant,SummaryTable_Efficacy_NeutRatio_SD_SEM[,c("Study","NeutRatio_Reported","SEM","TechnicalName")],by="Study")

Efficacy_Neut_Variant$NeutRatio_Reported<-Efficacy_Neut_Variant$NeutRatio_Reported+Efficacy_Neut_Variant$muL_changeFromWTL_Var
Efficacy_Neut_Variant$SEM_Combined<-sqrt(Efficacy_Neut_Variant$SEM^2+Efficacy_Neut_Variant$seL_changeFromWTL_Var^2)
Efficacy_Neut_Variant$RatioReported_UB<-Efficacy_Neut_Variant$NeutRatio_Reported+1.96*Efficacy_Neut_Variant$SEM_Combined
Efficacy_Neut_Variant$RatioReported_LB<-Efficacy_Neut_Variant$NeutRatio_Reported-1.96*Efficacy_Neut_Variant$SEM_Combined

#For plotting when lower CI of efficacy is less than zero set to zero.
Efficacy_Neut_Variant$Lower[Efficacy_Neut_Variant$Lower<0]=0
Efficacy_Neut_Variant$SD_fromWT<-SummaryTable_Efficacy_NeutRatio_SD_SEM$PooledSD[1]
Efficacy_Neut_Variant$SE_PooledSD<-SummaryTable_Efficacy_NeutRatio_SD_SEM$SE_PooledSD[1]

Efficacy_Neut_Variant_temp<-Efficacy_Neut_Variant

colnames(Efficacy_Neut_Variant_temp)[colnames(Efficacy_Neut_Variant_temp)=="NeutRatio_Reported"]="NeutRatio_Reported_Variant"
colnames(Efficacy_Neut_Variant_temp)[colnames(Efficacy_Neut_Variant_temp)=="RatioReported_LB"]="RatioReported_LB_Variant"
colnames(Efficacy_Neut_Variant_temp)[colnames(Efficacy_Neut_Variant_temp)=="RatioReported_UB"]="RatioReported_UB_Variant"

##Add the neut of each vaccine against ancestral as seperate columns - so neut against ancestral
# and predicted neut against variants are listed for each vaccine variant combination
Efficacy_Neut_Variant_temp<-join(Efficacy_Neut_Variant_temp,SummaryTable_Efficacy_NeutRatio_SD_SEM[,c("Study","NeutRatio_Reported","RatioReported_LB","RatioReported_UB")],by="Study")

###Make vaccine efficacy against ancestral table with same structure as variants table.
WTEfficacy<-Efficacy_Neut_Variant_temp[rep(1,nrow(SummaryTable_Efficacy_NeutRatio_SD_SEM)),]
WTEfficacy[,]<-NA
WTEfficacy$Study<-SummaryTable_Efficacy_NeutRatio_SD_SEM$Study
WTEfficacy$Efficacy<-SummaryTable_Efficacy_NeutRatio_SD_SEM$Efficacy
WTEfficacy$NeutRatio_Reported<-SummaryTable_Efficacy_NeutRatio_SD_SEM$NeutRatio_Reported
WTEfficacy$Lower<-SummaryTable_Efficacy_NeutRatio_SD_SEM$Lower/100
WTEfficacy$Upper<-SummaryTable_Efficacy_NeutRatio_SD_SEM$Upper/100
WTEfficacy$RatioReported_LB<-log10(SummaryTable_Efficacy_NeutRatio_SD_SEM$RatioReported_LB)
WTEfficacy$RatioReported_UB<-log10(SummaryTable_Efficacy_NeutRatio_SD_SEM$RatioReported_UB)
WTEfficacy$TechnicalName<-SummaryTable_Efficacy_NeutRatio_SD_SEM$TechnicalName

###Add Covaxin - which was not used for fitting original model.
WTEfficacy[nrow(WTEfficacy)+1,]<-NA
WTEfficacy$Study[nrow(WTEfficacy)]<-"Covaxin"
WTEfficacy$TechnicalName[nrow(WTEfficacy)]<-"Covaxin"
CovaxinNeutMean<-127.6773881
CovaxinNeutConv<-161.1587484
WTEfficacy$NeutRatio_Reported[nrow(WTEfficacy)]<-log10(CovaxinNeutMean/CovaxinNeutConv)
WTEfficacy$Efficacy[nrow(WTEfficacy)]<-0.806

##All Ancestral studies were all RCTs
WTEfficacy$Variant<-"Ancestral"
WTEfficacy$Method<-"RCT"

###Update name for Sputnik Vaccine
WTEfficacy$TechnicalName[WTEfficacy$TechnicalName==" rAd26-S+rAd5-S"]="Gam-COVID-Vac"

###Neut for against the ancestral variant is equal to neut against ancestral variant.
WTEfficacy$NeutRatio_Reported_Variant<-WTEfficacy$NeutRatio_Reported
WTEfficacy$RatioReported_LB_Variant<-WTEfficacy$RatioReported_LB
WTEfficacy$RatioReported_UB_Variant<-WTEfficacy$RatioReported_UB

###Combine WT 
Efficacy_Neut_Variant_Full<-rbind(Efficacy_Neut_Variant_temp,WTEfficacy)




### Fitting Model from Nature Med paper (exactly a reproduction of what was performed in
# Nature Med paper - need the hessian/covariance matrix for confidence bands.
LogisticEstimate=c("logk"=log(2.7),"C50"=log10(0.5))
FittedLogistic_RawEfficacy_MeanRept_SDPool<-nlm(function(p){FittingLogistic_Raw(p[1:sum(!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM$NumCont))],p[sum(!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM$NumCont))+1],p[sum(!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM$NumCont))+2],SummaryTable_Efficacy_NeutRatio_SD_SEM$NumCont[!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM$NumCont)],SummaryTable_Efficacy_NeutRatio_SD_SEM$NumVac[!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM$NumCont)],SummaryTable_Efficacy_NeutRatio_SD_SEM$InfCont[!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM$NumCont)],SummaryTable_Efficacy_NeutRatio_SD_SEM$InfVac[!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM$NumCont)],
                                                                                SummaryTable_Efficacy_NeutRatio_SD_SEM$NeutRatio_Reported[!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM$NumCont)],SummaryTable_Efficacy_NeutRatio_SD_SEM$PooledSD[!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM$NumCont)])},c(log(SummaryTable_Efficacy_NeutRatio_SD_SEM$InfCont[!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM$NumCont)]/SummaryTable_Efficacy_NeutRatio_SD_SEM$NumCont[!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM$NumCont)]),LogisticEstimate),hessian=TRUE)

#####Evaluating the fitted model
NeutValue=10^seq(log10(0.03),log10(10),by=0.1)
Efficacy_Logistic_Raw<-NULL

for (i in 1:length(NeutValue)) {
  
  Efficacy_Logistic_Raw[i]=LogisticModel_PercentUninfected(log10(NeutValue[i]),SummaryTable_Efficacy_NeutRatio_SD_SEM$PooledSD[1],tail(FittedLogistic_RawEfficacy_MeanRept_SDPool$estimate,2)[1],tail(FittedLogistic_RawEfficacy_MeanRept_SDPool$estimate,1))
}

LogisticModel_withPoolSD=data.frame("NeutRatio_Reported"=log10(NeutValue),"Efficacy"=Efficacy_Logistic_Raw)
LogisticModel_withPoolSD$Study<-rep("LogisticModel",length(NeutValue))


### Confidence bounds:
Cov<-solve(FittedLogistic_RawEfficacy_MeanRept_SDPool$hessian)[9:10,9:10]

#################################################### Adding 95% Confidence Intervals #########################################################
LogisticModel_withPoolSD$Variant="WT"

####Using Bootstrap to estimate CIs for the ribbon
N=10000
ModelParamtemp=rmvnorm(N,mean=c(tail(FittedLogistic_RawEfficacy_MeanRept_SDPool$estimate,2)),sigma=Cov)
SDrandom=rnorm(N,mean=SummaryTable_Efficacy_NeutRatio_SD_SEM$PooledSD[1],sd=SummaryTable_Efficacy_NeutRatio_SD_SEM$SE_PooledSD[1])
upperquant<-0.975
lowerquant<-0.025

for (i in 1:nrow(LogisticModel_withPoolSD)) {
  tempEvaluateFunction=LogisticModel_PercentUninfected(LogisticModel_withPoolSD$NeutRatio_Reported[i],SDrandom,ModelParamtemp[,1],ModelParamtemp[,2])
  LogisticModel_withPoolSD$Lower_boot[i]<-100*quantile(tempEvaluateFunction,lowerquant)
  LogisticModel_withPoolSD$Upper_boot[i]<-100*quantile(tempEvaluateFunction,upperquant)
}

LogisticModel_withPoolSD$Method="none"

###Update Variant labels to include "Alpha","Beta", etc naming convention.
Efficacy_Neut_Variant$VariantLabel=VariantLabel$Label[match(Efficacy_Neut_Variant$Variant,VariantLabel$Variant)]
Efficacy_Neut_Variant_Full$VariantLabel=VariantLabel$Label[match(Efficacy_Neut_Variant_Full$Variant,VariantLabel$Variant)]
Efficacy_Neut_Variant_Full$VariantLabel<-factor(Efficacy_Neut_Variant_Full$VariantLabel,levels=c("Ancestral","Alpha (B.1.1.7)","Beta (B.1.351)","Delta (B.1.617.2)"))


###Plot Figure
group.colors=c("gray85","dodgerblue","darkorange2","mediumpurple","black","turquoise","yellow","green","gray","black","purple")
group.colors.fill=c("gray85","dodgerblue","white","darkorange2","white","white","black","white","black",NA)

FigureS4<-ggplot(data=Efficacy_Neut_Variant_Full, aes(y=100*Efficacy,x=10^NeutRatio_Reported_Variant,color=VariantLabel)) +
  geom_point(aes(shape=TechnicalName,fill=interaction(Method,VariantLabel),alpha=Method),size=1) +
  geom_errorbar(data=Efficacy_Neut_Variant_Full[Efficacy_Neut_Variant_Full$Variant=="Ancestral",],aes(ymin=100*Lower,ymax=100*Upper),width=0.025,size=0.3,color="gray85") +
  geom_errorbarh(data=Efficacy_Neut_Variant_Full[Efficacy_Neut_Variant_Full$Variant=="Ancestral",],aes(xmin=10^RatioReported_LB_Variant,xmax=10^RatioReported_UB_Variant),height=1.5,size=0.3,color="gray85") +
  # geom_point(data=Efficacy_Neut_Variant_Full[Efficacy_Neut_Variant_Full$Variant=="Ancestral",],aes(shape=TechnicalName,alpha=Method,size=TechnicalName),fill="gray85",color="gray85") +
  geom_ribbon(data=LogisticModel_withPoolSD,aes(x=10^NeutRatio_Reported, ymin=Lower_boot, ymax=Upper_boot), fill = "black", alpha = 0.15,col=NA)+
  geom_line(data=LogisticModel_withPoolSD,aes(x=10^NeutRatio_Reported),color="black") +
  geom_errorbar(data=Efficacy_Neut_Variant_Full[Efficacy_Neut_Variant_Full$Variant!="Ancestral",],aes(ymin=100*Lower,ymax=100*Upper),width=0.05) +
  geom_errorbarh(data=Efficacy_Neut_Variant_Full[Efficacy_Neut_Variant_Full$Variant!="Ancestral",],aes(xmin=10^RatioReported_LB_Variant,xmax=10^RatioReported_UB_Variant),height=3) +
  geom_point(data=Efficacy_Neut_Variant_Full[Efficacy_Neut_Variant_Full$Variant!="Ancestral",],aes(shape=TechnicalName,fill=interaction(Method,VariantLabel),alpha=Method,size=TechnicalName)) +
  scale_x_log10(limit=c(0.01,11),breaks=c(0.0625/2,0.0625,0.125,0.25,0.5,1,2,4,8),labels=c(0.03125,0.0625,0.125,0.25,0.5,1,2,4,8)) +
  coord_cartesian(xlim=c(0.041,8)) +
  scale_y_continuous(lim=c(0,100)) +
  scale_color_manual(values=group.colors,guide=guide_legend(order=2)) +
  scale_size_manual(values=c(rep(1.5,6),4,rep(1.5,2)),guide="none") +
  scale_shape_manual(values=list(21,22,23,3,4,8,-as.hexmode("002A"),25,24),
                     guide=guide_legend(order=1,
                                        override.aes = list(fill=NA,size=c(rep(1.5,6),4,rep(1.5,2))))) +
  scale_fill_manual(values=group.colors.fill,guide="none") +
  scale_alpha_manual(values=c(1,1),
                     guide=guide_legend(order=3,override.aes = list(shape=c(16,1),linetype=c(0,0)) ) ) +
  theme_linedraw() +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_line(colour = "gray95"),
        panel.grid.minor = element_line(colour = "gray95"),
        panel.border = element_blank(),
        panel.background = element_blank(),
        legend.spacing.y = unit(0, 'cm'),
        legend.key = element_rect(size = 1),
        legend.key.size = unit(0.45, "cm"),
        legend.margin = margin(t=0.8,b=-0.4,unit="cm"),
        legend.background = element_blank()) +
  labs(x="Neutralisation level against ancestral\n(/convalescent plasma)",
       y="Reported efficacy (%)",
       shape="Vaccine/Serum",
       alpha="Study Design",
       color="Variant")


pdf("FigureS4.pdf",height=3.5,width=5.5)
FigureS4
dev.off()
FigureS4




###Correlation of neut and efficacy.
VariantCor<-cor.test(Efficacy_Neut_Variant$NeutRatio_Reported,Efficacy_Neut_Variant$Efficacy,method="spearman")
VariantCor



##############Using Model to predict Efficacy against each variant ##########
##### Using fold change in neut for each variant and model from Nature Med study
# we predict efficacy for each variant with Confidence Bounds.


ActualNeuts<-data.frame("Variant"=unique(Efficacy_Neut_Variant_Full[,c("Variant","NeutRatio_Reported_Variant")])[,1],
                        "ActualNeuts"=unique(Efficacy_Neut_Variant_Full[,c("Variant","NeutRatio_Reported_Variant")])[,2],
                        "LowerBound"=NA,
                        "UpperBound"=NA)
NeutValuelog10=seq(log10(0.0625),log10(70),by=0.1)
N=10000
ModelParamtemp=rmvnorm(N,mean=c(tail(FittedLogistic_RawEfficacy_MeanRept_SDPool$estimate,2)),sigma=Cov)
SDrandom=rnorm(N,mean=SummaryTable_Efficacy_NeutRatio_SD_SEM$PooledSD[1],sd=SummaryTable_Efficacy_NeutRatio_SD_SEM$SE_PooledSD[1])
standarderrorfromdata=max(SummaryTable_Efficacy_NeutRatio_SD_SEM$SEM)

UpperIntervalsVariant<-data.frame("Variant"=rep(c(unique(Efficacy_Neut_Variant$Variant),"Ancestral"),each=length(NeutValuelog10)),
                                  "NeutRatio_Reported"=rep(NeutValuelog10,length(unique(Efficacy_Neut_Variant$Variant))+1),
                                  "Lower"=NA)

UpperIntervalsVariant$ActualNeut<-UpperIntervalsVariant$NeutRatio_Reported+pmin(0,censoredChangeMeans_Var$muL_changeFromWTL_Var[match(UpperIntervalsVariant$Variant,censoredChangeMeans_Var$Variant)],na.rm=TRUE)
UpperIntervalsVariant$Efficacy=100*LogisticModel_PercentUninfected(UpperIntervalsVariant$ActualNeut,SummaryTable_Efficacy_NeutRatio_SD_SEM$PooledSD[1],tail(FittedLogistic_RawEfficacy_MeanRept_SDPool$estimate,2)[1],tail(FittedLogistic_RawEfficacy_MeanRept_SDPool$estimate,2)[2])


#### Bootstrap CIs of efficacy for each variant
# Upper95Interval<-NULL
UpperBound=0.975
LowerBound=0.025
for (i in 1:nrow(UpperIntervalsVariant)) {
  if (UpperIntervalsVariant$Variant[i]=="Ancestral") {
    MeanRandom=rnorm(N,mean=UpperIntervalsVariant$NeutRatio_Reported[i],sd=standarderrorfromdata)  
    tempEvaluateFunction=LogisticModel_PercentUninfected(MeanRandom,SDrandom,ModelParamtemp[,1],ModelParamtemp[,2])
    UpperIntervalsVariant$Lower[i]<-100*quantile(tempEvaluateFunction,LowerBound)
    UpperIntervalsVariant$Upper[i]<-100*quantile(tempEvaluateFunction,UpperBound)
  } else {
  MeanRandom=rnorm(N,mean=UpperIntervalsVariant$NeutRatio_Reported[i]+censoredChangeMeans_Var$muL_changeFromWTL_Var[censoredChangeMeans_Var$Variant==UpperIntervalsVariant$Variant[i]],sd=sqrt(standarderrorfromdata^2+censoredChangeMeans_Var$seL_changeFromWTL_Var[censoredChangeMeans_Var$Variant==UpperIntervalsVariant$Variant[i]]^2))
  tempEvaluateFunction=LogisticModel_PercentUninfected(MeanRandom,SDrandom,ModelParamtemp[,1],ModelParamtemp[,2])
  UpperIntervalsVariant$Lower[i]<-100*quantile(tempEvaluateFunction,LowerBound)
  UpperIntervalsVariant$Upper[i]<-100*quantile(tempEvaluateFunction,UpperBound)
}
}

for (i in 1:nrow(ActualNeuts)) {
if (ActualNeuts$Variant[i]=="Ancestral") {
  MeanRandom=rnorm(N,mean=ActualNeuts$ActualNeuts[i],sd=sqrt(standarderrorfromdata^2))
  tempEvaluateFunction=LogisticModel_PercentUninfected(MeanRandom,SDrandom,ModelParamtemp[,1],ModelParamtemp[,2])
  ActualNeuts$LowerBound[i]<-100*quantile(tempEvaluateFunction,LowerBound)
  ActualNeuts$UpperBound[i]<-100*quantile(tempEvaluateFunction,UpperBound)
  
} else {
    MeanRandom=rnorm(N,mean=ActualNeuts$ActualNeuts[i],sd=sqrt(standarderrorfromdata^2+censoredChangeMeans_Var$seL_changeFromWTL_Var[censoredChangeMeans_Var$Variant==ActualNeuts$Variant[i]]^2))
    tempEvaluateFunction=LogisticModel_PercentUninfected(MeanRandom,SDrandom,ModelParamtemp[,1],ModelParamtemp[,2])
    ActualNeuts$LowerBound[i]<-100*quantile(tempEvaluateFunction,LowerBound)
    ActualNeuts$UpperBound[i]<-100*quantile(tempEvaluateFunction,UpperBound)
}
  
}

Efficacy_Neut_Variant_Full$PredictedLower<-ActualNeuts$LowerBound[match(Efficacy_Neut_Variant_Full$NeutRatio_Reported_Variant,ActualNeuts$ActualNeuts)]
Efficacy_Neut_Variant_Full$PredictedUpper<-ActualNeuts$UpperBound[match(Efficacy_Neut_Variant_Full$NeutRatio_Reported_Variant,ActualNeuts$ActualNeuts)]
Efficacy_Neut_Variant_Full$PredictedEfficacy<-100*LogisticModel_PercentUninfected(Efficacy_Neut_Variant_Full$NeutRatio_Reported_Variant,SummaryTable_Efficacy_NeutRatio_SD_SEM$PooledSD[1],tail(FittedLogistic_RawEfficacy_MeanRept_SDPool$estimate,2)[1],tail(FittedLogistic_RawEfficacy_MeanRept_SDPool$estimate,2)[2])



###Setting some variables for plotting.
UpperIntervalsVariant$Study="Upper95Predict"
UpperIntervalsVariant$Method="RCT"
UpperIntervalsVariant_Predict<-UpperIntervalsVariant
UpperIntervalsVariant_Predict$Study="Upper95Predict"
UpperIntervalsVariant_Predict$Method="RCT"
UpperIntervalsVariant_Predict$VariantLabel=VariantLabel$Label[match(UpperIntervalsVariant_Predict$Variant,VariantLabel$Variant)]
UpperIntervalsVariant_Predict$VariantLabel<-factor(UpperIntervalsVariant_Predict$VariantLabel,levels=c("Ancestral","Alpha (B.1.1.7)","Delta (B.1.617.2)","Beta (B.1.351)"))
UpperIntervalsVariant_Predict$FillIndex=interaction(UpperIntervalsVariant_Predict$Method,UpperIntervalsVariant_Predict$VariantLabel)
UpperIntervalsVariant_Predict$FillIndex=factor(UpperIntervalsVariant_Predict$FillIndex,
                                               levels=c("RCT.Ancestral",
                                                        "RCT.Alpha (B.1.1.7)",
                                                        "RCT.Delta (B.1.617.2)",
                                                        "RCT.Beta (B.1.351)"))
Efficacy_Neut_Variant_Full$FillIndex=interaction(Efficacy_Neut_Variant_Full$Method,Efficacy_Neut_Variant_Full$VariantLabel)
Efficacy_Neut_Variant_Full$FillIndex=factor(Efficacy_Neut_Variant_Full$FillIndex,
                                            levels=c("RCT.Ancestral",
                                                     "RCT.Alpha (B.1.1.7)",
                                                     "RCT.Delta (B.1.617.2)",
                                                     "RCT.Beta (B.1.351)",
                                                     "TNCC.Alpha (B.1.1.7)",
                                                     "TNCC.Delta (B.1.617.2)",
                                                     "TNCC.Beta (B.1.351)"))


#####Plot
xticks=c(0.0625,0.125,0.25,0.5,1,2,4,8,16,32,64)
ticklabels<-xticks

group.colors=c("black","dodgerblue","mediumpurple","darkorange2","black","turquoise","yellow","green","gray","black","purple")
group.colors.fill=c("dodgerblue","black","darkorange2","mediumpurple","NA","NA","NA","NA","NA","NA","black",NA,"white")

Figure2a<-ggplot(data=UpperIntervalsVariant_Predict, aes(x=10^(NeutRatio_Reported),color=VariantLabel,fill=FillIndex)) +
  #adding the bands
  geom_line(aes(y=Efficacy)) +
  geom_ribbon(aes(ymin=Lower,
                  ymax=Upper),
              alpha = 0.15,
              col=NA,
              show.legend=FALSE) +
  scale_y_continuous(lim=c(0,100)) +
  scale_x_log10(breaks=c(xticks),labels=ticklabels) + 
  coord_cartesian(xlim=c(0.17,8)) +
  geom_point(data=Efficacy_Neut_Variant_Full,
             aes(x=10^NeutRatio_Reported,
                 y=100*Efficacy,
                 shape=TechnicalName,
                 size=TechnicalName,
                 alpha=Method)) +
  labs(x="Neutralisation of ancestral (fold of convalescent)",
       y="Efficacy (%)",
       shape="Vaccine/Serum",
       alpha="Study Design") +
  scale_color_manual(values=group.colors,guide="none") +
  scale_size_manual(values=c(rep(1.5,6),4,rep(1.5,2)),guide="none") +
  scale_shape_manual(values=list(21,22,23,3,4,8,-as.hexmode("002A"),25,24),
                     guide=guide_legend(order=1,
                                        override.aes = list(size=c(rep(1.5,6),4,rep(1.5,2))))) +
  scale_fill_manual(values=group.colors.fill,guide="none") +
  scale_alpha_manual(values=c(1,1),
                     guide=guide_legend(order=3,override.aes = list(shape=c(16,1),linetype=c(0,0)) ) ) +
  theme_linedraw() +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_line(colour = "gray95"),
        panel.grid.minor = element_line(colour = "gray95"),
        panel.border = element_blank(),
        panel.background = element_blank(),
        panel.spacing.y = unit(-1.5, "lines"),
        panel.spacing.x = unit(-1, "lines"),
        legend.spacing.y = unit(0, 'cm'),
        legend.key = element_rect(size = 0),
        legend.key.size = unit(0.45, "cm"),
        legend.margin = margin(b=0,unit="cm"),
        legend.background = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(color="black")
        ) +
  facet_rep_wrap(.~VariantLabel,ncol=4)
  

pdf("Figure2a.pdf",height=2.5,width=9)
Figure2a
dev.off()
Figure2a




###Boosting level data - for plot with boosting.
ModernaLevelPreBoost=SummaryTable_Efficacy_NeutRatio_SD_SEM$NeutRatio_Reported[SummaryTable_Efficacy_NeutRatio_SD_SEM$Study=="Moderna"]
ModernaLevelPreBoost_Pf=SummaryTable_Efficacy_NeutRatio_SD_SEM$NeutRatio_Reported[SummaryTable_Efficacy_NeutRatio_SD_SEM$Study=="Pfizer"]
AverageLevel=0.5*(ModernaLevelPreBoost_Pf+ModernaLevelPreBoost)
BoostingTable$ReferenceNormalisation<-SummaryTable_Efficacy_NeutRatio_SD_SEM$NeutRatio_Reported[match(BoostingTable$Reference,SummaryTable_Efficacy_NeutRatio_SD_SEM$Study)]
BoostingTable$ReferenceNormalisation[BoostingTable$Reference=="Pfizer/Moderna"]=AverageLevel
BoostingTable$EstimatedLevelAfter<-log10((10^BoostingTable$ReferenceNormalisation)*BoostingTable$FoldChange)

BoostingAverage=mean(BoostingTable$EstimatedLevelAfter[BoostingTable$PreviousImmune=="Convalescent"])
BoostingRatio_SE=sd(BoostingTable$EstimatedLevelAfter[BoostingTable$PreviousImmune=="Convalescent"])/sqrt(length(BoostingTable$EstimatedLevelAfter[BoostingTable$PreviousImmune=="Convalescent"]))



BoostingTable$PreviousImmune<-factor(BoostingTable$PreviousImmune)
levels(BoostingTable$PreviousImmune)<-c("Convalescent","CoronaVac","mRNA-1273")

BoostingTable$BoosterTechnical<-factor(BoostingTable$Booster)
levels(BoostingTable$BoosterTechnical)<-c("CoronaVac","mRNA-1273","BNT162b2","Pfizer/Moderna")


###Plotting boosting.
xticks=c(0.0625,0.125,0.25,0.5,1,2,4,8,16,32,64)
ticklabels<-xticks

group.colors=c("black","dodgerblue","mediumpurple","darkorange2","black","turquoise","yellow","green","gray","black","purple")
group.colors.fill=c("dodgerblue","dodgerblue","white","black","black","darkorange2","darkorange2","white","mediumpurple","white","white","black","white","black",NA,"white")

shiftfactor=0.05

Figure3a<-ggplot(data=UpperIntervalsVariant_Predict) +
  #adding the bands
  annotate(geom = "rect", 
           xmin=10^(min(BoostingTable$EstimatedLevelAfter[BoostingTable$PreviousImmune=="Convalescent"])),
           xmax=10^(max(BoostingTable$EstimatedLevelAfter[BoostingTable$PreviousImmune=="Convalescent"])),
           ymax=100,
           ymin=0,
           fill="red",
           alpha=0.2,
           color=NA) +
  annotate(geom = "rect", 
           xmin=10^(min(Efficacy_Neut_Variant_Full$NeutRatio_Reported[Efficacy_Neut_Variant_Full$Variant=="Ancestral"])),
           xmax=10^(max(Efficacy_Neut_Variant_Full$NeutRatio_Reported[Efficacy_Neut_Variant_Full$Variant=="Ancestral"])),
           ymax=100,
           ymin=0,
           fill="black",
           alpha=0.1,
           color=NA) +
  geom_vline(xintercept=10^BoostingAverage,linetype=2,size=0.3,color="red") +
  geom_vline(data=BoostingTable[BoostingTable$PreviousImmune=="Convalescent",],aes(xintercept = 10^(EstimatedLevelAfter)),linetype=1,size=0.08,color="pink") +
  geom_vline(data=BoostingTable[BoostingTable$PreviousImmune!="Convalescent",],aes(xintercept = 10^(EstimatedLevelAfter)),linetype=1,size=0.12,color="navyblue") +
  geom_vline(data=Efficacy_Neut_Variant_Full[Efficacy_Neut_Variant_Full$Variant=="Ancestral",],aes(xintercept = 10^(NeutRatio_Reported)),linetype=1,size=0.08,color="gray80") +
  geom_line(aes(x=10^(NeutRatio_Reported),y=Efficacy,color=VariantLabel)) +
  # geom_line(aes(x=10^(NeutRatio_Reported),y=Lower,color=Variant),size=0.08) +
  # geom_ribbon(aes(x=10^(NeutRatio_Reported),ymin=Lower, ymax=Upper,fill=Variant), alpha = 0.15,col=NA,show.legend=FALSE)+
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
  scale_y_continuous() +
  scale_x_log10(breaks=c(xticks),labels=ticklabels) + 
  coord_cartesian(xlim=c(0.2,37),ylim=c(0,103)) +
  labs(x="Neutralisation of ancestral (fold of convalescent)\n[International units, IU]",
       y="Efficacy (%)",
       color="Variant") +
  labs(x="Neutralisation of ancestral (fold of convalescent)\n[International units, IU]",
       y="Efficacy (%)",
       color="Variant",
       shape="Vaccine/Serum") +
  scale_color_manual(values=group.colors) +
  scale_size_manual(values=c(rep(1.5,5),4,rep(1.5,10)),guide="none") +
  scale_shape_manual(values=list(21,22,23,3,4,-as.hexmode("002A"),25,24,8,8,8,8),
                     guide=guide_legend(order=1,
                                        override.aes = list(size=c(rep(1.5,5),4,rep(1.5,2))))) +
  theme_linedraw() +
  theme(axis.line = element_line(colour = "black"),
        # axis.text = element_text(size=2),
        panel.grid.major = element_line(colour = "gray95"),
        panel.grid.minor = element_line(colour = "gray95"),
        panel.border = element_blank(),
        panel.background = element_blank(),
        panel.spacing.y = unit(-1.5, "lines"),
        panel.spacing.x = unit(-1, "lines"),
        legend.spacing.y = unit(0, 'cm'),
        legend.key = element_rect(size = 0),
        legend.key.size = unit(0.45, "cm"),
        legend.margin = margin(b=0,unit="cm"),
        legend.background = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(color="black")
        # plot.margin = margin(b=1,l=1,unit="cm")
  ) 


pdf("Figure3a.pdf",height=3,width=4.8)
Figure3a
dev.off()
Figure3a





#### Prediciting Efficacy (Figure 4)
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
                                                            tail(FittedLogistic_RawEfficacy_MeanRept_SDPool$estimate,2)[1],
                                                            tail(FittedLogistic_RawEfficacy_MeanRept_SDPool$estimate,2)[2])

DecayingProtection$Lower<-NA
DecayingProtection$VariantLabel<-VariantLabel$Label[match(DecayingProtection$Variant,VariantLabel$Variant)]


####Interpolating CI bounds for decay curve based on bootstrap CIs above.
for (i in 1:nrow(DecayingProtection)) {
  tempsubset=UpperIntervalsVariant_Predict[UpperIntervalsVariant_Predict$Variant==DecayingProtection$Variant[i],]

  SmallDiffindx=order(abs(DecayingProtection$NeutDecayAncest[i]-tempsubset$NeutRatio_Reported))[1:2]
  SmallDiffx=tempsubset$NeutRatio_Reported[SmallDiffindx]
  SmallDiffy=tempsubset$Lower[SmallDiffindx]
  SmallDiffy_UB=tempsubset$Upper[SmallDiffindx]
  tempgradient<-((SmallDiffy[2]-SmallDiffy[1])/(SmallDiffx[2]-SmallDiffx[1]))
  Interpolatey=tempgradient*(DecayingProtection$NeutDecayAncest[i]-SmallDiffx[1])+SmallDiffy[1]
  Interpolatey_UB=tempgradient*(DecayingProtection$NeutDecayAncest[i]-SmallDiffx[1])+SmallDiffy_UB[1]
  DecayingProtection$Lower[i]=Interpolatey
  DecayingProtection$Upper[i]=Interpolatey_UB

}


###plotting decay.

group.colors=c("black","dodgerblue","mediumpurple","darkorange2","black","turquoise","yellow","green","gray","black","purple")
cc1 <- scales::seq_gradient_pal("white", "black", "Lab")(seq(0,1,length.out=length(InitialEffVector)+2))
cc2 <- scales::seq_gradient_pal("white", "dodgerblue", "Lab")(seq(0,1,length.out=length(InitialEffVector)+2))
cc3 <- scales::seq_gradient_pal("white", "mediumpurple", "Lab")(seq(0,1,length.out=length(InitialEffVector)+2))
cc4 <- scales::seq_gradient_pal("white", "darkorange2", "Lab")(seq(0,1,length.out=length(InitialEffVector)+2))
cc<-c(tail(cc1,4),tail(cc2,4),tail(cc3,4),tail(cc4,4))
  
DecayingProtection$VariantLabel<-factor(DecayingProtection$VariantLabel,levels=c("Ancestral","Alpha (B.1.1.7)","Delta (B.1.617.2)","Beta (B.1.351)"))

Figure4a<-ggplot(data=DecayingProtection,aes(x=DaysSince,y=100*Efficacy,color=interaction(as.factor(InitialEff),VariantLabel))) +
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

pdf("Figure4a.pdf",height=2.5,width=10)
Figure4a
dev.off()
Figure4a




