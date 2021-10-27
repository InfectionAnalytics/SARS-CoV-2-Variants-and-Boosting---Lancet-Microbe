source('./dropPlot.R')
source('./plotValueByFactor.R')
source('./loadVariantDataFrame.R')
source('./getCensoredMeans.R')
source('./censoredRegression.R')
source('./runRegressionsForAllVariants_VariantNeutL.R')
source('./regressionSummaryHelperFunctions.R')

## Variant Regressions
xVarsUse=list(c(),c('Serum'),
              c('muL_WTneutL_Convalescent_Lab'),c('muL_WTneutL_Convalescent_Lab','Serum'),
              c('WTneutL'),c('WTneutL','Serum'),
              c('WTneutL','muL_WTneutL_Convalescent_Lab'),c('WTneutL','muL_WTneutL_Convalescent_Lab','Serum'),
              c('Laboratory'),c('Laboratory','Serum'),
              c('WTneutL','Laboratory'),c('WTneutL','Laboratory','Serum'))

WTdata=loadVariantDataFrame()

variantNeutTitreDropWTOffsetReg=runRegressionForAllVariants_VariantNeutL(WTdata,yVar='variantneutL', useConvalescent=T,xVarsList = xVarsUse,baseVars=NA, offsetVars=c('WTneutL','muL_changeFromWTL_Var'))
variantNeutTitreDropWTOffsetRegSummary=makeResultSummaryFromRegressionList(variantNeutTitreDropWTOffsetReg)
serumpvalforpaper=compareModels(variantNeutTitreDropWTOffsetRegSummary,11,12,FALSE)
print(paste('Serum pvals for paper:',serumpvalforpaper))
serumpvalforsupp=compareModels(variantNeutTitreDropWTOffsetRegSummary,5,6,FALSE)
print(paste('Serum pvals for supplement:',serumpvalforpaper,serumpvalforsupp))