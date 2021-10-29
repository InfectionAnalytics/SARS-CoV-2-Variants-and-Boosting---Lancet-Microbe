The following files are included in this directory

Data Files
variantNeutDataForImport.csv
manuscriptWTdata.csv
Both are similar, variantNeutDataForImport contains only some of the  columns of manuscriptWTdata

Code Files
Helper Files:
loadGlobalVariantParameters.R - Loads all the parameters needed for analysis and plotting
loadVariantDataFrame.R - loads the data frame used in all plotting and regression
dropPlot.R - Used to make panel B of Figure 1 and Supplementary Figure S2
plotValueByFactor.R - Used to make panel A of Figure 1 and all panels of Figure S1
findXVals.R - helper function to find correct positions along x axis to place data
getCensoredMeans.R - Function to calculate the censored means
censoredRegression.R - Includes functions to run the censored regressions
runRegressionForAllVariants_VariantNeutL.R - Helper functions for regressions
regressionSummaryHelperFunctions.R - more helper functions for regressions

Running Files:
Figure1andSuppFigs1-3.R - Makes the plots
RegressionValues.R - runs the regressions for the manuscript

It is easiest to source all the helper files at the start, however they are also loaded at the beginning of Figure1andSuppFigs1-3.R and RegressionValues.R.

The can be loaded separately using:
source('./loadGlobalVariantParameters.R')
source('./loadVariantDataFrame.R')
source('./getCensoredMeans.R')
source('./censoredRegression.R')
source('./regressionSummaryHelperFunctions.R')
source('./findXVals.R')
source('./dropPlot.R')
source('./plotValueByFactor.R')
source('./runRegressionForAllVariants_VariantNeutL.R')


In addition, to run all the code the following libraries are required:
reshape2,ggplot2,dplyr,stringr,lubridate,zoo,cowplot,ggpubr,latex2exp,doBy,RColorBrewer,ggnewscale,stats

They can be installed via the following commands:
packages_needed = c('reshape2','ggplot2','dplyr','stringr','lubridate','zoo','cowplot','ggpubr','latex2exp','doBy','RColorBrewer','ggnewscale','stats')
install.packages(packages_needed)

and then loaded via:
lapply(packages_needed, require, character.only = TRUE)

They are also loaded (not installed) in loadGlobalVariantParameters.R
