# pronghornsurvival
R code associated with Jones et al. Journal of Wildlife Management paper Annual Pronghorn Survival of a Partially Migratory Population. 

This document provides the R analysis and plotting code for the Bayesian time-to-event analysis used in Jones et al. in review. Note that due to MCMC simulation, results will differ slightly with each run of the same model. We tried to design the code to be flexible so that it can easily handle other time-to-event survival datasets, and we provided a set of custom functions that prepare the data for the JAGS models. However, priors in the actual models may need to be adjusted for the average intensity of mortality observed in your study. All Bayesian models are estimated using JAGS in the package jagsUI (Kellner 2018). Note that you will need to have JAGS installed before you can run the model code. The JAGS open-source software can be downloaded here: https://sourceforge.net/projects/mcmc-jags/files/. The R code runs all of the analyses and makes five figures in Jones et al. (Figures 2-4, Fig. S1, Fig. S2). The script basically follows the results in the paper and will produce all summary statistics and estimates in the text and within the tables. You will need to unzip this folder and open the .R file called: "AppendixA.R" and have JAGS installed. The two data files are stored as .rds files in the folder; one is an annual data set: "prongAnnual.rds", and one is a 4-season data set: "prongSeason.rds". Another .rds file holds the winter severity indices, temperatures, and snow depths for all years of the study and is called: tab.rds (note this is needed to make Figure S1). There is a folder of functions that is kind of like the beginnings of an R package; they are just functions we pieced together or borrowed and we use them throughout the code; note that they will be loaded by running the R code at the beginning, i.e., source("functionLoad.r"). The code will write the JAGS model files as text files to the folder before running each model, so you won't see any of these text files in the folder to begin with (or the figures as you have to create them at the bottom of the code).

http://doi.org/10.5281/zenodo.8475
