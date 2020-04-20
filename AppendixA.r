

#  Annual Pronghorn Survival of a Partially Migratory Population
#  Jones P. F., A. F. Jakes, D. R. Eacker, and M. Hebblewhite
#  Supporting Information
#  7 April, 2020


# Appendix A - Run summary statistics and analysis for pronghorn survival, cause-specific mortality, and risk regression


# * Notes *
# -One important difference in our BUGS model specification compared to the OpenBUGS example (http://www.openbugs.net/Examples/Leuk.html)
# is that we create the data before running the model in JAGS, whereas their example includes a data loop in the model to create the counting
# process data. We considered it a cleaner presentation to house this data creation step within R functions rather than in the model,
# but the results are identical between the two approaches.
#
# -Also, note that results will differ ever so slightly with each MCMC run since simulation is used to estimate the posterior distribution of parameters.


# BEGIN R CODE
source("functionLoad.R") # source .r file to load all functions

# install and load R packages using ipak function
pks = c("jagsUI","ggplot2","scales","dplyr","MASS","survival","survminer","viridis","bshazard","knitr","RColorBrewer")
ipak(pks) 

# read in data
pa = readRDS("prongAnnual.rds") # annual time-to-event dataset 
str(pa) # look at structure of data

# read in 4 season dataset
ps = readRDS("prongSeason.rds") # 4-season time-to-event dataset 
str(ps) # look at structure of data


# 1) BASIC DATA SUMMARIES

# total number of observations
nrow(pa)

# number of unique individuals in sample
length(unique(pa$id))

# Breakdown of migrants vs. residents within years (with individuals stacked across years)
with(pa, table(year, tactic))

# Breakdown of migrants vs. residents across years (with individuals stacked across years)
with(pa, table(tactic))

# number of observations with unknown tactic
with(pa, length(which(is.na(tactic)==TRUE)))

# individual level migration summaries
d=with(pa, table(id,tactic)) # store table of id and tactic frequencies
length(which(d[,1]>=1 & d[,2]>=1)) # number of individuals that switched tactics across years
length(which(d[,1]==0 & d[,2]==0)) # number of individuals with unknown migration status
length(which(d[,1]>=1 & d[,2]==0)) # number of individuals that were resident across all years
length(which(d[,1]==0 & d[,2]>=1)) # number of individuals that were migratory across all years

# number of years each individual was monitored
temp=with(pa, data.frame(id,diff=(exit-enter)/365.25))
d4=data.frame(temp %>% group_by(id) %>% summarise(tot=sum(diff))) # use dplyr package
sum(d4$tot) # total number of pronghorn risk-years
median(d4$tot) # median number of years at risk per individual
sd(d4$tot) # standard deviation of number of years at risk per individual



# 2) SURVIVORSHIP ESTIMATES

# 2a) Annual survival estimate
jags.data=prep.bsurv(data=pa, type="s") # use prep.bsurv function to prep data for survival estimates only (type="s") 

# set parameters to track
params <- c("S0","cum.haz") # parameters to track (here,S0=baseline survival function, cum.haz=baseline hazard function)

# specify model
sink("Bayesian_surv.txt") # JAGS version
cat("
# begin model loop
 model{
    for(j in 1:n.fail){
	for(i in 1 :n.ind){								
	   dN[i,j] ~ dpois(Idt[i,j]) # Poisson Likelihood
	   Idt[i,j] = Y[i,j] * h0[j]
	} # i	
     mu[j] = r * (t[j + 1] - t[j]) * c # prior mean hazard
	   h0[j] ~ dgamma(mu[j], c) # hazard rate
     S0[j] = exp(-cum.haz[j]) # baseline survivorship function
     cum.haz[j] = sum(h0[1:j]) # cumulative hazard function
    } # j
   # priors
    c = 0.001 # degree of confidence in guess of hazard rate
    r = 0.003 # guess at the daily hazard rate
# end model
}	
    ",fill = TRUE)
sink()

# run autojags() function from jagsUI package to get surivorship estimates and save as "outS"
outS=autojags(model.file="Bayesian_surv.txt",data=jags.data,parameters.to.save=params,parallel=TRUE,
             iter.increment=1000,n.chains=2, n.burnin=1000,n.adapt=1000,max.iter=100000,n.thin=1)

# get quantiles of survival estimate
round(quantile(outS$sims.list$S0[,ncol(outS$sims.list$S0)],c(0.5,0.025,0.975)),2)


# 2b) Annual survival difference bewteen migrants and residents

# set resident to reference group
pa$tactic = relevel(pa$tactic, ref="Resident")

split.data = split(pa, pa$tactic) # split data into a list for migratory and resident pronghorn (this will drop NA's for tactic)

# now use prep.bsurv fuction to create data for both tactics 
jags.data.list = lapply(split.data, function(x) prep.bsurv(data=x, type="s"))

# run autojags() function from jagsUI package to get surivorship estimates and save as "outTactic" for both tactics
outTactic=lapply(jags.data.list,function(x) autojags(model.file="Bayesian_surv.txt",data=x,parameters.to.save=params,parallel=TRUE,iter.increment=1000,n.chains=2, n.burnin=1000,n.adapt=1000,max.iter=100000,n.thin=1))

# annual survival for migrants
round(quantile(outTactic[[2]]$sims.list$S0[,ncol(outTactic[[2]]$sims.list$S0)],c(0.5,0.025,0.975)),2)

# annual survival for residents
round(quantile(outTactic[[1]]$sims.list$S0[,ncol(outTactic[[1]]$sims.list$S0)],c(0.5,0.025,0.975)),2)

# difference in survival between migrants and residents
round(quantile(outTactic[[2]]$sims.list$S0[,ncol(outTactic[[2]]$sims.list$S0)]-outTactic[[1]]$sims.list$S0[,ncol(outTactic[[1]]$sims.list$S0)],c(0.5,0.025,0.975)),2)


# 2c) Table 1 results - Annual yearly survival for resident, migratory, annual pooled (included individuals with unknown migratory status) 

# produce multi-year survival estimates for migrants
sub.data = pa[pa$tactic=="Migratory",] # first subset the data

split.data = split(sub.data, sub.data$year) # now split the data into a list baesd on year

(year.index=unlist(lapply(split.data, function(x) sum(x$event)))) # look at sum of events in each year

years=unique(na.omit(sub.data$year))# store years for table below

# remove years when no events occur (i.e., Survival = 1.00 )
  if(any(year.index==0)){
    split.data=split.data[which(year.index>0)] # remove years when no events occur (i.e., Survival = 1.00 )
    years=unique(na.omit(sub.data$year))[which(year.index>0)] # store years for table below
  }

# use prep.bsurv function to get data sets for migrants for all years
jags.data.list = lapply(split.data, function(x) prep.bsurv(data=x, type="s"))

# run jags model for all years by applying autojags to the data list using lapply
all.years = lapply(jags.data.list, function(x) autojags(model.file="Bayesian_surv.txt",data=x,parameters.to.save=params,parallel=TRUE,
             iter.increment=1000,n.chains=2, n.burnin=1000,n.adapt=1000,max.iter=100000,n.thin=1))

names(all.years)=years # name jags output list

out=data.frame(year=years,survival=unlist(lapply(all.years, function(x) tail(x$q50$S0,1))),
           lowerBCI=unlist(lapply(all.years, function(x) tail(x$q2.5$S0,1))),
           upperBCI=unlist(lapply(all.years, function(x) tail(x$q97.5$S0,1))))
# adjust survival in 2011 for the fact that it only included 109 days (use adj.surv() function)
out[which(out$year==2011),2:4]=adj.surv(surv=out[which(out$year==2011),2:4],days=109,totDays=365)
out[,2:4]=apply(out[,2:4],2,function(x)round(x,2))
mig.out=out # store migratory estimates

mig.out

# produce multi-year survival estimates for residents
sub.data = pa[pa$tactic=="Resident",] # first subset the data

split.data = split(sub.data, sub.data$year) # now split the data into a list baesd on year

(year.index=unlist(lapply(split.data, function(x) sum(x$event)))) # look at sum of events in each year

years=unique(na.omit(sub.data$year))# store years for table below

# remove years when no events occur (i.e., Survival = 1.00 )
  if(any(year.index==0)){
    split.data=split.data[which(year.index>0)] # remove years when no events occur (i.e., Survival = 1.00 )
    years=unique(na.omit(sub.data$year))[which(year.index>0)] # store years for table below
  }

# use prep.bsurv function to get data sets for residents for all years
jags.data.list = lapply(split.data, function(x) prep.bsurv(data=x, type="s"))

# run jags model for all years by applying autojags to the data list using lapply
all.years = lapply(jags.data.list, function(x) autojags(model.file="Bayesian_surv.txt",data=x,parameters.to.save=params,parallel=TRUE,iter.increment=1000,n.chains=2,n.burnin=1000,n.adapt=1000,max.iter=100000,n.thin=1))

names(all.years)=years # name jags output list

out=data.frame(year=years,survival=unlist(lapply(all.years, function(x) tail(x$q50$S0,1))),
           lowerBCI=unlist(lapply(all.years, function(x) tail(x$q2.5$S0,1))),
           upperBCI=unlist(lapply(all.years, function(x) tail(x$q97.5$S0,1))))
# adjust survival in 2011 for the fact that it only included 109 days (use adj.surv() function)
out[which(out$year==2011),2:4]=adj.surv(surv=out[which(out$year==2011),2:4],days=109,totDays=365)
out[,2:4]=apply(out[,2:4],2,function(x)round(x,2))

res.out=out # store resident estimates

res.out

# produce multi-year survival estimates for all data pooled
sub.data = pa # use the entire annual dataset

split.data = split(sub.data, sub.data$year) # now split the data into a list baesd on year

(year.index=unlist(lapply(split.data, function(x) sum(x$event)))) # look at sum of events in each year

years=unique(na.omit(sub.data$year))# store years for table below

# remove years when no events occur (i.e., Survival = 1.00 )
  if(any(year.index==0)){
    split.data=split.data[which(year.index>0)] # remove years when no events occur (i.e., Survival = 1.00 )
    years=unique(na.omit(sub.data$year))[which(year.index>0)] # store years for table below
  }

# use prep.bsurv function to get data sets for all years with individuals pooled
jags.data.list = lapply(split.data, function(x) prep.bsurv(data=x, type="s"))

# run jags model for all years by applying autojags to the data list using lapply
all.years = lapply(jags.data.list, function(x) autojags(model.file="Bayesian_surv.txt",data=x,parameters.to.save=params,parallel=TRUE,
             iter.increment=1000,n.chains=2, n.burnin=1000,n.adapt=1000,max.iter=100000,n.thin=1))

names(all.years)=years # name jags output list

out=data.frame(year=years,survival=unlist(lapply(all.years, function(x) tail(x$q50$S0,1))),
           lowerBCI=unlist(lapply(all.years, function(x) tail(x$q2.5$S0,1))),
           upperBCI=unlist(lapply(all.years, function(x) tail(x$q97.5$S0,1))))
# adjust survival in 2011 for the fact that it only included 109 days (use adj.surv() function)
out[which(out$year==2011),2:4]=adj.surv(surv=out[which(out$year==2011),2:4],days=109,totDays=365)
out[,2:4]=apply(out[,2:4],2,function(x)round(x,2))
all.out=out # store pooled yearly estimates

all.out


# 2d) Seasonal survival estimates for resident, migratory, annual pooled (included individuals with unknown migratory status)

# migratory individuals
mig.data=ps[ps$tactic=="Migratory",]
mig.data$season=factor(mig.data$season)
split.data = split(mig.data, mig.data$season) # split data into a list seasons

(season.index=unlist(lapply(split.data, function(x) sum(x$event)))) # look at sum of events in each season

seasons=levels(mig.data$season) # store years for table below

# remove seasons when no events occur (i.e., Survival = 1.00 )
  if(any(season.index==0)){
    split.data=split.data[which(season.index>0)] # remove seasons when no events occur (i.e., Survival = 1.00 )
    seasons=levels(mig.data$season)[which(season.index>0)] # store seaons for table below
  }

# now use prep.bsurv fuction to create data for all seasons
jags.data.list = lapply(split.data, function(x) prep.bsurv(data=x, type="s"))

# run autojags() function from jagsUI package to get surivorship estimates and save as "outSeasonM" for 3 seasons
outSeasonM=lapply(jags.data.list,function(x) autojags(model.file="Bayesian_surv.txt",data=x,parameters.to.save=params,parallel=TRUE,iter.increment=1000,n.chains=2, n.burnin=1000,n.adapt=1000,max.iter=100000,n.thin=1))

out=data.frame(season=seasons,survival=unlist(lapply(outSeasonM, function(x) tail(x$q50$S0,1))),
           lowerBCI=unlist(lapply(outSeasonM, function(x) tail(x$q2.5$S0,1))),
           upperBCI=unlist(lapply(outSeasonM, function(x) tail(x$q97.5$S0,1))))
out[,2:4]=apply(out[,2:4],2,function(x)round(x,2))
season.mig=out # store pooled estimates

season.mig

# resident individuals
res.data=ps[ps$tactic=="Resident",]
res.data$season=factor(res.data$season)
split.data = split(res.data, res.data$season) # split data into a list seasons

(season.index=unlist(lapply(split.data, function(x) sum(x$event)))) # look at sum of events in each season

seasons=levels(res.data$season) # store years for table below

# remove seasons when no events occur (i.e., Survival = 1.00 )
  if(any(season.index==0)){
    split.data=split.data[which(season.index>0)] # remove seasons when no events occur (i.e., Survival = 1.00 )
    seasons=levels(res.data$season)[which(season.index>0)] # store seaons for table below
  }

# now use prep.bsurv fuction to create data for all seasons
jags.data.list = lapply(split.data, function(x) prep.bsurv(data=x, type="s"))

# run autojags() function from jagsUI package to get surivorship estimates and save as "outSeasonR" for 3 seasons
outSeasonR=lapply(jags.data.list,function(x) autojags(model.file="Bayesian_surv.txt",data=x,parameters.to.save=params,parallel=TRUE,iter.increment=1000,n.chains=2, n.burnin=1000,n.adapt=1000,max.iter=100000,n.thin=1))

out=data.frame(season=seasons,survival=unlist(lapply(outSeasonR, function(x) tail(x$q50$S0,1))),
           lowerBCI=unlist(lapply(outSeasonR, function(x) tail(x$q2.5$S0,1))),
           upperBCI=unlist(lapply(outSeasonR, function(x) tail(x$q97.5$S0,1))))
out[,2:4]=apply(out[,2:4],2,function(x)round(x,2))
season.res=out # store pooled estimates

season.res

# all individuals pooled (including unknown migration status)

split.data = split(ps, ps$season) # split data into a list seasons

(season.index=unlist(lapply(split.data, function(x) sum(x$event)))) # look at sum of events in each season

seasons=levels(ps$season) # store years for table below

# remove seasons when no events occur (i.e., Survival = 1.00 )
  if(any(season.index==0)){
    split.data=split.data[which(season.index>0)] # remove seasons when no events occur (i.e., Survival = 1.00 )
    seasons=levels(ps$season)[which(season.index>0)] # store seaons for table below
  }

# now use prep.bsurv fuction to create data for all seasons
jags.data.list = lapply(split.data, function(x) prep.bsurv(data=x, type="s"))

# run autojags() function from jagsUI package to get surivorship estimates and save as "outSeason" for 3 seasons
outSeason=lapply(jags.data.list,function(x) autojags(model.file="Bayesian_surv.txt",data=x,parameters.to.save=params,parallel=TRUE,iter.increment=1000,n.chains=2, n.burnin=1000,n.adapt=1000,max.iter=100000,n.thin=1))

out=data.frame(season=seasons,survival=unlist(lapply(outSeason, function(x) tail(x$q50$S0,1))),
           lowerBCI=unlist(lapply(outSeason, function(x) tail(x$q2.5$S0,1))),
           upperBCI=unlist(lapply(outSeason, function(x) tail(x$q97.5$S0,1))))
out[,2:4]=apply(out[,2:4],2,function(x)round(x,2))
season.out=out # store seasonal pooled estimates

season.out


# 3) cause-specific mortality CIF summary and estimates

# total number of mortalities observed
sum(pa$event)

# migration tactic of mortalities (raw numbers)
with(pa, table(tactic, event))[,2]

# number of mortalities with unknown tactic
length(which(pa$event==1 & is.na(pa$tactic)==TRUE))

# raw proportion of mortalities attributable to each cause (note that this doesn't account for the risk set)
table(pa[pa$event==1,]$cause)[c(1,3:5)]/sum(pa$event)

# 3a) get annual CIFs for each cause for migrants, residents and pooled

# annual migrant CIFs
csm.dataM=prep.bsurv(data=pa[pa$tactic=="Migratory",], type="a") # use prep.bsurv function to prep data for both survival and CIF estimates (type="a")

# set parameters to track
params.csm <- c("h0","S0","cum.haz","h0.all","CIF") # parameters to track (here,S0=baseline survival function,  cum.haz=baseline hazard function, h0=hazard rate, h0.all=cause-specific hazard rates, CIF=cumulative incidence functions for each cause)

# specify model
sink("Bayesian_csm.txt") # JAGS version
cat("
# begin model loop
 model{
  for(j in 1:n.fail){
	for(i in 1 :n.ind){								
	   dN[i,j] ~ dpois(Idt[i,j]) # Poisson Likelihood
	   Idt[i,j] = Y[i,j] * h0[j]
	} # i	
     mu[j] = r * (t[j + 1] - t[j]) * c # prior mean hazard
	   h0[j] ~ dgamma(mu[j], c) # hazard rate
     S0[j] = exp(-cum.haz[j]) # baseline survivorship function
     cum.haz[j] = sum(h0[1:j]) # cumulative hazard function
    } # j
   # set survivorship to 1 at t=0
      S = c(1,S0) # combine 1 into vector with baseline survivorship
   # now loop over each cause and estimate cause-specific hazards
 for(m in 1:n.causes){
  for(j in 1:n.fail.all[m]){
	for(i in 1:n.ind){								
	   dN.all[i,j,m] ~ dpois(Idt.all[i,j,m]) # Poisson Likelihood
	   Idt.all[i,j,m] = Y.all[i,j,m] * h0.all[j,m]
	} # i	
     mu.all[j,m] = r.all * (t.all[j + 1,m] - t.all[j,m]) * c # prior mean hazard
	   h0.all[j,m] ~ dgamma(mu.all[j,m], c)
    # now calculate CIFs
     haz.all[j,m] = h0.all[j,m] * S[times.all[j,m]]
     CIF[j,m] = sum(haz.all[1:j,m])
    } # j
  } # m
   # priors
    c = 0.001 # degree of confidence in guess of hazard rate
    r = 0.003 # guess at the daily hazard rate
    r.all = 0.0005 # guess at daily cause-specific hazard rate
# end model
}		
    ",fill = TRUE)
sink()


# run autojags() function from jagsUI package to get CIF estimates and save as "out.csmM"
out.csmM=autojags(model.file="Bayesian_csm.txt",data=csm.dataM,parameters.to.save=params.csm,parallel=TRUE, iter.increment=1000,n.chains=2, n.burnin=1000,n.adapt=1000,max.iter=100000,n.thin=1)

# organize output for summary

causes=c("Anthro","Censored","Predation","Unknown","Winterkill")

out1a=out.csmM$q50$CIF
colnames(out1a)=causes[c(1,3:5)]
out1a=round(apply(out1a, 2, function(x) tail(na.omit(x),1)),3)

out2a=out.csmM$q2.5$CIF
colnames(out2a)=causes[c(1,3:5)]
out2a=round(apply(out2a, 2, function(x) tail(na.omit(x),1)),3)

out3a=out.csmM$q97.5$CIF
colnames(out3a)=causes[c(1,3:5)]
out3a=round(apply(out3a, 2, function(x) tail(na.omit(x),1)),3)

cif.mig = data.frame(out1a, out2a, out3a)
colnames(cif.mig)=c("CIF","lowerBCI","upperBCI")

cif.mig

# annual resident CIFs
csm.dataR=prep.bsurv(data=pa[pa$tactic=="Resident",], type="a") # use prep.bsurv function to prep data for both survival and CIF estimates (type="a")

# set parameters to track
params.csm <- c("h0","S0","cum.haz","h0.all","CIF") # parameters to track (here,S0=baseline survival function,  cum.haz=baseline hazard function, h0=hazard rate, h0.all=cause-specific hazard rates, CIF=cumulative incidence functions for each cause)

# run autojags() function from jagsUI package to get CIF estimates and save as "out.csmR"
out.csmR=autojags(model.file="Bayesian_csm.txt",data=csm.dataR,parameters.to.save=params.csm,parallel=TRUE, iter.increment=1000,n.chains=2, n.burnin=1000,n.adapt=1000,max.iter=100000,n.thin=1)

# organize output for summary

causes=c("Anthro","Censored","Predation","Unknown","Winterkill")

out1a=out.csmR$q50$CIF
colnames(out1a)=causes[c(1,3:5)]
out1a=round(apply(out1a, 2, function(x) tail(na.omit(x),1)),3)

out2a=out.csmR$q2.5$CIF
colnames(out2a)=causes[c(1,3:5)]
out2a=round(apply(out2a, 2, function(x) tail(na.omit(x),1)),3)

out3a=out.csmR$q97.5$CIF
colnames(out3a)=causes[c(1,3:5)]
out3a=round(apply(out3a, 2, function(x) tail(na.omit(x),1)),3)

cif.res = data.frame(out1a, out2a, out3a)
colnames(cif.res)=c("CIF","lowerBCI","upperBCI")

cif.res

# annual CIFs for all individuals pooled

csm.data=prep.bsurv(data=pa, type="a") # use prep.bsurv function to prep data for both survival and CIF estimates (type="a")

# set parameters to track
params.csm <- c("h0","S0","cum.haz","h0.all","CIF") # parameters to track (here,S0=baseline survival function,  cum.haz=baseline hazard function, h0=hazard rate, h0.all=cause-specific hazard rates, CIF=cumulative incidence functions for each cause)

# run autojags() function from jagsUI package to get CIF estimates and save as "out.csm"
out.csm=autojags(model.file="Bayesian_csm.txt",data=csm.data,parameters.to.save=params.csm,parallel=TRUE, iter.increment=1000,n.chains=2, n.burnin=1000,n.adapt=1000,max.iter=100000,n.thin=1)

# organize output for summary

causes=c("Anthro","Censored","Predation","Unknown","Winterkill")

out1a=out.csm$q50$CIF
colnames(out1a)=causes[c(1,3:5)]
out1a=round(apply(out1a, 2, function(x) tail(na.omit(x),1)),3)

out2a=out.csm$q2.5$CIF
colnames(out2a)=causes[c(1,3:5)]
out2a=round(apply(out2a, 2, function(x) tail(na.omit(x),1)),3)

out3a=out.csm$q97.5$CIF
colnames(out3a)=causes[c(1,3:5)]
out3a=round(apply(out3a, 2, function(x) tail(na.omit(x),1)),3)

cif.all = data.frame(out1a, out2a, out3a)
colnames(cif.all)=c("CIF","lowerBCI","upperBCI")

cif.all

# 3b) Table 2 results - seasonal CIF estimates for migrants, residents, and pooled individuals

# seaonal CIFs migratory individuals
with(mig.data, table(season,cause)) # look at table of causes by seasons

split.data = split(mig.data, mig.data$season) # again, split data into a list seasons

(season.index=unlist(lapply(split.data, function(x) sum(x$event)))) # look at sum of events in each season


# remove seasons when no events occur (i.e., Survival = 1.00 )
  if(any(season.index==0)){
    split.data=split.data[which(season.index>0)] # remove seasons when no events occur (i.e., Survival = 1.00 )
    seasons=levels(mig.data$season)[which(season.index>0)] # store seaons for table below
  }

# now use prep.bsurv fuction to create data for all seasons
jags.data.list = lapply(split.data, function(x) prep.bsurv(data=x, type="a"))

# run autojags() function from jagsUI package to get CIF estimates and save as "outSeason.csmM" for 3 seasons
outSeason.csmM=lapply(jags.data.list,function(x) autojags(model.file="Bayesian_csm.txt",data=x,parameters.to.save=params.csm,parallel=TRUE,iter.increment=1000,n.chains=2, n.burnin=1000,n.adapt=1000,max.iter=100000,n.thin=1))


# organize output for summer
out1b=outSeason.csmM[[1]]$q50$CIF
colnames(out1b)=causes[c(1,3:4)]
out1b=round(apply(out1b, 2, function(x) tail(na.omit(x),1)),3)

out2b=outSeason.csmM[[1]]$q2.5$CIF
colnames(out2b)=causes[c(1,3:4)]
out2b=round(apply(out2b, 2, function(x) tail(na.omit(x),1)),3)

out3b=outSeason.csmM[[1]]$q97.5$CIF
colnames(out3b)=causes[c(1,3:4)]
out3b=round(apply(out3b, 2, function(x) tail(na.omit(x),1)),3)

cif.summer = data.frame(out1b, out2b, out3b)
colnames(cif.summer)=c("CIF","lowerBCI","upperBCI")

# organize output for winter
out1c=outSeason.csmM[[2]]$q50$CIF
colnames(out1c)=causes[c(1,3,5)]
out1c=round(apply(out1c, 2, function(x) tail(na.omit(x),1)),3)

out2c=outSeason.csmM[[2]]$q2.5$CIF
colnames(out2c)=causes[c(1,3,5)]
out2c=round(apply(out2c, 2, function(x) tail(na.omit(x),1)),3)

out3c=outSeason.csmM[[2]]$q97.5$CIF
colnames(out3c)=causes[c(1,3,5)]
out3c=round(apply(out3c, 2, function(x) tail(na.omit(x),1)),3)

cif.winter = data.frame(out1c, out2c, out3c)
colnames(cif.winter)=c("CIF","lowerBCI","upperBCI")

# Migratory CIF seasonal estimates
cif.winter
cif.summer

# seaonal CIFs resident individuals
with(res.data, table(season,cause)) # look at table of causes by seasons

split.data = split(res.data, res.data$season) # again, split data into a list seasons

(season.index=unlist(lapply(split.data, function(x) sum(x$event)))) # look at sum of events in each season


# remove seasons when no events occur (i.e., Survival = 1.00 )
  if(any(season.index==0)){
    split.data=split.data[which(season.index>0)] # remove seasons when no events occur (i.e., Survival = 1.00 )
    seasons=levels(res.data$season)[which(season.index>0)] # store seaons for table below
  }

# now use prep.bsurv fuction to create data for all seasons
jags.data.list = lapply(split.data, function(x) prep.bsurv(data=x, type="a"))

# run autojags() function from jagsUI package to get CIF estimates and save as "outSeason.csmR" for 3 seasons
outSeason.csmR=lapply(jags.data.list,function(x) autojags(model.file="Bayesian_csm.txt",data=x,parameters.to.save=params.csm,parallel=TRUE,iter.increment=1000,n.chains=2, n.burnin=1000,n.adapt=1000,max.iter=100000,n.thin=1))

# organize output for spring for residents
out1a=data.frame(outSeason.csmR[[1]]$q50$CIF)
colnames(out1a)=causes[c(5)]
out1a=round(apply(out1a, 2, function(x) tail(na.omit(x),1)),3)

out2a=data.frame(outSeason.csmR[[1]]$q2.5$CIF)
colnames(out2a)=causes[c(5)]
out2a=round(apply(out2a, 2, function(x) tail(na.omit(x),1)),3)

out3a=data.frame(outSeason.csmR[[1]]$q97.5$CIF)
colnames(out3a)=causes[c(5)]
out3a=round(apply(out3a, 2, function(x) tail(na.omit(x),1)),3)

# organize output for summer for residents
out1b=outSeason.csmR[[2]]$q50$CIF
colnames(out1b)=causes[c(1,4)]
out1b=round(apply(out1b, 2, function(x) tail(na.omit(x),1)),3)

out2b=outSeason.csmR[[2]]$q2.5$CIF
colnames(out2b)=causes[c(1,4)]
out2b=round(apply(out2b, 2, function(x) tail(na.omit(x),1)),3)

out3b=outSeason.csmR[[2]]$q97.5$CIF
colnames(out3b)=causes[c(1,4)]
out3b=round(apply(out3b, 2, function(x) tail(na.omit(x),1)),3)

# organize output for winter for residents
out1c=outSeason.csmR[[3]]$q50$CIF
colnames(out1c)=causes[c(3:5)]
out1c=round(apply(out1c, 2, function(x) tail(na.omit(x),1)),3)

out2c=outSeason.csmR[[3]]$q2.5$CIF
colnames(out2c)=causes[c(3:5)]
out2c=round(apply(out2c, 2, function(x) tail(na.omit(x),1)),3)

out3c=outSeason.csmR[[3]]$q97.5$CIF
colnames(out3c)=causes[c(3:5)]
out3c=round(apply(out3c, 2, function(x) tail(na.omit(x),1)),3)

cif.res.spring = cbind(out1a, out2a, out3a)
colnames(cif.res.spring)=c("CIF","lowerBCI","upperBCI")

cif.res.summer = data.frame(out1b, out2b, out3b)
colnames(cif.res.summer)=c("CIF","lowerBCI","upperBCI")

cif.res.winter = data.frame(out1c, out2c, out3c)
colnames(cif.res.winter)=c("CIF","lowerBCI","upperBCI")

# Resident CIF seasonal estimates
cif.res.winter
cif.res.spring
cif.res.summer


# seaonal CIFs all individuals pooled
with(ps, table(season,cause)) # look at table of causes by seasons

split.data = split(ps, ps$season) # again, split data into a list seasons

(season.index=unlist(lapply(split.data, function(x) sum(x$event)))) # look at sum of events in each season


# remove seasons when no events occur (i.e., Survival = 1.00 )
  if(any(season.index==0)){
    split.data=split.data[which(season.index>0)] # remove seasons when no events occur (i.e., Survival = 1.00 )
    seasons=levels(ps$season)[which(season.index>0)] # store seaons for table below
  }

# now use prep.bsurv fuction to create data for all seasons
jags.data.list = lapply(split.data, function(x) prep.bsurv(data=x, type="a"))

# run autojags() function from jagsUI package to get CIF estimates and save as "outSeason.csm" for 3 seasons
outSeason.csm=lapply(jags.data.list,function(x) autojags(model.file="Bayesian_csm.txt",data=x,parameters.to.save=params.csm,parallel=TRUE,iter.increment=1000,n.chains=2, n.burnin=1000,n.adapt=1000,max.iter=100000,n.thin=1))

# organize output for spring
out1a=outSeason.csm[[1]]$q50$CIF
colnames(out1a)=causes[c(3:5)]
out1a=round(apply(out1a, 2, function(x) tail(na.omit(x),1)),3)

out2a=outSeason.csm[[1]]$q2.5$CIF
colnames(out2a)=causes[c(3:5)]
out2a=round(apply(out2a, 2, function(x) tail(na.omit(x),1)),3)

out3a=outSeason.csm[[1]]$q97.5$CIF
colnames(out3a)=causes[c(3:5)]
out3a=round(apply(out3a, 2, function(x) tail(na.omit(x),1)),3)

cif.spring = data.frame(out1a, out2a, out3a)
colnames(cif.spring)=c("CIF","lowerBCI","upperBCI")

# organize output for summer
out1b=outSeason.csm[[2]]$q50$CIF
colnames(out1b)=causes[c(1,3:4)]
out1b=round(apply(out1b, 2, function(x) tail(na.omit(x),1)),3)

out2b=outSeason.csm[[2]]$q2.5$CIF
colnames(out2b)=causes[c(1,3:4)]
out2b=round(apply(out2b, 2, function(x) tail(na.omit(x),1)),3)

out3b=outSeason.csm[[2]]$q97.5$CIF
colnames(out3b)=causes[c(1,3:4)]
out3b=round(apply(out3b, 2, function(x) tail(na.omit(x),1)),3)

cif.summer = data.frame(out1b, out2b, out3b)
colnames(cif.summer)=c("CIF","lowerBCI","upperBCI")

# organize output for winter
out1c=outSeason.csm[[3]]$q50$CIF
colnames(out1c)=causes[c(1,3:5)]
out1c=round(apply(out1c, 2, function(x) tail(na.omit(x),1)),3)

out2c=outSeason.csm[[3]]$q2.5$CIF
colnames(out2c)=causes[c(1,3:5)]
out2c=round(apply(out2c, 2, function(x) tail(na.omit(x),1)),3)

out3c=outSeason.csm[[3]]$q97.5$CIF
colnames(out3c)=causes[c(1,3:5)]
out3c=round(apply(out3c, 2, function(x) tail(na.omit(x),1)),3)

cif.winter = data.frame(out1c, out2c, out3c)
colnames(cif.winter)=c("CIF","lowerBCI","upperBCI")

# Pooled CIF seasonal estimates
cif.winter
cif.spring
cif.summer


# 4) Factors influencing mortity risk

# 4a) Table S1 results - get Hazard Ratios from Bayesian and frequentist models (full model fit)

# estimate Hazard ratio coefficients 
ph.data = prep.bsurv(data=pa, type="s")

# prep variables for full model
ph.data$tactic = as.vector(factor2ind(pa$tactic, "Resident"))
ph.data$wsi = pa$wsi 
ph.data$pred.wsi = seq(range(pa$wsi)[1],range(pa$wsi)[2],length.out=120) # makes 120 predictions for wsi
ph.data$n.wsi = length(ph.data$pred.wsi)
ph.data$year = pa$year-min(pa$year)+1
ph.data$nyears=length(unique(pa$year))

# set parameters to track
params.bph <- c("psi","h0","beta","hr","cum.haz","S0","S.res","S.mig","S.res.wsi","S.mig.wsi")

# specify model
sink("Bayesian_ph_reg_psi.txt") # JAGS version
cat("
 # begin model loop
 model{
    for(j in 1:n.fail){
	  for(i in 1 :n.ind){								
	   dN[i,j] ~ dpois(Idt[i,j]) # Poisson Likelihood
	   Idt[i,j] = Y[i,j] * h0[j] * exp(beta[1]*tactic[i]+beta[2]*wsi[i]+beta[3]*tactic[i]*wsi[i])
	} # i	
     mu[j] = r * (t[j + 1] - t[j]) * c # prior mean hazard
	   h0[j] ~ dgamma(mu[j], c) # hazard rate
     S0[j] = exp(-cum.haz[j]) # baseline survivorship function
     cum.haz[j] = sum(h0[1:j]) # cumulative hazard function
     # get predicted surviavl curves for migrants and residents adjusted for mean winter severity
     S.res[j] = S0[j]^exp(beta[1]*0+beta[2]*mean(wsi)+beta[3]*0*mean(wsi))
     S.mig[j] = S0[j]^exp(beta[1]*1+beta[2]*mean(wsi)+beta[3]*1*mean(wsi))
    } # j
   # get predicted survival curves for migrants and residents over the range of observed winter severity indices
    for(k in 1:n.wsi){  
    for(j in 1:n.fail){
     S.res.wsi[j,k] = S0[j]^exp(beta[1]*0+beta[2]*pred.wsi[k]+beta[3]*0*pred.wsi[k])
     S.mig.wsi[j,k] = S0[j]^exp(beta[1]*1+beta[2]*pred.wsi[k]+beta[3]*1*pred.wsi[k])
    } # j
    } # k

   # priors
      c = 0.001 # degree of confidence in guess of hazard rate
      r = 0.003 # guess at the daily hazard rate
   for(i in 1 :n.ind){	
      tactic[i]~dbern(psi) #  probability of being a migratory individual in the sample
   }
    # priors for log risk coefficients
    for(k in 1:3){
      beta[k]~dnorm(0,1e-04) # log hazard ratios
      hr[k]=exp(beta[k]) # hazard ratios
    }
      psi~dbeta(1,1)
# end model
 }	
    ",fill = TRUE)
sink()


# estimate the HR for the interaction of unstandardized winter severity index and migration tactic Bayesian PH model
out.bph.m3 = autojags(model.file="Bayesian_ph_reg_psi.txt",data=ph.data,parameters.to.save=params.bph,parallel=TRUE, iter.increment=1000,n.chains=2,n.burnin=1000,n.adapt=1000,max.iter=100000,n.thin=1,codaOnly=c("cum.haz","S0","S.res","S.mig","S.res.wsi","S.mig.wsi"))

# relevel tactic variable so residents are the reference
pa$tactic = relevel(pa$tactic, ref="Resident")

# fit frequentists model (without individuals having unknown migration tacts)
reg.out = coxph(Surv(enter,exit,event)~tactic*wsi,data=pa)

# get hazard ratios
# frequentist
cbind(exp(reg.out$coefficients),exp(confint(reg.out)))

# Bayesian
rbind(round(quantile(out.bph.m3$sims.list$hr[,1],c(0.5,0.025,0.975)),4),round(quantile(out.bph.m3$sims.list$hr[,2],c(0.5,0.025,0.975)),4),round(quantile(out.bph.m3$sims.list$hr[,3],c(0.5,0.025,0.975)),4))


# 4b) Table 3 results - Gibbs variable selection (used standardized winter severity index)

# set up dataframe to hold Gibbs variable selecton results
GVSout=as.data.frame(matrix(0, ncol=5, nrow=5))
row.names(GVSout)=c("000", "100", "010", "110", "111")
colnames(GVSout)=c("10","100","1000","10000","100000")

ph.data = prep.bsurv(data=pa, type="s")

# prep variables for full model
ph.data$tactic = as.vector(factor2ind(pa$tactic, "Resident"))
ph.data$wsi = std2(pa$wsi) # use std2 function to standardize wsi for model selection

ph.data$pred.wsi = seq(range(std2(pa$wsi))[1],range(std2(pa$wsi))[2],length.out=120) # makes 120 predictions for wsi

# number of wsi to predict
ph.data$n.wsi = length(ph.data$pred.wsi)

# get log hazard coefficients with winter severity index standardized
out.bph.m4 = autojags(model.file="Bayesian_ph_reg_psi.txt",data=ph.data,parameters.to.save=params.bph,parallel=TRUE, iter.increment=1000,n.chains=2,n.burnin=1000,n.adapt=1000,max.iter=100000,n.thin=1,codaOnly=c("cum.haz","S0","S.res","S.mig","S.res.wsi","S.mig.wsi"))

# get pseduopriors from previous full cox model fit
ph.data$mu.beta.ps = out.bph.m4$q50$beta # get log coefficients from standardize full model for pseudo-priors

ph.data$tau.beta.ps = 1/((out.bph.m4$sd$beta)^2)# get log coefficients from standardize full model for psuedo-priors

# set parameters to track
params.bph2 <- c("psi","h0","beta","hr","cum.haz","S0","S.res","S.mig","S.res.wsi","S.mig.wsi","g")
# g is an indicator variable for parameter under selection

# specify model
sink("Bayesian_ph_regGVS_psi.txt") # JAGS version
cat("
# begin model loop
 model{
  for(j in 1:n.fail){
	for(i in 1 :n.ind){								
	   dN[i,j] ~ dpois(Idt[i,j]) # Poisson Likelihood
	   Idt[i,j] = Y[i,j] * h0[j] * exp(g[1]*beta[1]*tactic[i]+g[2]*beta[2]*wsi[i]+g[3]*beta[3]*tactic[i]*wsi[i])
	} # i	
     mu[j] = r * (t[j + 1] - t[j]) * c # prior mean hazard
	   h0[j] ~ dgamma(mu[j], c) # hazard rate
     S0[j] = exp(-cum.haz[j]) # baseline survivorship function
     cum.haz[j] = sum(h0[1:j]) # cumulative hazard function
    # get predicted surviavl curves for migrants and residents adjusted for mean winter severity
     S.res[j] = S0[j]^exp(beta[1]*0+beta[2]*mean(wsi)+beta[3]*0*mean(wsi))
     S.mig[j] = S0[j]^exp(beta[1]*1+beta[2]*mean(wsi)+beta[3]*1*mean(wsi))
    } # j
   # get predicted survival curves for migrants and residents over the range of observed winter severity indices
  for(k in 1:n.wsi){  
  for(j in 1:n.fail){
     S.res.wsi[j,k] = S0[j]^exp(beta[1]*0+beta[2]*pred.wsi[k]+beta[3]*0*pred.wsi[k])
     S.mig.wsi[j,k] = S0[j]^exp(beta[1]*1+beta[2]*pred.wsi[k]+beta[3]*1*pred.wsi[k])
      } # j
      } # k
   # priors
    c = 0.001 # degree of confidence in guess of hazard rate
    r = 0.003 # guess at the daily hazard rate

    # prior for latent unobserved migration tactic
      psi~dbeta(1,1) # probability of migrating
    for(i in 1 :n.ind){	
      tactic[i]~dbern(psi) # 55% chance of migrating 
    } # i
    # coefficient estimates
     for(k in 1:3){
      beta[k]~dnorm(beta.mean.prior[k], beta.tau.prior[k]) # log hazard ratios
      hr[k]=exp(beta[k]) # hazard ratios
     }#k       
     # Gibbs Variables Selection (need to code so that only 5 models are considered)
     g[3] ~ dbern(0.20) # given 5 models (1 model has no explanatory variables and is just the baseline hazard)
     pi = g[3] + 0.5*(1-g[3])
     for(j in 1:2){
     g[j]~dbern(pi) # model coefficient indicator variable
     }
    # priors for log risk coefficients
    for(k in 1:3){
      # GVS priors for risk coefficients
      beta.mean.prior[k] = (1-g[k]) * mu.beta.ps[k]
      beta.tau.prior[k] = (1-g[k]) * tau.beta.ps[k] + g[k] * prec.vague
    } # k
# end model
 }		
    ",fill = TRUE)
sink()

# set up a for loop across all 5 priors
for(i in 1:5){
  
  startTime=Sys.time()
  
  ph.data$prec.vague=1/as.numeric(colnames(GVSout)[i]) # set precision across beta coefficients
  
  # full Bayesian PH model with Gibbs variable selection
  out.bph.gvs = autojags(model.file="Bayesian_ph_regGVS_psi.txt",data=ph.data,parameters.to.save=params.bph2,parallel=TRUE,iter.increment=2500,n.chains=2,n.burnin=1000,n.adapt=2500,max.iter=100000,n.thin=1,codaOnly=c("cum.haz","S0","S.res","S.mig","S.res.wsi","S.mig.wsi"))
  
  # save results at each loop
 # saveRDS(out.bph.gvs, paste0("GVS_var",colnames(GVSout)[i],"_bayesph.rds"))
  
  # calculate posterior model probabilities
  g = out.bph.gvs$sims.list$g
  g = paste(g[,1],g[,2],g[,3],sep="")
  
  # fill in GVSout table with estimates
  GVSout[match(names(round(table(g)/length(g),3)),row.names(GVSout)),i]=as.vector(round(table(g)/length(g), 6))
  
  endTime=Sys.time()
  
  cat("Time elapsed: ",paste0(format(endTime-startTime),","),"model run",i,"finished")
  
}

# Examine model selection results
GVSout


# 4c) Goodness-of-fit test using Bayesian hazard model without any variables

gof.data = prep.bsurv(data=pa, type="s") # get just the survival data (not cause-specific mortality)

params.gof = c("h0","Tobs","Tnew","d.new","d.obs","d.exp")

# specify model
sink("Bayesian_ph_reg_postCheck.txt") # JAGS version
cat("
# begin model loop
 model{
  for(j in 1:n.fail){
	for(i in 1 :n.ind){								
	   dN[i,j] ~ dpois(Idt[i,j]) # Poisson Likelihood
	   Idt[i,j] = Y[i,j] * h0[j] 
     dN.new[i,j] ~ dpois(Idt[i,j]) 
	} # i	
     mu[j] = r * (t[j + 1] - t[j]) * c # prior mean hazard
	   h0[j] ~ dgamma(mu[j], c) 
    } # j
   # priors
    c = 0.001 # degree of confidence in guess of hazard rate
    r = 0.003 # guess at the daily hazard rate
   # goodness-of-fit (collapse counts across individuals)
   for(j in 1:n.fail){   
    d.new[j]=sum(dN.new[1:n.ind,j])
    d.obs[j]=sum(dN[1:n.ind,j])
    d.exp[j]=sum(Idt[1:n.ind,j])    
    err[j]=pow(pow(d.obs[j],0.5) - pow(d.exp[j],0.5),2)
    errnew[j]=pow(pow(d.new[j],0.5) - pow(d.exp[j],0.5),2)
     }
   # calculate discrepancy statistic
    Tobs = sum(err[])
    Tnew = sum(errnew[])
    # end model
     }		
    ",fill = TRUE)
sink()

# fit model without covariates to test goodness-of-fit
out.ph.gof = autojags(model.file="Bayesian_ph_reg_postCheck.txt",data=gof.data,parameters.to.save=params.gof,parallel=TRUE,
                      iter.increment=1000,n.chains=2, n.burnin=1000,n.adapt=1000,max.iter=100000,n.thin=1,codaOnly = c("d.new","d.obs","d.exp"))

# see figure creation at bottom of code


# 4d) Get unstandardized hazard ratio from top model that only includes winter severity index and predict survival under range of observed winter conditions
ph.data = prep.bsurv(data=pa, type="s")

# prep varaible for univariate model
ph.data$wsi=pa$wsi

ph.data$pred.wsi = seq(range(pa$wsi)[1],range(pa$wsi)[2],1)
ph.data$n.wsi = length(ph.data$pred.wsi)

# set parameters to track
params.uni <- c("psi","h0","beta","hr","cum.haz","S0","S.wsi") # hr = hazard ratio; beta=log hazard ratio

# specify model
sink("Bayesian_ph_reg_univar_wsi_pred.txt") # JAGS version
cat("
 # begin model loop
 model{
  for(j in 1:n.fail){
	 for(i in 1 :n.ind){								
	   dN[i,j] ~ dpois(Idt[i,j]) # Poisson Likelihood
	   Idt[i,j] = Y[i,j] * h0[j] * exp(beta*wsi[i])
	 } # i	
     mu[j] = r * (t[j + 1] - t[j]) * c # prior mean hazard
	   h0[j] ~ dgamma(mu[j], c) # hazard rate 
     S0[j] = exp(-cum.haz[j]) # baseline survivorship function
     cum.haz[j] = sum(h0[1:j]) # cumulative hazard function
   } # j
   # priors
     c = 0.001 # degree of confidence in guess of hazard rate
     r = 0.003 # guess at the daily hazard rate
    # get predicted survival curves for migrants and residents over the range of observed winter severity indices
  for(k in 1:n.wsi){  
   for(j in 1:n.fail){
     S.wsi[j,k] = S0[j]^exp(beta*pred.wsi[k])
   } # j
   } # k
   # priors for log risk coefficients
     beta~dnorm(0,1e-04) # log hazard ratios
     hr=exp(beta) # hazard ratios
# end model
}	
    ",fill = TRUE)
sink()

# estimate the HR for winter severity index in univariate Bayesian PH model
# * note that we do not use the model with "psi" in the name
out.bph.m1 <- autojags(model.file="Bayesian_ph_reg_univar_wsi_pred.txt",data=ph.data,parameters.to.save=params.uni,parallel=TRUE,iter.increment=1000,n.chains=2, n.burnin=1000,n.adapt=1000,max.iter=100000,n.thin=1)

# give unstandardized estimate of beta coefficient for winter severity index
round(quantile(out.bph.m1$sims.list$hr,c(0.5,0.025,0.975)),4)

# calculate predicted decline in survival over the least to most severe winter severity conditions observed
out.bph.m1$q50$S.wsi[45,1]-out.bph.m1$q50$S.wsi[45,120]

# predicted survival for lowest observed winter severity
quantile(out.bph.m1$sims.list$S.wsi[,45,1],c(0.5,0.025,0.975))

# predicted survival for most severe observed winter severity
quantile(out.bph.m1$sims.list$S.wsi[,45,120],c(0.5,0.025,0.975))

# 4e) Report estimate of interaction hazard ratio from model with unstandardized winter severity index (out.bph.m3 above)
round(quantile(out.bph.m3$sims.list$hr[,3],c(0.5,0.025,0.975)),4)

# 4f) Get unstandardized hazard ratio from top model includes full model and predict survival under range of observed winter conditions
rbind(round(quantile(out.bph.m3$sims.list$hr[,1],c(0.5,0.025,0.975)),4),round(quantile(out.bph.m3$sims.list$hr[,2],c(0.5,0.025,0.975)),4),round(quantile(out.bph.m3$sims.list$hr[,3],c(0.5,0.025,0.975)),4))

# see figure creation at bottom of code to get predictions from model results out.bph.m3



# 5) report competing risk regression estimates (effect of winter severity stratified by mortality sources) that were based on top model (i.e., only with WSI)

csm.data = prep.bsurv(data=pa, type="a") # process cause-specific mortality data for JAGS

csm.data$wsi = pa$wsi

params.csm2 = c("beta","hr")

# specify model
sink("Bayesian_ph_csm_reg.txt") # JAGS version
cat("
 # begin model loop
 model{
   # now loop over each cause and estimate cause-specific regression coefficients for WSI
 for(m in 1:n.causes){
  for(j in 1:n.fail.all[m]){
   for(i in 1:n.ind){								
	   dN.all[i,j,m] ~ dpois(Idt.all[i,j,m]) # Poisson Likelihood
	   Idt.all[i,j,m] = Y.all[i,j,m] * h0.all[j,m] * exp(beta[m]*wsi[i])
	 } # i	
     mu.all[j,m] = r.all * (t.all[j + 1,m] - t.all[j,m]) * c # prior mean hazard
	   h0.all[j,m] ~ dgamma(mu.all[j,m], c)
   } # j
   } # m
   # priors
     c = 0.001 # degree of confidence in guess of hazard rate
     r.all = 0.001 # guess at daily cause-specific hazard rate
    # priors for log risk coefficients
    for(k in 1:n.causes){
     beta[k]~dnorm(0,1e-04) # log hazard ratios
     hr[k]=exp(beta[k]) # hazard ratios
    }
# end model
}	
    ",fill = TRUE)
sink()

# run Bayesian competing risks regression model using jagsUI
out.reg.csm = autojags(model.file="Bayesian_ph_csm_reg.txt",data=csm.data,parameters.to.save=params.csm2,parallel=TRUE, iter.increment=2000,n.chains=2,n.burnin=1000,n.adapt=10000,max.iter=100000,n.thin=1)

# Anthropogenic
cat("Anthropogenic",round(quantile(out.reg.csm$sims.list$hr[,1],c(0.5,0.025,0.975)),4))

# Predation
cat("Predation",round(quantile(out.reg.csm$sims.list$hr[,2],c(0.5,0.025,0.975)),4))

# Unknown
cat("Unknown",round(quantile(out.reg.csm$sims.list$hr[,3],c(0.5,0.025,0.975)),4))

# Winterkill
cat("Winterkill",round(quantile(out.reg.csm$sims.list$hr[,4],c(0.5,0.025,0.975)),4))

# proportion of pronghorn observations that were resident in first half of study (2004-2007)
with(pa[pa$year<=2007,], table(tactic))/sum(with(pa[pa$year<=2007,], table(tactic)))

# proportion of pronghorn observations that were resident in second half of study (2008-2011)
with(pa[pa$year>2007,], table(tactic))/sum(with(pa[pa$year>2007,], table(tactic)))




# 6) code to make figures

#load Windows font Times New Roman
windowsFonts(Times = windowsFont("Times New Roman"))

## MAKE FIGURE 2

# use survfit function to create time,n.event,n.censor for survivorhsip curve
reg.out5<-survfit(Surv(enter,exit,event)~tactic,conf.type="log-log",data=pa)

temp.surv=numeric(length(reg.out5$surv)) # empty vector to Bayesian hold survival estimates
temp.surv[which(reg.out5$n.event>0)]=c(outTactic[[1]]$q50$S0,outTactic[[2]]$q50$S0)
lower=numeric(length(reg.out5$surv))
upper=numeric(length(reg.out5$surv))
lower[which(reg.out5$n.event>0)]=c(outTactic[[1]]$q2.5$S0,outTactic[[2]]$q2.5$S0)
upper[which(reg.out5$n.event>0)]=c(outTactic[[1]]$q97.5$S0,outTactic[[2]]$q97.5$S0)

for(i in 1:length(temp.surv)){
   temp.surv[i]=ifelse(temp.surv[i]==0,temp.surv[i-1],temp.surv[i])
   upper[i]=ifelse(upper[i]==0,upper[i-1],upper[i])
   lower[i]=ifelse(lower[i]==0,lower[i-1],lower[i])
}

# overwrite survival and confidence intervals from frequentist fit with Bayesian estimates
reg.out5$surv=temp.surv
reg.out5$lower=lower
reg.out5$upper=upper

surv.data=surv_summary(reg.out5, data = pa)
levels(surv.data$strata)=levels(surv.data$tactic)
surv.data=rbind(surv.data[surv.data$tactic=="Migratory",],surv.data[surv.data$tactic=="Resident",])

# use functions that were altered from survminer package in R
p6=ggsurvplot_df(fit=surv.data,color=NULL,palette=c("black","darkgray"),ylim=c(0,1),xlim=c(0,365),legend=c(0.15,0.15),legend.labs=c("Resident","Migratory"),
   font.x=14,font.y=14,font.legend=14,legend.title="",break.x.by=60, ylab="Survivorship", xlab="Days since 11 Nov",
   ggtheme = theme_survminer(base_family="Times",base_size=14),censor.shape=124 # Change ggplot2 theme
  ) + scale_color_manual(breaks=c("Migratory","Resident"),values=c("Migratory"="black","Resident"="darkgray"))
p6 

tiff("Figure2.tif", units="in", width=6, height=5, res=600, compression = "lzw")
p6
dev.off()   


## MAKE FIGURE 3

# first get smoothed cause-specific hazard curves for migrants and residents for each mortality cause (using a coefficient for migration tactic rather than fitting curves separately)

# anthro

d.anthro = pa[is.na(pa$tactic)==FALSE,]
d.anthro$event=ifelse(d.anthro$cause=="Anthro",1,0) # set for anthropogenic mortality
d.anthro.mig=d.anthro[d.anthro$tactic=="Migratory",]
d.anthro.mig$tactic=factor(d.anthro.mig$tactic)
d.anthro.res=d.anthro[d.anthro$tactic=="Resident",]
d.anthro.res$tactic=factor(d.anthro.res$tactic)

# now prep the splines for both tactics
bs.list1=prep.bspline(data=d.anthro.mig,nbin=20,nbasis=5,npred=50,degree=3,penalize=TRUE,diffOrder=2)
bs.list2=prep.bspline(data=d.anthro.res,nbin=20,nbasis=5,npred=50,degree=3,penalize=TRUE,diffOrder=2)

raw.list=list(bs.list1[[2]],bs.list2[[2]]) # store raw data

# create data list for JAGS by aggregating migratory and resident data
jags.list=list(y=c(bs.list1[[1]]$y,bs.list2[[1]]$y),P=c(bs.list1[[1]]$P,bs.list2[[1]]$P),nobs=length(c(bs.list1[[1]]$y,bs.list2[[1]]$y)),
   time = c(rep(1:length(bs.list1[[1]]$y),2)),tactic = c(rep(1,length(bs.list1[[1]]$y)),rep(0,length(bs.list2[[1]]$y))),
   X = bs.list1[[1]]$X,X.pred = bs.list1[[1]]$X.pred,beta0 = bs.list1[[1]]$beta0,
   Q = bs.list1[[1]]$Q,m = bs.list1[[1]]$m,K = bs.list1[[1]]$K)

# save parameters to track
bs.params2=c("mu.rep.mig","mu.rep.res","mu","gamma","lambda","beta1","beta00")
# here beta1 is the offset for the migration tactic (Migratory=1,0 otherwise)

# provide a set of inital values
jags.inits=function()
    list(beta1=0.5,beta00=-16,lambda=100)

# specify model
sink("Bayesian_pspline_wtactic.txt") # JAGS version
cat("
model{
  for(i in 1:nobs){
       y[i] ~ dpois(mu[i]*P[i])
       log(mu[i]) = inprod(X[time[i],],beta[]) + beta1*tactic[i] + log(P[i])
    } # i
  # priors
      beta[1:K] ~ dmnorm(beta0[1:K,1] + beta00, lambda*Q[1:K,1:K])
      gamma = beta - beta00 # derive spline coefficient offsets from the overall mean
      lambda ~ dgamma(0.001, 0.001) # gamma prior on dispersion parameter
      beta00 ~ dnorm(0, 1e-6) # shared prior for mean of coefficients
      beta1 ~ dnorm(0, 1e-6) # prior on effect of migration tactic
  # predict smoothed hazards over finer range of times
  for(j in 1:m){
      # migrants
      mu.rep.mig[j] = exp(inprod(X.pred[j,],gamma[]) + beta00 + beta1*1 + log(mean(P)))
      # residents
      mu.rep.res[j] = exp(inprod(X.pred[j,],gamma[]) + beta00 + beta1*0 + log(mean(P)))
  } # j
# end model
}
    ",fill = TRUE)
sink()

# fit Penalized B-spline mixed effects model in JAGS for both tactics
out.bs.mig = autojags(model.file="Bayesian_pspline_wtactic.txt",data=jags.list,inits=jags.inits,parameters.to.save=bs.params2,parallel=TRUE,iter.increment=5000,n.chains=2,n.burnin=1000,n.adapt=20000,max.iter=100000,n.thin=1, codaOnly=c("gamma","mu","mu.rep.mig","mu.rep.res"))

# store output for migrants in a dataframe
bayes.mig=data.frame(x=bs.list1[[3]],y=apply(out.bs.mig$sims.list$mu.rep.mig,2,median),                
                     lci=apply(out.bs.mig$sims.list$mu.rep.mig,2,function(x) quantile(x,0.025)),uci=apply(out.bs.mig$sims.list$mu.rep.mig,2,function(x) quantile(x,0.975)))

# store output for residents in a dataframe
bayes.res=data.frame(x=bs.list2[[3]],y=apply(out.bs.mig$sims.list$mu.rep.res,2,median),              
                     lci=apply(out.bs.mig$sims.list$mu.rep.res,2,function(x) quantile(x,0.025)),uci=apply(out.bs.mig$sims.list$mu.rep.res,2,function(x) quantile(x,0.975)))

plot.data=cbind(rbind(bayes.mig,bayes.res),group=c(rep("Migratory",length(bayes.mig$x)),rep("Resident",length(bayes.res$x))))

# Store results for Anthro
anthro.res = plot.data[plot.data$group=="Resident",]
anthro.mig = plot.data[plot.data$group=="Migratory",]


# winterkill

d.winter = pa[is.na(pa$tactic)==FALSE,]
d.winter$tactic = relevel(d.winter$tactic, ref="Resident")
d.winter$event=ifelse(d.winter$cause=="Winter",1,0) # set for anthropogenic mortality
d.wint.mig = d.winter[d.winter$tactic=="Migratory",]
d.wint.mig$tactic=factor(d.wint.mig$tactic)
d.wint.res = d.winter[d.winter$tactic=="Resident",]
d.wint.res$tactic=factor(d.wint.res$tactic)

# now prep the splines for both tactics
bs.list1=prep.bspline(data=d.wint.mig,nbin=25,nbasis=5,npred=30,degree=3,penalize=TRUE,diffOrder=2)
bs.list2=prep.bspline(data=d.wint.res,nbin=25,nbasis=5,npred=30,degree=3,penalize=TRUE,diffOrder=2)


raw.list=list(bs.list1[[2]],bs.list2[[2]]) # store raw data

# create data list for JAGS by aggregating migratory and resident data
jags.list=list(y=c(bs.list1[[1]]$y,bs.list2[[1]]$y),P=c(bs.list1[[1]]$P,bs.list2[[1]]$P),nobs=length(c(bs.list1[[1]]$y,bs.list2[[1]]$y)),
   time = c(rep(1:length(bs.list1[[1]]$y),2)),tactic = c(rep(1,length(bs.list1[[1]]$y)),rep(0,length(bs.list2[[1]]$y))),
   X = bs.list1[[1]]$X,X.pred = bs.list1[[1]]$X.pred,beta0 = bs.list1[[1]]$beta0,
   Q = bs.list1[[1]]$Q,m = bs.list1[[1]]$m,K = bs.list1[[1]]$K)

# provide a set of inital values 
jags.inits=function()
    list(beta1=-1,beta00=-28,lambda=0.002)

# fit Penalized B-spline mixed effects model in JAGS for both tactics # inits=jags.inits,
out.bs.mig = autojags(model.file="Bayesian_pspline_wtactic.txt",data=jags.list,inits=jags.inits,parameters.to.save=bs.params2,parallel=TRUE,iter.increment=5000,n.chains=2,n.burnin=1000,n.adapt=20000,max.iter=100000,n.thin=1, codaOnly=c("gamma","mu","mu.rep.mig","mu.rep.res"))

# store output for migrants in a dataframe
bayes.mig=data.frame(x=bs.list1[[3]],y=apply(out.bs.mig$sims.list$mu.rep.mig,2,median),                
                     lci=apply(out.bs.mig$sims.list$mu.rep.mig,2,function(x) quantile(x,0.025)),uci=apply(out.bs.mig$sims.list$mu.rep.mig,2,function(x) quantile(x,0.975)))

# store output for residents in a dataframe
bayes.res=data.frame(x=bs.list2[[3]],y=apply(out.bs.mig$sims.list$mu.rep.res,2,median),              
                     lci=apply(out.bs.mig$sims.list$mu.rep.res,2,function(x) quantile(x,0.025)),uci=apply(out.bs.mig$sims.list$mu.rep.res,2,function(x) quantile(x,0.975)))

plot.data=cbind(rbind(bayes.mig,bayes.res),group=c(rep("Migratory",length(bayes.mig$x)),rep("Resident",length(bayes.res$x))))

# Store results for winterkill
winter.res = plot.data[plot.data$group=="Resident",]
winter.mig = plot.data[plot.data$group=="Migratory",]


# unknown

d.unknown = pa[is.na(pa$tactic)==FALSE,]
d.unknown$event=ifelse(d.unknown$cause=="Unknown",1,0) # set for unknownpogenic mortality
d.unknown.mig=d.unknown[d.unknown$tactic=="Migratory",]
d.unknown.mig$tactic=factor(d.unknown.mig$tactic)
d.unknown.res=d.unknown[d.unknown$tactic=="Resident",]
d.unknown.res$tactic=factor(d.unknown.res$tactic)

# now prep the splines for both tactics
bs.list1=prep.bspline(data=d.unknown.mig,nbin=20,nbasis=5,npred=50,degree=3,penalize=TRUE,diffOrder=2)
bs.list2=prep.bspline(data=d.unknown.res,nbin=20,nbasis=5,npred=50,degree=3,penalize=TRUE,diffOrder=2)

raw.list=list(bs.list1[[2]],bs.list2[[2]]) # store raw data

# create data list for JAGS by aggregating migratory and resident data
jags.list=list(y=c(bs.list1[[1]]$y,bs.list2[[1]]$y),P=c(bs.list1[[1]]$P,bs.list2[[1]]$P),nobs=length(c(bs.list1[[1]]$y,bs.list2[[1]]$y)),
   time = c(rep(1:length(bs.list1[[1]]$y),2)),tactic = c(rep(1,length(bs.list1[[1]]$y)),rep(0,length(bs.list2[[1]]$y))),
   X = bs.list1[[1]]$X,X.pred = bs.list1[[1]]$X.pred,beta0 = bs.list1[[1]]$beta0,
   Q = bs.list1[[1]]$Q,m = bs.list1[[1]]$m,K = bs.list1[[1]]$K)

# provide a set of inital values
jags.inits=function()
    list(beta1=-2,beta00=-15,lambda=100)

# fit Penalized B-spline mixed effects model in JAGS for both tactics
out.bs.mig = autojags(model.file="Bayesian_pspline_wtactic.txt",data=jags.list,inits=jags.inits,parameters.to.save=bs.params2,parallel=TRUE,iter.increment=5000,n.chains=2,n.burnin=1000,n.adapt=20000,max.iter=100000,n.thin=1, codaOnly=c("gamma","mu","mu.rep.mig","mu.rep.res"))

# store output for migrants in a dataframe
bayes.mig=data.frame(x=bs.list1[[3]],y=apply(out.bs.mig$sims.list$mu.rep.mig,2,median),                
                     lci=apply(out.bs.mig$sims.list$mu.rep.mig,2,function(x) quantile(x,0.025)),uci=apply(out.bs.mig$sims.list$mu.rep.mig,2,function(x) quantile(x,0.975)))

# store output for residents in a dataframe
bayes.res=data.frame(x=bs.list2[[3]],y=apply(out.bs.mig$sims.list$mu.rep.res,2,median),              
                     lci=apply(out.bs.mig$sims.list$mu.rep.res,2,function(x) quantile(x,0.025)),uci=apply(out.bs.mig$sims.list$mu.rep.res,2,function(x) quantile(x,0.975)))

plot.data=cbind(rbind(bayes.mig,bayes.res),group=c(rep("Migratory",length(bayes.mig$x)),rep("Resident",length(bayes.res$x))))

# Store results for unknown
unknown.res = plot.data[plot.data$group=="Resident",]
unknown.mig = plot.data[plot.data$group=="Migratory",]


# predation

d.pred = pa[is.na(pa$tactic)==FALSE,]
d.pred$event=ifelse(d.pred$cause=="Predation",1,0) # set for predpogenic mortality
d.pred.mig=d.pred[d.pred$tactic=="Migratory",]
d.pred.mig$tactic=factor(d.pred.mig$tactic)
d.pred.res=d.pred[d.pred$tactic=="Resident",]
d.pred.res$tactic=factor(d.pred.res$tactic)

# now prep the splines for both tactics
bs.list1=prep.bspline(data=d.pred.mig,nbin=20,nbasis=5,npred=50,degree=3,penalize=TRUE,diffOrder=2)
bs.list2=prep.bspline(data=d.pred.res,nbin=20,nbasis=5,npred=50,degree=3,penalize=TRUE,diffOrder=2)

raw.list=list(bs.list1[[2]],bs.list2[[2]]) # store raw data

# create data list for JAGS by aggregating migratory and resident data
jags.list=list(y=c(bs.list1[[1]]$y,bs.list2[[1]]$y),P=c(bs.list1[[1]]$P,bs.list2[[1]]$P),nobs=length(c(bs.list1[[1]]$y,bs.list2[[1]]$y)),
   time = c(rep(1:length(bs.list1[[1]]$y),2)),tactic = c(rep(1,length(bs.list1[[1]]$y)),rep(0,length(bs.list2[[1]]$y))),
   X = bs.list1[[1]]$X,X.pred = bs.list1[[1]]$X.pred,beta0 = bs.list1[[1]]$beta0,
   Q = bs.list1[[1]]$Q,m = bs.list1[[1]]$m,K = bs.list1[[1]]$K)

# provide a set of inital values
jags.inits=function()
    list(beta1=1,beta00=-15,lambda=15)

# fit Penalized B-spline mixed effects model in JAGS for both tactics
out.bs.mig = autojags(model.file="Bayesian_pspline_wtactic.txt",data=jags.list,inits=jags.inits,parameters.to.save=bs.params2,parallel=TRUE,iter.increment=5000,n.chains=2,n.burnin=1000,n.adapt=20000,max.iter=100000,n.thin=1, codaOnly=c("gamma","mu","mu.rep.mig","mu.rep.res"))

# store output for migrants in a dataframe
bayes.mig=data.frame(x=bs.list1[[3]],y=apply(out.bs.mig$sims.list$mu.rep.mig,2,median),                
                     lci=apply(out.bs.mig$sims.list$mu.rep.mig,2,function(x) quantile(x,0.025)),uci=apply(out.bs.mig$sims.list$mu.rep.mig,2,function(x) quantile(x,0.975)))

# store output for residents in a dataframe
bayes.res=data.frame(x=bs.list2[[3]],y=apply(out.bs.mig$sims.list$mu.rep.res,2,median),              
                     lci=apply(out.bs.mig$sims.list$mu.rep.res,2,function(x) quantile(x,0.025)),uci=apply(out.bs.mig$sims.list$mu.rep.res,2,function(x) quantile(x,0.975)))

plot.data=cbind(rbind(bayes.mig,bayes.res),group=c(rep("Migratory",length(bayes.mig$x)),rep("Resident",length(bayes.res$x))))

# Store results for predation
predation.res = plot.data[plot.data$group=="Resident",]
predation.mig = plot.data[plot.data$group=="Migratory",]


#set column names for migrant dataset
colnames(anthro.mig)<-c("time","hazard","lower.ci","upper.ci")
colnames(predation.mig)<-c("time","hazard","lower.ci","upper.ci")
colnames(winter.mig)<-c("time","hazard","lower.ci","upper.ci")
colnames(unknown.mig)<-c("time","hazard","lower.ci","upper.ci")


#set column names for resident dataset
colnames(anthro.res)<-c("time","hazard","lower.ci","upper.ci")
colnames(predation.res)<-c("time","hazard","lower.ci","upper.ci")
colnames(winter.res)<-c("time","hazard","lower.ci","upper.ci")
colnames(unknown.res)<-c("time","hazard","lower.ci","upper.ci")

# now make figure 3 as a 2 x 2 panel
tiff("Figure3.tiff", res=600, compression = "lzw", height=10, width=10, units="in")

layout(matrix(c(1,2,3,4), ncol = 2), widths = c(5,5),
       heights = c(5,5), respect = FALSE)


#winter  for mig vs. res
options(digits=3)

par(mai=c(0.22, 0.82, 0.82, 0.30),family= "Times", xpd=NA)

plot(anthro.res$time,anthro.res$hazard,xlab="", ylab="", ylim=c(0, 0.0035),xlim=c(0, 366), type="n",cex.lab=1.9,cex.axis=1.3, family="Cambria",axes=F)
lines(anthro.mig$time,anthro.mig$hazard, lwd=7, lty=1)
lines(anthro.mig$time,anthro.mig$lower.ci, lty=2,lwd=3)
lines(anthro.mig$time,anthro.mig$upper.ci, lty=2,lwd=3)

lines(anthro.res$time,anthro.res$hazard, lwd=7, lty=1, col="darkgray")
lines(anthro.res$time,anthro.res$lower.ci, lty=2,lwd=3, col="darkgray")
lines(anthro.res$time,anthro.res$upper.ci, lty=2,lwd=3, col="darkgray")

text(280, 0.00325,labels="Anthropogenic", cex=2.1)
axis(side=1, at = c(0,90,180,270,360), labels = TRUE, padj=1,pos=0,cex.axis=2)
axis(side=2, at = c(0,0.0005,0.001,0.0015,0.002,0.0025,0.003,0.0035), labels = c("0.000","","0.001","","0.002","","0.003",""),pos=0,cex.axis=2)
mtext('A', side=3, line=-3, at=c(22.5),cex=1.9)

legend(270, 0.004, c("Migratory","Resident"), ncol=2, col = c(1,"darkgray"),xpd=NA,
       text.col = "black", lty = c(1, 1),lwd=c(6,6),cex=2.1, bty="n", seg.len=1,x.intersp = 0.1)

#predation for mig vs. res
par(mai=c(1.02, 0.82 ,0.22, 0.30),family= "Times")

plot(unknown.res$time,unknown.res$hazard,xlab="", ylab="", ylim=c(0, 0.0035),xlim=c(0, 366), type="n",cex.lab=1.9,cex.axis=1.3, family="Times",axes=F)

lines(unknown.mig$time,unknown.mig$hazard, lwd=7, lty=1)
lines(unknown.mig$time,unknown.mig$lower.ci, lty=2,lwd=3)
lines(unknown.mig$time,unknown.mig$upper.ci, lty=2,lwd=3)

lines(unknown.res$time,unknown.res$hazard, lwd=7, lty=1, col="darkgray")
lines(unknown.res$time,unknown.res$lower.ci, lty=2,lwd=3, col="darkgray")
lines(unknown.res$time,unknown.res$upper.ci, lty=2,lwd=3, col="darkgray")


text(305, 0.00325, labels="Unknown", cex=2.1)
axis(side=1, at = c(0,90,180,270,360), labels = TRUE, padj=1,pos=0,cex.axis=2)
axis(side=2, at = c(0,0.0005,0.001,0.0015,0.002,0.0025,0.003,0.0035), labels = c("0.000","","0.001","","0.002","","0.003",""),pos=0,cex.axis=2)
mtext('C', side=3, line=-3, at=c(22.5),cex=1.9)


par(mai=c(0.22, 0.30, 0.82, 0.82),family= "Times")
plot(predation.res$time,predation.res$hazard,xlab="", ylab="", ylim=c(0, 0.0035),xlim=c(0, 366), type="n",cex.lab=1.9,cex.axis=1.3, family="Cambria",axes=F)

lines(predation.mig$time,predation.mig$hazard, lwd=7, lty=1)
lines(predation.mig$time,predation.mig$lower.ci, lty=2,lwd=3)
lines(predation.mig$time,predation.mig$upper.ci, lty=2,lwd=3)

lines(predation.res$time,predation.res$hazard, lwd=7, lty=1, col="darkgray")
lines(predation.res$time,predation.res$lower.ci, lty=2,lwd=3, col="darkgray")
lines(predation.res$time,predation.res$upper.ci, lty=2,lwd=3, col="darkgray")

text(308, 0.00325, labels="Predation", cex=2.1)
axis(side=1, at = c(0,90,180,270,360), labels = TRUE, padj=1,pos=0,cex.axis=2)
axis(side=2, at = c(0,0.0005,0.001,0.0015,0.002,0.0025,0.003,0.0035), labels = c("0.000","","0.001","","0.002","","0.003",""),pos=0,cex.axis=2)
mtext('B', side=3, line=-3, at=c(22.5),cex=1.9)

par(mai=c(1.02, 0.30, 0.22, 0.82),family= "Times")
plot(winter.res$time,winter.res$hazard,xlab="", ylab="", ylim=c(0, 0.0035),xlim=c(0, 366), type="n",cex.lab=1.9,cex.axis=1.3, family="Cambria",axes=F)

lines(winter.mig$time,winter.mig$hazard, lwd=7, lty=1)
lines(winter.mig$time,winter.mig$lower.ci, lty=2,lwd=3)
lines(winter.mig$time,winter.mig$upper.ci, lty=2,lwd=3)

lines(winter.res$time,winter.res$hazard, lwd=7, lty=1, col="darkgray")
lines(winter.res$time,winter.res$lower.ci, lty=2,lwd=3, col="darkgray")
lines(winter.res$time,winter.res$upper.ci, lty=2,lwd=3, col="darkgray")

text(300, 0.00325, labels="Winterkill", cex=2.1)
axis(side=1, at = c(0,90,180,270,360), labels = TRUE, padj=1,pos=0,cex.axis=2)
axis(side=2, at = c(0,0.0005,0.001,0.0015,0.002,0.0025,0.003,0.0035), labels = c("0.000","","0.001","","0.002","","0.003",""),pos=0,cex.axis=2)
mtext('D', side=3, line=-3, at=c(22.5),cex=1.9)

mtext("          Mortality rate/day", side = 2, cex=2.1, outer=T, adj=0.5, line=-2.25)	  
mtext("Days since 11 Nov", side = 1, cex=2.1, outer=T, adj=0.5, line=-2)	

dev.off()


## MAKE FIGURE 4

# uses full model fit from: out.ph.m3
time.cut=dim(out.bph.m3$sims.list$S.res.wsi)[2]

datasurv = data.frame(WSI=rep(seq(range(pa$wsi)[1],range(pa$wsi)[2],length.out=120),2),
                Survivorship=c(apply(out.bph.m3$sims.list$S.mig.wsi, c(2,3), median)[time.cut,],
                apply(out.bph.m3$sims.list$S.res.wsi, c(2,3), median)[time.cut,]),
                lower=c(apply(out.bph.m3$sims.list$S.mig.wsi, c(2,3), function(x) quantile(x,0.025))[time.cut,],
                apply(out.bph.m3$sims.list$S.res.wsi, c(2,3), function(x) quantile(x,0.025))[time.cut,]),
                upper=c(apply(out.bph.m3$sims.list$S.mig.wsi, c(2,3), function(x) quantile(x,0.975))[time.cut,],
                apply(out.bph.m3$sims.list$S.res.wsi, c(2,3), function(x) quantile(x,0.975))[time.cut,]),
                Tactic=as.factor(c(rep("Migratory",120),
                                  rep("Resident",120))))

#plot in ggplot2 
p4<-ggplot(data=datasurv, aes(x = WSI, y = Survivorship,color=Tactic,fill=Tactic,ymin=lower,ymax=upper)) +
   geom_ribbon(aes(ymin=lower, ymax=upper),alpha=0.3, colour = NA) + geom_line(size=2) 
p5<-p4+ theme_classic() +scale_y_continuous(limits=c(0,1))+scale_x_continuous(breaks=c(30,60,90,120,150),limits=c(30,150))+xlab("Winter severity index")+ylab("Predicted survivorship")
p6<- p5 + theme(axis.text.y = element_text(size=14, family="Times"),axis.text.x = element_text(size=14, family="Times"),axis.title.x=element_text(size=14, family="Times"),
                axis.title.y=element_text(size=14, family="Times"),legend.title=element_text(size=14, family="Times"),
                legend.text=element_text(size=14, family="Times"),legend.position=c(0.15,0.15))
p7<- p6 + scale_fill_manual(name="",labels=c("Migratory","Resident"),values = c("darkgray","gray"))  + scale_color_manual(name="",labels=c("Migratory","Resident"),values = c("black","darkgray")) 
p7

tiff("Figure4.tif", units="in", width=6, height=5, res=600, compression = "lzw")
p7
dev.off()  


## MAKE FIGURE S1

tab =  readRDS("tab.rds")
tab$relativeSnow =  tab$snow/sum(tab$snow) *100
tab$relativeTemp =  abs( tab$temp)/sum(abs( tab$temp)) *100

tab2=data.frame(Year=rep(tab$year,2),Variable=c(rep("Snow",length(tab$relativeSnow)),rep("Temp",length(tab$relativeSnow))),
               Relative.index=c(tab$relativeSnow, tab$relativeTemp))

x6=ggplot() + geom_bar(mapping = aes(x = tab$year, y = tab$wsi), stat = "identity", fill = "lightgrey") +
   geom_line(data=tab2, aes(x = Year, y = Relative.index, group=Variable, color=Variable),size=1.75) +
  geom_point(data=tab2, aes(x = Year, y = Relative.index, group=Variable, color=Variable),size=3) +
  scale_y_continuous(name = "Winter severity index", 
     sec.axis = sec_axis(~./150, name = "Relative severity index")) + theme_bw() +
   theme(axis.text=element_text(size=14, family="Times"),axis.title=element_text(size=14, family="Times"),
        legend.position=c(0.12, 0.88),legend.text=element_text(size=14, family="Times")) + scale_x_continuous(breaks=seq(2004,2011,1)) +
        scale_colour_manual(name = "",values=c(plasma(100)[1],plasma(100)[100])) + xlab("Year")

x6

tiff("FigS1.tif", units="in", width=6, height=5, res=600, compression = "lzw")
x6
dev.off()  

## MAKE FIGURE S2

# plot goodness of fit results
discData=data.frame(x=out.ph.gof$sims.list$Tobs,y=out.ph.gof$sims.list$Tnew)

# make label
lab=as.expression(parse(text=paste(expression(Bayesian),"~",expression(italic(P)),"~",expression("-value  =="),"~",as.character(round(mean(out.ph.gof$sims.list$Tnew>out.ph.gof$sims.list$Tobs),2)),sep="")))

# make plot in ggplot2
x2=ggplot()  +
  theme_bw() + geom_point(data=discData, aes(x=x,y=y),shape=21) + geom_abline(intercept = 0, slope = 1) +
  scale_y_continuous(breaks=seq(0,30,5),limits=c(0,30))+scale_x_continuous(breaks=seq(0,30,5),limits=c(0,30))+xlab("Discrepancy for observed data")+ylab("Discrepancy for simulated data") +
  theme(text = element_text(size=14),legend.position="none") +
  annotate("text", x=24, y=1, label=lab)

x2

tiff("FigS2.tif", units="in", width=5, height=5, res=600, compression = "lzw")
x2
dev.off()   

