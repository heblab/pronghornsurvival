#
#   prep.bsurv function documentation
#   author: Daniel R. Eacker, created 11 Nov 2019
#
#   Function to format standard time-to-event data for Bayesian survival modeling in JAGS
#
#   --- Function definitions ---
#   
#   prep.bsurv(data=data,type="s") # usage
#
#   "data" is a time-to-event data set with at least the fields: "enter","exit","event" (possibly "cause")
#
#   "type" is "s" is for only survival rate estimates; type = "a" is for both survival and cause-specific mortality
#
# * note that the field: "cause" is needed to create the cause-specific mortality data for jags when "a" is selected
# * and that any right-censored individuals should have "Censored" for cause; capture-related mortalities should be removed
#


prep.bsurv=function(data=data,type="s"){
	  
    # create survival data regardless of type input
    if(type=="s" | type=="a"){
      
    # test for basic time-to-event fields (enter should be zero when no left-truncating is present)
    if(length(which((names(data) %in% c("enter","exit","event"))==TRUE))!=3){
      stop("missing one or more of enter, exit, or event fields in data")
    } 
   
    # drop unused levels of cause (due to stratifying datasets) when there are no observed motalities
     data=data[data$cause %in% names(which(table(data$cause)>0)),] 
     data$cause=factor(data$cause)
 
    # store objects to create data for jags list
    times=c(sort(unique(data[data$event==1,]$exit)),max(data$exit)) # get exit times for all mortality events (and max observed time)
     # add 1 to max time of the dataset if the last observed exit time is equal to the max time     
    if(any(times[1:(length(times)-1)]==max(times))){
      times[length(times)] =  times[length(times)] + 1
     }
      obs.t = data$exit # get exit times for all morts and censored cases
	    ent.t= data$enter # get enter times for all records
	    eps = 1.0E-10 # numerical offset
	    n.ind= nrow(data) # number of observations
	    fail = data$event # store fates (1=died, 0) 
	    n.fail=length(unique(data[data$event==1,]$exit))
	    dN=matrix(0, ncol=n.fail, nrow=n.ind) # create empty matrix to hold died=1, lived/cencored=0
      Y=matrix(0, ncol=n.fail, nrow=n.ind) # create empty matrix for indicator in risk set=1, not in risk set=0
      csm=matrix(0, ncol=n.fail, nrow=n.ind) # create empty matrix for cause-specific mortality data
      levels(data$cause)[which(levels(data$cause)=="Censored")]=NA  # set "Censored" in cause-specific mortality to NA's
  
  # loop over observations to create data for JAGS
    for(j in 1 :n.fail){
	    for(i in 1 :n.ind){
	  		Y[i,j] <- step(obs.t[i] - times[j] + eps) * step((times[j] - 1) - ent.t[i]  + eps)
        dN[i,j] <- Y[i, j] * step(times[(j + 1)] - obs.t[i] - eps) * fail[i]
        csm[i,j] = ifelse(dN[i,j]==1,data$cause[i],dN[i,j]) # store cause-specific mortality for use below
	    }
    }
      # return jags list
      jags.list=list(t=times,n.ind=nrow(dN),n.fail=ncol(dN),dN=dN,Y=Y)
       
      # end for type: type=="s" | type=="a"
    }
      
    if(type=="a"){ # logical indicator to also create cause-specific mortality data
  
      # test for basic time-to-event fields (enter should be zero when no left-truncating is present)
    if(length(which((names(data) %in% c("enter","exit","event","cause"))==TRUE))!=4){
      stop("missing one or more of enter, exit, event, caues fields in data")
    } 

      data.split = split(data, data$cause) # split data into list by cause
      times.list = lapply(data.split, function(x) c(sort(unique(x[x$event==1,]$exit)),max(data$exit))) # get exit times for all mortality events (and max observed time)
      max.check.list = lapply(times.list, function(x) any(x[1:length(x)-1]==max(x))) # check if the last observed exit time is equal to the max time for any of the data sets 
  
    # add 1 to max time of the dataset if the last observed exit time is equal to the max time     
    if(any(unlist(max.check.list)==TRUE)){
      
      # store index for times check list
      max.check.index = which(unlist(max.check.list)==TRUE)
       
      # if only 1 strata have the issue
      if(length(max.check.index)==1){ 
         times.temp = times.list[[max.check.index]]
         times.temp[length(times.temp)] = times.temp[length(times.temp)] + 1
         times.list[[max.check.index]] = times.temp
      }else
        # if only multiple strata have the issue
        if(length(max.check.index)>1){
         times.temp = times.list[max.check.index]
         times.temp = lapply(times.temp, function(x) c(x[1:(length(x)-1)],(x[length(x)]+1)))
         times.list[max.check.index] = times.temp
    }  
  }  
      
    # create objects for to make survival data for jags data list
      obs.t.list = lapply(data.split, function(x) x$exit) # all observed exit times
      ent.t.list = lapply(data.split, function(x) x$enter) # all observed enter times
      n.fail.all = as.numeric(unlist(lapply(data.split, function(x) length(unique(x[x$event==1,]$exit))))) # number of unique failure times
      n.causes = length(data.split)
      dN.all = array(0, c(n.ind,max(n.fail.all),n.causes)) # create array to hold survival indicator (1=died from cause k, 0=other
      Y.all=array(0, c(n.ind,max(n.fail.all),n.causes))  # create array to hold survival indicator (1=died from cause k, 0=other
      times.all = matrix(NA, nrow=max(n.fail.all),ncol=n.causes) # matrix of unique event times for each cause
      t.all = matrix(NA, nrow=max(n.fail.all)+1,ncol=n.causes)

  # now loop over each mortality data to create data needed for JAGS
  for(k in 1:length(data.split)){
      dN = ifelse(csm==k,1,0) # set csm to 1 if cause loop is for specific cause
      temp.csm = data.split[[k]]
      Y = jags.list$Y # use same risk indicator for each level of dat
    for(j in 1:n.fail.all[k]){
	    for(i in 1:nrow(dN)){
          dN.all[,1:n.fail.all[k],k]=dN[,which(apply(dN, 2, sum)!=0)]
          Y.all[,1:n.fail.all[k],k]=Y[,which(apply(dN, 2, sum)!=0)]
	     }
    } 
      t.all[1:(n.fail.all[k]+1),k]=c(unique(sort(temp.csm$exit)),ifelse(any((unique(sort(temp.csm$exit)))==max(data$exit))==TRUE,max(data$exit)+1,max(data$exit)))
      times.all[1:n.fail.all[k],k]=match(t.all[1:(n.fail.all[k]),k],times)
    }
   # end for type: type=="a"     
}    
    # logical indicator if only survival is to be estimated (most likely case so test first)
    if(type=="s"){
       return(jags.list)
    } else
    if(type=="a"){
       return(list(t=jags.list$t,n.ind=jags.list$n.ind,n.fail=jags.list$n.fail,dN=jags.list$dN,Y=jags.list$Y,
            dN.all=dN.all,Y.all=Y.all,t.all=t.all,times.all=times.all,n.fail.all=n.fail.all,n.causes=n.causes))
    }
  # end function
}    