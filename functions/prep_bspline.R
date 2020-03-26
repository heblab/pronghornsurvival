#
#   prep.bspline function documentation
#   author: Daniel R. Eacker, created 11 Nov 2019
#
#   Function to format standard time-to-event data for smoothed hazard rate estimation in JAGS
#
#   --- Function definitions ---
#   
#   prep.bspline(data=data,nbin=NULL,nbasis=20,npred=1000,degree=3,penalize=TRUE,diffOrder=2) # usage
#
#   "data" is a time-to-event data set with at least the fields: "enter","exit","event" (possibly "cause")
#
#   "type" is "s" is for only survival rate estimates; type = "a" is for both survival and cause-specific mortality
#
#   "nbin" is number of time bins to split survival data into for smoothing
#
#   "nbasis" is the number of basis functions to use
#
#   "npred" is the number of time points to predict the response over
#
#   "degree" is the degree of the polynomial (default is 3)
#
#   "penalize" is a logical argument to use the penalty matrix (Q) or the indentity matrix (I) (default is TRUE)
#
#   "diffOrder" is order of differences to use (higher order results in more smoothing)
#
# * note that degrees is set to 3 for cubic splines in the create.bspline function 
#
# * note that the bins will be equally spaced (unless set to NULL), and knots and Boundary.knots will automatically be determined
# * by the data
#

prep.bspline=function(data=data,nbin=NULL,nbasis=20,npred=1000,degree=3,penalize=TRUE,diffOrder=2){

  # test for basic time-to-event fields (enter should be zero when no left-truncating is present)
  if(length(which((names(data) %in% c("enter","exit","event"))==TRUE))!=3){
    stop("missing one or more of enter, exit, or event fields in data")
  } 
     enter = data$enter # store enter times
     exit = data$exit # store exit times
     
  if (!is.null(nbin)){
      times = seq(min(enter),max(exit),len=nbin+1) # aportions survival times into equal bins
      } else{
  # if nbin is not specified split the data at each time at exit in the data
      times<- sort(unique(exit)) # sort unique exit times (for event=1 and event=0)
    
      }
    
     # save mid points of time bins depending on how nbin is set
     if (!is.null(nbin)){
     mid.times =  times[-length(times)] + diff(times)/2
     }else{
     mid.times = times 
     }

     x = rescale(mid.times, c(0,1)) # rescale times between 0 and one
	 
     if(is.null(nbasis)){
       N = length(times)
       d <- max(4, floor(N/35))
       nbasis <- floor(N/d - 1)
     } 
     
     # create coefficient matrix for bsplines
     X=create.bspline(x, K = nbasis, bdeg=degree, cyclic=FALSE, xl=0, xr=1)
     matplot(x,X,type="l",xlab="Time",ylab="Basis function, Bj(X)",cex.lab=1.5,cex.axis=1.5,lwd=2)
     
    # test for too few of basis functions that cause NA's in X matrix
    if(any(is.na(X))){
       stop("Too few of basis functions...minimum is 4")
     }

     # make n-order difference matrix for penalties
     if(penalize==TRUE){
        Q = makeQ(diffOrder, dim(X)[2]) 
     }else{
        Q = diag(nbasis)
     }
     
     # store vector to share a common prior for beta coefficients
     beta0 = matrix(0, nrow=nbasis,ncol=1)
 
    # create Poisson count data and raw dataset (omit na's from dataset)
    temp.data = data[,c("enter","exit","event")]
    temp.d=survSplit(Surv(enter,exit,event) ~.,temp.data,cut=unique(c(min(enter),times,max(exit))),episode="bin")     # use survSplit from survival package to split data into the time bins
    person.time=unlist(lapply(split(temp.d, temp.d$bin), function(x) sum((x$exit-x$enter)))) # calculate the total person-time at risk
    y=unlist(lapply(split(temp.d, temp.d$bin), function(x) sum((x$event)))) # calculate total number of events in each interval
    raw.hazard=y/person.time  # calculate raw harzard rate
  
    # store raw data
    raw.data=data.frame(time=mid.times,n.event=as.numeric(y),person.time=as.numeric(person.time),raw.hazard=as.numeric(round(raw.hazard,8)))

    # now predict fine time scale for basis function
    x.pred <- seq(0, 1, length.out=npred+1)
    X.pred=create.bspline(x.pred, K = nbasis, bdeg=degree, cyclic=FALSE, xl=min(x.pred), xr=max(x.pred))
    
    # store data to output to JAGS
    bs.list=list(y=raw.data$n.event,P=raw.data$person.time,X=as.matrix(X),X.pred=as.matrix(X.pred),m=length(x.pred),Q=Q,K=nbasis,nbin=length(raw.data$n.event),beta0=beta0)
    
    # display jags data list (bs.list) and raw data
    return(list(bs.list=bs.list,raw.data=raw.data,pred.time=rescale(x.pred, c(min(times),max(times)))))
   }

