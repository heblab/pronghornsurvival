  # load all functions

  # load ipak() function	
  source("functions/ipak.r")

  # load factor2ind() function	
  source("functions/factor2ind.r")

  # source cause.survival function (Heisey and Patterson 2006)
  source("functions/cause_survival.r")

  # load function that creates time-to-event data for Bayesian PH model
  source("functions/prep_bsurv.R")

	# load function to prep B-spline smoothing of daily hazard rate
  source("functions/prep_bspline.R")
  
  # load function to create B-spline basis functions
  source("functions/create_bspline.r")
  
  # make penalty matrix with small numeric offset
  makeQ = function(degree, basis, epsilon=1e-3){
    x <- diag(basis)
    E <- diff(x, differences=degree)
    return( t(E) %*% E + x*epsilon)
  } 
  
  # create negate %in% function
	`%!in%` = Negate(`%in%`)
	
  # simple function to standardize variables
  std2=function(x){
  (x - mean(x))/(2*sd(x))
  }
  
  # load step function
  step=function(x){
  y = ifelse(x>=0,1,0)
  return(y)
  }
  
  