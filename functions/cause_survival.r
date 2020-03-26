###FOR CI's log-log transformation_Chourdhury (2002) Non-parametric confidence interval estimation for competing risks analysis: application to contraceptive data
#Run this code first
######CIF ESTIMATES############################################################################################################

###Heisey and Patterson CIF estimator
"cause.survival" = function(table, p)
{
	assign("p", p)
	
	# create tables to hold the results of GKM survival estimates for all events and for cause specific event
	temp.all <- summary(survfit(Surv(enter, exit, event) ~ 1, conf.type="log-log",data = table))
	temp.s <- summary(survfit(Surv(enter, exit, cause == p) ~ 1, conf.type="log-log",data = table))
	
	#combine the two tables so survival of all events can be combined with those of the cause specific events
	s.df <- data.frame(time = temp.s$time, n.event = temp.s$n.event, n.risk = temp.s$n.risk, survival = temp.s$surv)
	all.df <- data.frame(time = temp.all$time, n.event = temp.all$n.event, n.risk = temp.all$n.risk, survival = temp.all$surv)
	all.s.df <- merge(all.df, s.df, by.x = "time", by.y = "time", all.x = T, suffixes = c(".all", ".s"))
	assign("n", all.s.df)
	x <- length(n[,1])
	
	#create temporary placeholders for the calculation of the mortality rate and the cause-specific cumulative incidence func.
	tmp.string <- numeric(x)
	tmp.string2 <- numeric(x)
	t <- 1

	#cycle through the records of the table, including all events to calculate mortality rate and CIF
	while(t <= x) {
		tmp.string[1] <- n$n.event.s[1]/n$n.risk.s[1]
		if (t == 1) tmp.string[t] <- NA else tmp.string[t] <- (n$survival.all[t-1] * n$n.event.s[t])/n$n.risk.s[t]
		if(is.na(tmp.string[t])) tmp.string2[t] <- NA else tmp.string2[t] <- sum(tmp.string[1:t], na.rm = T)
		t = t + 1
	}
	
	MORT <- data.frame(mort.rate = tmp.string)
	CIF2 <- data.frame(CIF = tmp.string2)
	CIF.s.all <- cbind(all.s.df, MORT, CIF2)

# Calculate the variance, standard error and the Confidence Intervals around CIF

SE <- numeric(x)
totvar.t <- numeric(x)

#Reset all temporary variables
t <- 1
j <- 1
Ij <- 0
cumvar.p1 <- 0
cumvar.p2 <- 0
cumvar.p3 <- 0

#loop for the total number of records
while (t <= x) 
{
   It <- CIF.s.all$CIF[t]
   if(is.na(It)) {		
		CIF.s.all$cumvar[t] <- "NA"
		CIF.s.all$StdErr[t] <- "NA"
		CIF.s.all$CI.l[t] <- "NA"
		CIF.s.all$CI.u[t] <- "NA" 
		
		t = t + 1
	}
   else 
	{
		while (j < t) 
		{
		if(is.na(CIF.s.all$CIF[j]))
			Ij <- Ij			
	     else	
			Ij <- CIF.s.all$CIF[j]		
       cumvar.p1 <- cumvar.p1 + (It - Ij)^2 * (CIF.s.all$n.event.all[j]/(CIF.s.all$n.risk.all[j] * (CIF.s.all$n.risk.all[j] - CIF.s.all$n.event.all[j])))

		if(!is.na(CIF.s.all$CIF[j]))
		{
			if(j == 1)
				Sj3 <- 1
			else
			Sj3 <- CIF.s.all$survival.all[j-1]		
			Ijc <- CIF.s.all$CIF[j]
			cumvar.p3 <- cumvar.p3 + (It - Ijc)*(Sj3)*(CIF.s.all$n.event.all[j] / (CIF.s.all$n.risk.all[j])^2)
		}
		j <- j + 1
   		}

		if (t == 1) 
			Sj2 <- 1  
		else
			Sj2 <- CIF.s.all$survival.all[t-1]  
					
  		cumvar.p2 <- (Sj2)^2 * (((CIF.s.all$n.event.all[t])*(CIF.s.all$n.risk.all[t] - CIF.s.all$n.event.all[t]))/(CIF.s.all$n.risk.all[t])^3) + cumvar.p2

		#total all three components of the variance equation to get the final variance,  generate std. err and confidence intervals using log-log transformation
		#Assign all results to the output table
		totvar.t[t] <- cumvar.p1 + cumvar.p2 - (2 * cumvar.p3)
		CIF.s.all$cumvar[t] <- totvar.t[t] 
		SE[t] <- sqrt(totvar.t[t])
		CIF.s.all$StdErr[t] <- SE[t]
		CIF.s.all$CI.l[t] <- CIF.s.all$CIF[t]^(exp((-1.96 * SE[t])/(CIF.s.all$CIF[t]*log(CIF.s.all$CIF[t]))))
		CIF.s.all$CI.u[t] <- CIF.s.all$CIF[t]^(exp((1.96 * SE[t])/(CIF.s.all$CIF[t]*log(CIF.s.all$CIF[t]))))
		 
		t = t + 1
		j <- 1
   }

cumvar.p1 <- 0
cumvar.p3 <- 0
Ij <- 0
It <- 0
}

#Variance calculations end here ----------



	
	return(CIF.s.all)
}
#####End of Heisey and Patterson estimator




