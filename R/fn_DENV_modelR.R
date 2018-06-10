#' Deterministic SEIR-SEI function with no age structure 
#' 
#' This function solves an SEIR-SEI vector borne disease model.
#' @param agestructure Binary indicator variable if model is age structured (1) or not (0) between children and adults. Defaults to NULL in which case no sampling happens 
#' @param theta Vector of parameters for model
#' @param theta_init Vector of initial conditions for model
#' @param locationI Name of location of data
#' @param seroposdates Vector with dates of seroprevalence surveys
#' @param episeason Vector with start and end point of epidemic case data
#' @param include.count True or False - whether to include count data in likelihood. Defaults to True
#' @keywords deterministic 
#' @export
#' @examples 
#' Deterministic_modelR(c(thetatab[1,],thetaAlltab[1,iiH,]), theta_initAlltab[1,iiH,],locationI=locationtab[iiH],seroposdates= c(as.Date("2013-11-25"), as.Date("2015-11-02")), episeason=c(as.Date("2016-02-22"), as.Date("2016-06-06")))

#Deterministic_modelR_final(agestructure=0,c(thetatab[1,],thetaAlltab[1,iiH,]), theta_initAlltab[1,iiH,], locationtab[iiH], seroposdates, episeason)
#theta=c(thetatab[1,],thetaAlltab[1,iiH,])
#theta_init=theta_initAlltab[1,iiH,]
#locationI=locationtab[iiH]

#Deterministic_modelR_final(agestructure=0,c(thetatab[1,],thetaAlltab[1,iiH,]), theta_initAlltab[1,iiH,], locationtab[iiH], seroposdates, episeason)
#theta=c(theta_star,thetaA_star)
#theta_init=theta_init_star
#locationI=locationtab[iiH]

DENV_modelR<-function(agestructure=NULL,theta, theta_init, locationI, seroposdates, episeason, include.count=T){
  # These values tell how to match states of compartment with data points
  sim.vals <- seq(0,max(time.vals)-min(time.vals),7) + 7 
  time.vals.sim <-    seq(0,max(sim.vals),dt)
  
  if(agestructure==0){
    init1=c(
      s_init=theta_init[["s_init"]],i_init=theta_init[["i1_init"]],r_init=theta_init[["r_init"]],c_init=0)
    
    # Output simulation data
    output <- DENV_ode(theta, init1, time.vals.sim)
    
    # Match compartment states at sim.vals time
    S_traj <- output[match(time.vals.sim,output$time),"s_init"]
    R_traj <- output[match(time.vals.sim,output$time),"r_init"]
    cases1 <- output[match(time.vals.sim,output$time),"c_init"]
    casecount <- cases1-c(0,cases1[1:(length(time.vals.sim)-1)])
      if(casecount[1]==0){casecount[1]=1e-10}
    
    # Calculate seropositivity at pre-specified dates and corresponding likelihood
    i=1; seroP=NULL; binom.lik=NULL
    for(date in seroposdates){
      if(date <= max(date.vals) & date >= min(date.vals)){
        seroP[i] <- min(R_traj[date.vals<date+3.5 & date.vals>date-3.5])/theta[["npop"]]
        binom.lik[i] <- (dbinom(nLUM[i], size=nPOP[i], prob=seroP[i], log = T))
        i <- i+1
      }
    } 
    
    if(include.count==T){
      likelihood <- sum(binom.lik) + sum(log(dnbinom(y.vals,mu=theta[["rep"]]*(casecount),size=1/theta[["repvol"]])))
    }else{
      likelihood <- sum(binom.lik)
    }
    
    if(likelihood == -Inf){likelihood=-1e10}
    if(is.nan(likelihood)){likelihood=-1e10}
    
    # Return results ##CHANGED
    output1=list(C_trace=casecount,
                 S_trace=S_traj,R_trace=R_traj,
                 lik=likelihood)
  }
}  