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

Deterministic_modelR_mosqpop<-function(theta, theta_init, locationI, seroposdates, episeason, include.count=T){
  # These values tell how to match states of compartment with data points
  sim.vals <- seq(0,max(time.vals)-min(time.vals),7) + 7 
  time.vals.sim <-    seq(0,max(sim.vals),dt)
  
    init1=c(
      s_init=theta_init[["s_init"]],e_init=theta_init[["i1_init"]],i_init=theta_init[["i1_init"]],r_init=theta_init[["r_init"]],c_init=0,
      sm_init=theta_init[["sm_init"]],em_init=theta_init[["em_init"]],im_init=theta_init[["im_init"]])
    
    # Output simulation data
    output <- simulate_deterministic_mosqPop(theta, init1, time.vals.sim)
    #plot(output$sm_init,type='l',ylim=c(0,max(output$sm_init*1.1)))
    #plot(output$em_init,type='l',ylim=c(0,max(output$em_init*1.1)))
    #plot(output$im_init,type='l',ylim=c(0,max(output$im_init*1.1)))
    #plot(output$sm_init + output$em_init + output$im_init,type='l')
    #plot(time.vals.sim,output$c_init,type='l')
    #plot(time.vals, time.interventions, axes=F, xlab="", ylab="")
    
    # Match compartment states at sim.vals time
    S_traj <- output[match(time.vals.sim,output$time),"s_init"]
    X_traj <- output[match(time.vals.sim,output$time),"sm_init"]
    R_traj <- output[match(time.vals.sim,output$time),"r_init"]
    cases1 <- output[match(time.vals.sim,output$time),"c_init"]
    casecount <- cases1-c(0,cases1[1:(length(time.vals.sim)-1)])
    casecount[casecount<0] <- 0
    #plot(casecount, type='l')
    
    # Calculate seropositivity at pre-specified dates and corresponding likelihood
    i=1; seroP=NULL; binom.lik=NULL
    for(date in seroposdates){
      if(date <= max(date.vals) & date >= min(date.vals)){
        seroP[i] <- min(R_traj[date.vals<date+3.5 & date.vals>date-3.5])/theta[["npop"]]
        binom.lik[i] <- (dbinom(nLUM[i], size=nPOP[i], prob=seroP[i], log = T))
        i <- i+1
      }
    } 
    
    # Calculate the likelihood 
    #likelihood <- sum(binom.lik) + sum(log(dnbinom(y.vals[date.vals>=episeason[1] & date.vals<=episeason[2]],
    #                                                mu=theta[["rep"]]*(casecount[date.vals>=episeason[1] & date.vals<=episeason[2]]),
    #                                                size=1/theta[["repvol"]])))
    if(include.count==T){
      likelihood <- sum(binom.lik) + sum(log(dnbinom(y.vals,
                                                     mu=theta[["rep"]]*(casecount),
                                                     size=1/theta[["repvol"]])))
    }else{
      likelihood <- sum(binom.lik)
    }
    
    likelihood=max(-1e10, likelihood)
    if(is.null(likelihood)){likelihood=-1e10}
    if(is.nan(likelihood)){likelihood=-1e10}
    if(is.na(likelihood)){likelihood=-1e10}
    if(length(likelihood)==0){likelihood=-1e10}
    if(likelihood == -Inf){likelihood=-1e10}
    
    # Return results ##CHANGED
    output1=list(C_trace=casecount,
                 S_trace=S_traj,R_trace=R_traj,X_trace=X_traj,
                 lik=likelihood)
    # Return results
   return(list(C_trace=casecount,
                   S_trace=S_traj,
                   R_trace=R_traj,
                   X_trace=X_traj,
                   lik=likelihood))
}
  
