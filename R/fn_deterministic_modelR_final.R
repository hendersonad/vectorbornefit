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

Deterministic_modelR_final<-function(agestructure=NULL,theta, theta_init, locationI, seroposdates, episeason, include.count=T){
  # These values tell how to match states of compartment with data points
  sim.vals <- seq(0,max(time.vals)-min(time.vals),7) + 7 
  time.vals.sim <-    seq(0,max(sim.vals),dt)
  
  if(agestructure==0){
  init1=c(
    s_init=theta_init[["s_init"]],e_init=theta_init[["i1_init"]],i_init=theta_init[["i1_init"]],r_init=theta_init[["r_init"]],c_init=0,
    sm_init=theta_init[["sm_init"]],em_init=theta_init[["em_init"]],im_init=theta_init[["im_init"]])
  
  # Output simulation data
  output <- simulate_deterministic_noage(theta, init1, time.vals.sim)
  
  # Match compartment states at sim.vals time
  S_traj <- output[match(time.vals.sim,output$time),"s_init"]
  X_traj <- output[match(time.vals.sim,output$time),"sm_init"]
  R_traj <- output[match(time.vals.sim,output$time),"r_init"]
  cases1 <- output[match(time.vals.sim,output$time),"c_init"]
  casecount <- cases1-c(0,cases1[1:(length(time.vals.sim)-1)])
  #casecount[casecount<0] <- 0
  
  # Calculate seropositivity at pre-specified dates and corresponding likelihood
  i=1; seroP=NULL; binom.lik=NULL
  for(date in seroposdates){
    if(date <= max(date.vals) & date >= min(date.vals)){
      seroP[i] <- R_traj[date.vals<date+3.5 & date.vals>date-3.5]/theta[["npop"]]
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
  
  if(likelihood == -Inf){likelihood=-1e10}
  if(is.nan(likelihood)){likelihood=-1e10}
  
  # Return results
  return(list(C_trace=casecount,C_trace=casecount,
              S_trace=S_traj,R_trace=R_traj,X_trace=X_traj,
              lik=likelihood))
  }else{
    init1=c(
      s_initC=theta_init[["s_initC"]],e_initC=theta_init[["i1_initC"]],i_initC=theta_init[["i1_initC"]],r_initC=theta_init[["r_initC"]],c_initC=0,
      sm_initC=theta_init[["sm_initC"]],em_initC=theta_init[["em_initC"]],im_initC=theta_init[["im_initC"]],
      s_initA=theta_init[["s_initA"]],e_initA=theta_init[["i1_initA"]],i_initA=theta_init[["i1_initA"]],r_initA=theta_init[["r_initA"]],c_initA=0,
      sm_initA=theta_init[["sm_initA"]],em_initA=theta_init[["em_initA"]],im_initA=theta_init[["im_initA"]])
  
  # Output simulation data
  output <- simulate_deterministic(theta, init1, time.vals.sim)
  
  # Match compartment states at sim.vals time
  S_trajC <- output[match(time.vals.sim,output$time),"s_initC"]
  S_trajA <- output[match(time.vals.sim,output$time),"s_initA"]
  X_trajC <- output[match(time.vals.sim,output$time),"sm_initC"]
  X_trajA <- output[match(time.vals.sim,output$time),"sm_initA"]
  R_trajC <- output[match(time.vals.sim,output$time),"r_initC"]
  R_trajA <- output[match(time.vals.sim,output$time),"r_initA"]
  cases1 <- output[match(time.vals.sim,output$time),"c_initC"]
  cases2 <- output[match(time.vals.sim,output$time),"c_initA"]
  casecountC <- cases1-c(0,cases1[1:(length(time.vals.sim)-1)])
  casecountA <- cases2-c(0,cases2[1:(length(time.vals.sim)-1)])
  casecount <- casecountC + casecountA
  #casecount[casecount<0] <- 0
  
  # Calculate seropositivity at pre-specified dates and corresponding likelihood
  i=1; seroPC=NULL;seroPA=NULL; binom.lik=NULL
  for(date in seroposdates){
    if(date <= max(date.vals)){
      seroPC[i] <- R_trajC[date.vals<date+3.5 & date.vals>date-3.5]/theta[["npopC"]]
      seroPA[i] <- R_trajA[date.vals<date+3.5 & date.vals>date-3.5]/theta[["npopA"]]
      binom.lik[i] <- (dbinom(nLUM[i], size=nPOP[i], prob=seroPC[i], log = T) +
                         dbinom(nLUM[(length(nLUM)/2)+1], size=nPOP[(length(nLUM)/2)+1], prob=seroPA[i], log = T))
      i <- i+1
    }
  }
  
  # Calculate the likelihood 
  if(include.count==T){
    likelihood <- sum(binom.lik) + sum(log(dnbinom(y.vals,
                                                   mu=theta[["rep"]]*(casecount),
                                                   size=1/theta[["repvol"]])))
  }else{
    likelihood <- sum(binom.lik)
  }
  
  if(likelihood == -Inf){likelihood=-1e10}
  if(is.nan(likelihood)){likelihood=-1e10}
  
  # Return results
  return(list(C_trace=casecount,C_traceC=casecountC,C_traceA=casecountA,
              S_traceC=S_trajC,S_traceA=S_trajA,R_traceC=R_trajC,R_traceA=R_trajA,
              X_traceC=X_trajC,X_traceA=X_trajA,
              lik=likelihood))
  }
}

