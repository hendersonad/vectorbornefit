#' Deterministic SEIR-SEI function with no age structure 
#' 
#' This function solves an SEIR-SEI vector borne disease model.
#' @param theta Vector of parameters for model
#' @param theta_init Vector of initial conditions for model
#' @param locationI Name of location of data
#' @param seroposdates Vector with dates of seroprevalence surveys
#' @param episeason Vector with start and end point of epidemic case data
#' @param include.count True or False - whether to include count data in likelihood. Defaults to True
#' @keywords deterministic 
#' @export
#' @examples 

Deterministic_modelR_final_DENVimmmunity <- function(theta, theta_init, locationI, seroposdates, episeason, include.count=T){
  # These values tell how to match states of compartment with data points
  sim.vals <- seq(0,max(time.vals)-min(time.vals),7) + 7 
  time.vals.sim <-    seq(0,max(sim.vals),dt)
  
    # DENV epidemic has fixed start date
    denv.intro <- as.Date("2013-11-01") 
    # Set indicator on when DENV epidemic begins in context of Zika timeline
    if(sum(date.vals<denv.intro+3.5 & date.vals>denv.intro-3.5)==0){
      theta[["denv_start"]] <- 0
    }else{
      theta[["denv_start"]] <- time.vals[date.vals<denv.intro+3.5 & date.vals>denv.intro-3.5]
    }
    
    # set initial conditions
    init1=c(
      s_init=theta_init[["s_init"]],e_init=theta_init[["i1_init"]],i_init=theta_init[["i1_init"]],r_init=theta_init[["r_init"]],c_init=0,
      sd_init=theta_init[["sd_init"]],ed_init=theta_init[["ed_init"]],id_init=theta_init[["id_init"]],rd_init=theta_init[["rd_init"]],cd_init=0,
      sm_init=theta_init[["sm_init"]],em_init=theta_init[["em_init"]],im_init=theta_init[["im_init"]])
    
    # Output simulation data
    output <- simulate_deterministic_noage_DENVimm(theta, init1, time.vals.sim)
    
    # Match compartment states at sim.vals time
    S_traj <- output[match(time.vals.sim,output$time),"s_init"]
    X_traj <- output[match(time.vals.sim,output$time),"sm_init"]
    R_traj <- output[match(time.vals.sim,output$time),"r_init"]
    I_traj <- output[match(time.vals.sim,output$time),"i_init"]
    cases1 <- output[match(time.vals.sim,output$time),"c_init"]
    casesD <- output[match(time.vals.sim,output$time),"cd_init"]
    casecount <- cases1-c(0,cases1[1:(length(time.vals.sim)-1)])
    casecountD <- casesD-c(0,casesD[1:(length(time.vals.sim)-1)])
    casecount[casecount<0] <- 0
      #plot(casecountD,type='l')
      #lines(casecount,type='l',col=2)
    
    # Calculate seropositivity at pre-specified dates and corresponding likelihood
    i=1; seroP=NULL; binom.lik=NULL
    sero.years <- format(as.Date(seroposdates, format="%d/%m/%Y"),"%Y")
    sero.y <- substr(sero.years,3,4)
    lum.y <- c("13","15","17")
    for(date in seroposdates){
      if(date < min(date.vals)){ # if seroprevalence date is before Zika timeline { seroprevalence = Luminex data + epislon}
        seroP[i] <- nLUM[lum.y==sero.y[i]]/nPOP[lum.y==sero.y[i]] + theta[['epsilon']]
        binom.lik[i] <- (dbinom(nLUM[lum.y==sero.y[i]], size=nPOP[lum.y==sero.y[i]], prob=seroP[i], log = T))
        }else{ # else { seroprevalence = model predicted recovered as a proportion of pop + epsilon }
          seroP[i] <-  min(R_traj[date.vals<date+3.5 & date.vals>date-3.5])/theta[["npop"]] + theta[['epsilon']]
          binom.lik[i] <- (dbinom(nLUM[lum.y==sero.y[i]], size=nPOP[lum.y==sero.y[i]], prob=seroP[i], log = T))
          }
        i <- i+1
      }
    
      likelihood <- sum(binom.lik) + sum(log(dnbinom(y.vals,
                                                     mu=theta[["rep"]]*(casecount),
                                                     size=1/theta[["repvol"]])))
    likelihood=max(-1e10, likelihood)
      if(is.null(likelihood)){likelihood=-1e10}
      if(is.nan(likelihood)){likelihood=-1e10}
      if(is.na(likelihood)){likelihood=-1e10}
      if(length(likelihood)==0){likelihood=-1e10}
      if(likelihood == -Inf){likelihood=-1e10}
    
    # Return results
    output1=list(C_trace=casecount,CD_trace=casecountD,
                 S_trace=S_traj,R_trace=R_traj,X_trace=X_traj,
                 lik=likelihood)
}  
