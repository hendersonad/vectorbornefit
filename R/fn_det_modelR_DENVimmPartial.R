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
# theta=c(theta_star,thetaA_star,theta_denv); theta_init =theta_init_star; locationI=locationtab[iiH];
# theta=c(thetaA_star,theta_denv); theta_init =theta_init_star; locationI=locationtab[iiH];
# theta=c(thetaMed,theta_star,thetaA_star,theta_denv); theta_init =theta_init_star; locationI=locationtab[iiH];

Deterministic_modelR_final_DENVimmmunityPartial <- function(theta, theta_init, locationI, seroposdates, episeason, include.count=T){
    model.start.date <- as.Date(theta[["model_st"]],origin="1970-01-01") 
    # DENV epidemic has fixed start date
    denv.intro <- as.Date("2013-11-01") 
    
    # Set indicator on when DENV epidemic begins in context of Zika timeline
    if(model.start.date>=denv.intro){
      theta[["denv_start"]] <- 0
        if(model.start.date<min(date.vals)){
          theta[["zika_start"]] <- 0
        }else{
          theta[["zika_start"]] <- time.vals[date.vals<model.start.date+3.5 & date.vals>model.start.date-3.5]
        }
    }else{
      theta[["zika_start"]] <- 0
      data <- load.data.multistart(add.nulls = 0, startdate=model.start.date, virusTab[iiH], dataTab[iiH], serology.excel, init.conditions.excel)
        list2env(data,globalenv())
        if(denv.intro<min(date.vals)){
          theta[["denv_start"]] <- 0
        }else{
          theta[["denv_start"]] <- time.vals[date.vals<=denv.intro+3.5 & date.vals>=denv.intro-3.5]
        }
    }
    
    denv_end <- as.Date("2014-07-13")
    theta[["denv_end"]] = time.vals[date.vals<=denv_end+3.5 & date.vals>=denv_end-3.5]
    
    # These values tell how to match states of compartment with data points
    sim.vals <- seq(0,max(time.vals)-min(time.vals),7) + 7 
    time.vals.sim <- seq(0,max(sim.vals),dt)
    
    # set initial conditions
    init1=c(
      #s_init=theta_init[["s_init"]],e_init=theta_init[["i1_init"]],i_init=theta_init[["i1_init"]],r_init=theta_init[["r_init"]],
      s_init=theta_init[["s_init"]],e_init=0,i_init=0,r_init=theta_init[["r_init"]],
      s2_init=theta_init[["s2_init"]],e2_init=theta_init[["e2_init"]],i2_init=theta_init[["i1_2_init"]],r2_init=theta_init[["r2_init"]],c_init=0,
      #sm_init=theta_init[["sm_init"]],em_init=theta_init[["em_init"]],im_init=theta_init[["im_init"]])
      sm_init=1,em_init=0,im_init=0,
      "D3_init"=0,"FP_init"=0)
    
    # Output simulation data
    #init1[["e_init"]]=1;init1[["i_init"]]=0
    #theta[["psi"]] <- 0.001
    #theta[["chi"]] <- 0.5
    #theta[["beta_v_amp"]] <- 0.1
    #theta[["beta_h"]] <- 0.08
    #theta[["chi"]] <- 0.3
    #theta[["psi"]]=3e-5
    
    output <- simulate_deterministic_noage_DENVimm_partial(theta, init1, time.vals.sim)
    
    # Match compartment states at sim.vals time
    S_traj <- output[match(time.vals.sim,output$time),"s_init"]
    X_traj <- output[match(time.vals.sim,output$time),"sm_init"]
    R_traj <- output[match(time.vals.sim,output$time),"r_init"]
    I_traj <- output[match(time.vals.sim,output$time),"i_init"]
    cases1 <- output[match(time.vals.sim,output$time),"c_init"]
    S2 <- output[match(time.vals.sim,output$time),"s2_init"]
    E2 <- output[match(time.vals.sim,output$time),"e2_init"]
    I2 <- output[match(time.vals.sim,output$time),"i2_init"]
    R2 <- output[match(time.vals.sim,output$time),"r2_init"]
    D3 <- output[match(time.vals.sim,output$time),"D3_init"]
    FP <- output[match(time.vals.sim,output$time),"FP_init"]
    casecount <- cases1-c(0,cases1[1:(length(time.vals.sim)-1)])
    casecount[casecount<0] <- 0
      #plot(date.vals[1:length(casecount)],casecount,type='l',col=2)
    
    #plot(FP,type='l')  
    #plot(D3,type='l')  
    #plot(S2, type='l')
    #plot(S_traj, type='l')
    #plot(S_traj+S2, type='l')
    #plot(R_traj, type='l')
    #plot(I_traj, type='l')
    #plot(E2, type='l')
    #plot(I2, type='l')
    #plot(R2, type='l')
      
      #lines(date.vals[1:length(casecount)],casecountD,type='l')
      ##
      #plot(date.vals[1:length(casecount)],casecount,type='l',col=2)
      #par(new=T)
      #plot(date.vals[1:length(casecount)],R_traj/theta[["npop"]],type='l',col=4,yaxt='n',xaxt='n',ylim=c(0,1))
      #axis(side=4)
      
    # Calculate seropositivity at pre-specified dates and corresponding likelihood
    i=1; seroP=NULL; binom.lik=NULL
    sero.years <- format(as.Date(seroposdates, format="%d/%m/%Y"),"%Y")
    sero.y <- substr(sero.years,3,4)
    lum.y <- c("13","15","17")
    for(date in seroposdates){
      if(date < min(date.vals)){ # if seroprevalence date is before Zika timeline { seroprevalence = Luminex data + epislon}
        seroP[i] <- theta[['epsilon']]
        binom.lik[i] <- (dbinom(nLUM[lum.y==sero.y[i]], size=nPOP[lum.y==sero.y[i]], prob=seroP[i], log = T))
      }else{ # else { seroprevalence = model predicted recovered as a proportion of pop + epsilon }
          seroP[i] <-  (min(R_traj[date.vals<date+3.5 & date.vals>date-3.5])/theta[["npop"]]) + 
                        (1 - min(R_traj[date.vals<date+3.5 & date.vals>date-3.5])/theta[["npop"]])*theta[['epsilon']]
          binom.lik[i] <- (dbinom(nLUM[lum.y==sero.y[i]], size=nPOP[lum.y==sero.y[i]], prob=seroP[i], log = T))
          }
        i <- i+1
      }
    
    ln.denv <- length(denv.timeseries)
    ln.full <- length(y.vals)
      likelihood <- sum(binom.lik) + sum(log(dnbinom(y.vals,
                                                      mu=theta[["rep"]]*(casecount),
                                                     size=1/theta[["repvol"]]))) #+
                                    #sum(log(dnbinom(round(denv.timeseries*theta[["iota"]]),
                                    #                  mu=theta[["rep"]]*(casecount[1:ln.denv]),
                                    #                  size=1/theta[["repvol"]])))
    likelihood=max(-1e10, likelihood)
      if(is.null(likelihood)){likelihood=-1e10}
      if(is.nan(likelihood)){likelihood=-1e10}
      if(is.na(likelihood)){likelihood=-1e10}
      if(length(likelihood)==0){likelihood=-1e10}
      if(likelihood == -Inf){likelihood=-1e10}
    
    # Return results
    output1=list(C_trace=casecount,I_trace=I_traj,
                 S_trace=S_traj,R_trace=R_traj,X_trace=X_traj,
                 lik=likelihood, newDates=date.vals)
}  
