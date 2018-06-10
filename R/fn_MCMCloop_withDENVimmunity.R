##' MCMC function for a single iteration using MH algorithm
##' 
##' MCMC function for a single iteration using MH algorithm
#' @param agestructure Binary indicator variable if model is age structured (1) or not (0) between children and adults. Defaults to NULL in which case no sampling happens 
#' @param sample.start.point Whether or not to sample start time 
#' @param all.Priors Whether to estimate all priors - only if running multistart!
#' @param include.count True or False - whether to include count data in likelihood. Defaults to True
#' @param startdate original start date in the data series
##' @export 

MCMCloop_withDENVimmunity <- function(agestructure, sample.start.point=T, include.count=T, startdate=as.Date("2013-10-28")){
  for (m in 1:MCMC.runs){
    # Scale COV matrices for resampling using error term epsilon
    if(m==1){
      epsilon0 = 0.001
      cov_matrix_theta=epsilon0*cov_matrix_theta0
      cov_matrix_thetaA=epsilon0*cov_matrix_thetaAll
      cov_matrix_theta_init=epsilon0*cov_matrix_theta_initAll
        # set current values for vectors to be updated in MH algroithm 
        thetatab_current = thetatab[m,]
        thetaAlltab_current = thetaAlltab[m,,]
        theta_initAlltab_current = theta_initAlltab[m,,]
        c_trace_tab_current = c_trace_tab[m,,]
        s_trace_tab_current = s_trace_tab[m,,]
        r_trace_tab_current = r_trace_tab[m,,]
        x_trace_tab_current =  x_trace_tab[m,,]
        sim_liktab_current = sim_liktab[m]
        prior_current = prior[m]
        #prior.theta <- ComputePrior(iiH, c(thetatab_current,thetaAlltab_current[iiH,]), c(thetatab_current,thetaAlltab_current[iiH,]))
        #prior_current <- prior.theta$prior.star
        # initialise counter for storing results (m/thining parameter)
        j=1
    }else{
      epsilon0 = max(min(0.1,exp(log(epsilon0)+(accept_rate-0.234)*0.999^m)),1e-6) # Stop epsilon getting too big or small
      cov_matrix_theta=epsilon0*cov_matrix_theta0
      cov_matrix_thetaA=epsilon0*cov_matrix_thetaAll
      cov_matrix_theta_init=epsilon0*cov_matrix_theta_initAll
    }
    
    ## Resample global theta every X step
    if(m %% 2==0){
      output_theta = SampleTheta(thetatab_current,theta_initAlltab_current[iiH,],cov_matrix_theta,cov_matrix_theta_init,agestructure,global=1)
      theta_star=output_theta$thetaS
    }else{
      theta_star=thetatab_current
    }
    
    ## Resample local parameters every step
    prior.star=1
    prior.current=1
    sim_marg_lik_star=0
    thetaAllstar=0*thetaAlltab_current
    theta_initAllstar=0*theta_initAlltab_current
    cTraceStar=0*c_trace_tab_current
    sTraceStar=0*s_trace_tab_current
    rTraceStar=0*r_trace_tab_current
    xTraceStar=0*x_trace_tab_current
    sTraceCStar=0*c_trace_tab_current; cTraceAStar=0*c_trace_tab_current
    sTraceCStar=0*s_trace_tab_current; sTraceAStar=0*s_trace_tab_current
    rTraceCStar=0*r_trace_tab_current; rTraceAStar=0*r_trace_tab_current
    xTraceCStar=0*x_trace_tab_current; xTraceAStar=0*x_trace_tab_current
    #kk=1
    for(kk in itertab){ 
      iiH=kk
      data <- load.data.multistart(agestructure, add.nulls = 0, startdate, virusTab[iiH], dataTab[iiH], serology.excel, init.conditions.excel)
        list2env(data,globalenv())
      
      if(m==1){ # Don't resample on 1st step - check the zeroes!
        output_H = SampleTheta(thetaAlltab_current[iiH,],theta_initAlltab_current[iiH,],0*cov_matrix_thetaA,0*cov_matrix_theta_init,agestructure,global=0)
      }else{
        output_H = SampleTheta(thetaAlltab_current[iiH,],theta_initAlltab_current[iiH,],cov_matrix_thetaA,cov_matrix_theta_init,agestructure,global=0)
        #output_H$thetaS[['t0']]; log(output_H$thetaS[['t0']])
      } 
      thetaA_star=output_H$thetaS
      theta_init_star=output_H$theta_initS
      
      # Adjust time and date series if start point is flexible
      if(sample.start.point==T){
        Sample.StartTime <- (log(thetaA_star[['t0']]))
        if(Sample.StartTime>0){
          new.start.time <- startdate + (Sample.StartTime*365)
          data <- load.data.multistart(agestructure, add.nulls=0, new.start.time, virusTab[iiH], dataTab[iiH], serology.excel, init.conditions.excel)
        }else if(Sample.StartTime<=0){
          new.start.time <- startdate + (Sample.StartTime*365)
          if(log(thetaA_star[['t0']])==-Inf){new.start.time=startdate}
          data <- load.data.multistart(agestructure, add.nulls=0, new.start.time, virusTab[iiH], dataTab[iiH], serology.excel, init.conditions.excel)
        }
      }else{
        data <- load.data.multistart(agestructure,add.nulls=0, startdate, virusTab[iiH], dataTab[iiH], serology.excel, init.conditions.excel)
      }
      list2env(data, globalenv())

      # Run model simulation
      output1 = Deterministic_modelR_final_DENVimmmunity(agestructure,c(theta_star,thetaA_star,theta_denv), theta_init_star, locationI=locationtab[iiH], seroposdates=seroposdates, episeason=episeason, include.count=include.count)
      sim_marg_lik_star=sim_marg_lik_star + output1$lik

      #Store vales
      thetaAllstar[iiH,]=thetaA_star
      theta_initAllstar[iiH,]=theta_init_star
  
      # choose selection region so results vector only stores from original STARTDATE onwards - so all results are the same length
      #if(sum(round(date.vals) <= startdate)==0){
      if(min(date.vals)>start.output.date){
        start.of.output1 <- 1
        ##get end point of time series
        length.of.output1 <- length(output1$C_trace) - max(0,length(output1$C_trace) - length(cTraceStar[iiH,]))
        extra.vals <- max(0,length(cTraceStar[iiH,]) - length(output1$C_trace))
        extra.zero <- rep(0,extra.vals)
          cTraceStar[iiH,]=c(extra.zero,output1$C_trace[(start.of.output1:length.of.output1)])
          sTraceStar[iiH,]=c(extra.zero,output1$S_trace[(start.of.output1:length.of.output1)])
          rTraceStar[iiH,]=c(extra.zero,output1$R_trace[(start.of.output1:length.of.output1)])
          xTraceStar[iiH,]=c(extra.zero,output1$X_trace[(start.of.output1:length.of.output1)])
            if(agestructure==1){
              sTraceCStar[iiH,]=c(extra.zero,output1$S_traceC[(start.of.output1:length.of.output1)])
              sTraceAStar[iiH,]=c(extra.zero,output1$S_traceA[(start.of.output1:length.of.output1)])
              cTraceStar[iiH,]=c(extra.zero,output1$C_trace[(start.of.output1:length.of.output1)])
              rTraceCStar[iiH,]=c(extra.zero,output1$R_traceC[(start.of.output1:length.of.output1)])
              rTraceAStar[iiH,]=c(extra.zero,output1$R_traceA[(start.of.output1:length.of.output1)])
              cTraceCStar[iiH,]=c(extra.zero,output1$C_traceC[(start.of.output1:length.of.output1)])
              cTraceAStar[iiH,]=c(extra.zero,output1$C_traceA[(start.of.output1:length.of.output1)])
              xTraceCStar[iiH,]=c(extra.zero,output1$X_traceC[(start.of.output1:length.of.output1)])
              xTraceAStar[iiH,]=c(extra.zero,output1$X_traceA[(start.of.output1:length.of.output1)])
            }
      }else{
        start.of.output1 <- length(output1$C_trace) - length(cTraceStar[iiH,]) + 1  #max(which(date.vals <= startdate))+1 #inclusive so add 1
        length.of.output1 <- length(output1$C_trace)
          cTraceStar[iiH,]=output1$C_trace[(start.of.output1:length.of.output1)]
          sTraceStar[iiH,]=output1$S_trace[(start.of.output1:length.of.output1)]
          rTraceStar[iiH,]=output1$R_trace[(start.of.output1:length.of.output1)]
          xTraceStar[iiH,]=output1$X_trace[(start.of.output1:length.of.output1)]
            if(agestructure==1){
              sTraceCStar[iiH,]=output1$S_traceC[(start.of.output1:length.of.output1)]
              sTraceAStar[iiH,]=output1$S_traceA[(start.of.output1:length.of.output1)]
              cTraceStar[iiH,]=output1$C_trace[(start.of.output1:length.of.output1)]
              rTraceCStar[iiH,]=output1$R_traceC[(start.of.output1:length.of.output1)]
              rTraceAStar[iiH,]=output1$R_traceA[(start.of.output1:length.of.output1)]
              cTraceCStar[iiH,]=output1$C_traceC[(start.of.output1:length.of.output1)]
              cTraceAStar[iiH,]=output1$C_traceA[(start.of.output1:length.of.output1)]
              xTraceCStar[iiH,]=output1$X_traceC[(start.of.output1:length.of.output1)]
              xTraceAStar[iiH,]=output1$X_traceA[(start.of.output1:length.of.output1)]
            }
      }
      #message(paste0("T/F=",sum(date.vals <= startdate)==0," ln dataframe=",length(cTraceStar[iiH,])," // ln y.vals=",length(y.vals)," // start out=",start.of.output1, " len= ", length.of.output1, " = ", length.of.output1-start.of.output1))
      prior.theta <- ComputePrior(iiH, c(thetatab_current,thetaAlltab_current[iiH,]), c(theta_star,thetaA_star))
      prior.star <- prior.theta$prior.star*prior.star
      prior.current <- prior.theta$prior*prior.current
    } # end loop over regions
    
    # Calculate probability function - MH algorithm
    q_theta_given_theta_star = 1
    q_theta_star_given_theta = 1
    
    if(m == 1){
      prior.star = 1 
    }
    
    val = exp((sim_marg_lik_star-sim_liktab_current))*(prior.star/prior.current)*(q_theta_given_theta_star/q_theta_star_given_theta) 
    if(is.na(val)){
      output_prob=0}else if(is.nan(val)){
        output_prob=0}else if(is.null(val)){
          output_prob=0}else if(length(val)==0){
            output_prob=0}else{
              output_prob = min(val, 1)}
    
    #print(c(thetaAlltab_current[iiH,'epsilon'],output1$lik,prior.star,output_prob))
    
    # Update parameter values
    MH_random_unif <- runif(1,0,1)
    if(m %% thinning.parameter == 0){
      if(MH_random_unif < output_prob){
        thetatab[j+1,] = theta_star
        thetaAlltab[j+1,,] = thetaAllstar
        theta_initAlltab[j+1,,] = theta_initAllstar
        c_trace_tab[j+1,,]=cTraceStar
        s_trace_tab[j+1,,]=sTraceStar
        r_trace_tab[j+1,,]=rTraceStar
        x_trace_tab[j+1,,]=xTraceStar
        sim_liktab[j+1] = sim_marg_lik_star
        accepttab[j]=1
        prior[j+1] = prior.star
        #if(agestructure==1){
        #  s_trace_tabC[m+1,,]=sTraceCStar; s_trace_tabA[m+1,,]=sTraceAStar
        #  r_trace_tabC[m+1,,]=rTraceCStar; r_trace_tabA[m+1,,]=rTraceAStar
        #  c_trace_tabC[m+1,,]=cTraceCStar; c_trace_tabA[m+1,,]=cTraceAStar
        #  x_trace_tabC[m+1,,]=xTraceCStar; x_trace_tabA[m+1,,]=xTraceAStar
        #}
      }else{
        thetatab[j+1,] = thetatab[j,]
        thetaAlltab[j+1,,] = thetaAlltab[j,,]
        theta_initAlltab[j+1,,] = theta_initAlltab[j,,]
        c_trace_tab[j+1,,]=c_trace_tab[j,,]
        s_trace_tab[j+1,,]=s_trace_tab[j,,]
        r_trace_tab[j+1,,]=r_trace_tab[j,,]
        x_trace_tab[j+1,,]=x_trace_tab[j,,]
        sim_liktab[j+1] = sim_liktab[j]
        accepttab[j]=0
        prior[j+1] = prior[j] 
        #if(agestructure==1){
        #  s_trace_tabC[m+1,,]=s_trace_tabC[m,,]; s_trace_tabA[m+1,,]=s_trace_tabA[m,,]
        #  r_trace_tabC[m+1,,]=r_trace_tabC[m,,]; r_trace_tabA[m+1,,]=r_trace_tabA[m,,]
        #  c_trace_tabC[m+1,,]=c_trace_tabC[m,,]; c_trace_tabA[m+1,,]=c_trace_tabA[m,,]
        #  x_trace_tabC[m+1,,]=x_trace_tabC[m,,]; x_trace_tabA[m+1,,]=x_trace_tabA[m,,]
        #}
      }
      accept_rate=sum(accepttab[1:j])/j
      j <- j+1
    }
      if(MH_random_unif < output_prob){
        thetatab_current = theta_star
        thetaAlltab_current = thetaAllstar
        theta_initAlltab_current = theta_initAllstar
        c_trace_tab_current = cTraceStar
        s_trace_tab_current = sTraceStar
        r_trace_tab_current = rTraceStar
        x_trace_tab_current = xTraceStar
        sim_liktab_current = sim_marg_lik_star
        prior_current = prior.star
      }
    
    if(m<20){
      accept_rate=0.234
    }
    #print(c(m, j, val, prior.star,accept_rate))
    
  } # End MCMC loop
  
  
   if(agestructure==1){
    return(list(sim_liktab=sim_liktab,
                prior=prior,
                accepttab=accepttab,
                c_trace_tab=c_trace_tab,
                s_trace_tab=s_trace_tab,
                r_trace_tab=r_trace_tab,
                x_trace_tab=x_trace_tab,
                c_trace_tabC=c_trace_tabC,
                s_trace_tabC=s_trace_tabC,
                r_trace_tabC=r_trace_tabC,
                x_trace_tabC=x_trace_tabC,
                c_trace_tabA=c_trace_tabA,
                s_trace_tabA=s_trace_tabA,
                r_trace_tabA=r_trace_tabA,
                x_trace_tabA=x_trace_tabA,
                thetatab=thetatab,
                thetaAlltab=thetaAlltab,
                theta_initAlltab=theta_initAlltab))
  }else if(agestructure==0){
    return(list(sim_liktab=sim_liktab,
                prior=prior,
                accepttab=accepttab,
                c_trace_tab=c_trace_tab,
                s_trace_tab=s_trace_tab,
                r_trace_tab=r_trace_tab,
                x_trace_tab=x_trace_tab,
                thetatab=thetatab,
                thetaAlltab=thetaAlltab,
                theta_initAlltab=theta_initAlltab))
  }
  
} #end function
