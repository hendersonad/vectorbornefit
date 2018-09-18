##' MCMC function for fitting Zika model using MH algorithm
##' 
##' MCMC function for fitting Zika model using MH algorithm
#' @param sammple.start.point Whether or not to sample start time 
#' @param startdate Original start date in the data series
##' @export 

AdaptMCMCloop_withDENVimmunity <- function(sample.start.point=T, startdate=as.Date("2013-10-28"), 
                                      adapt.shape.start=MCMC.runs/10){
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
          covmat.empirical <- cov_matrix_thetaA
          theta.mean <- thetaAlltab_current[1,]
          adapting.shape <- 0
        # initialise counter for storing results (m/thining parameter)
        j=1
    }else{
      epsilon0 = max(min(0.1,exp(log(epsilon0)+(accept_rate-0.234)*0.999^m)),1e-6) # Stop epsilon getting too big or small
      
    if (!is.null(adapt.shape.start) && accept_rate*m >= adapt.shape.start) {
            adapting.shape <- m
        scaling.sd <- 2.38/sqrt(length(thetaAlltab_current[1,]))
        cov_matrix_thetaA <- scaling.sd^2 * covmat.empirical
    }
      cov_matrix_theta=epsilon0*cov_matrix_theta0
      cov_matrix_thetaA=epsilon0*cov_matrix_thetaAll
      cov_matrix_theta_init=epsilon0*cov_matrix_theta_initAll
    }
    
    ## Resample global theta every 2nd step
    if(m %% 2==0){
      output_theta = SampleTheta(thetatab_current,theta_initAlltab_current[iiH,],cov_matrix_theta,cov_matrix_theta_init,global=1)
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
    for(kk in itertab){ 
      iiH=kk
      #iiH=1
      data <- load.data.multistart(add.nulls = 0, startdate, virusTab[iiH], dataTab[iiH], serology.excel, init.conditions.excel)
        list2env(data,globalenv())
      
      if(m==1){ # Don't resample on 1st step - check the zeroes!
        output_H = SampleTheta(thetaAlltab_current[iiH,],theta_initAlltab_current[iiH,],0*cov_matrix_thetaA,0*cov_matrix_theta_init,global=0)
      }else{
        output_H = SampleTheta(thetaAlltab_current[iiH,],theta_initAlltab_current[iiH,],cov_matrix_thetaA,cov_matrix_theta_init,global=0)
        #output_H$thetaS[['t0']]; log(output_H$thetaS[['t0']])
      } 
      thetaA_star=output_H$thetaS
      theta_init_star=output_H$theta_initS
      
      # Adjust time and date series if start point is flexible
      if(sample.start.point==T){
        thetaA_star[['t0']] <- max(thetaA_star[['t0']],0)
        Sample.StartTime <- (log(thetaA_star[['t0']]))
        if(Sample.StartTime>0){
          #new.start.time <- startdate + (Sample.StartTime*365)
          thetaA_star[["model_st"]] <- startdate + (Sample.StartTime*365)
          # data <- load.data.multistart(add.nulls=0, new.start.time, virusTab[iiH], dataTab[iiH], serology.excel, init.conditions.excel)
        }else if(Sample.StartTime<=0){
          #new.start.time <- startdate + (Sample.StartTime*365)
          thetaA_star[["model_st"]] <- startdate + (Sample.StartTime*365)
          if(log(thetaA_star[['t0']])==-Inf){thetaA_star[["model_st"]]=startdate}
          #data <- load.data.multistart(add.nulls=0, new.start.time, virusTab[iiH], dataTab[iiH], serology.excel, init.conditions.excel)
        }
      }else{
        data <- load.data.multistart(add.nulls=0, startdate, virusTab[iiH], dataTab[iiH], serology.excel, init.conditions.excel)
        list2env(data, globalenv())
        thetaA_star[["model_st"]] <- startdate
      }
      
      # Run model simulation
      output1 = Deterministic_modelR_final_DENVimmmunity(theta=c(theta_star,thetaA_star,theta_denv), theta_init_star, locationI=locationtab[iiH], seroposdates=seroposdates, episeason=episeason, include.count=include.count)
      sim_marg_lik_star=sim_marg_lik_star + output1$lik 
      
      # drop model_st parameter
      thetaA_star <- thetaA_star[names(thetaA_star)!="model_st"]
      
      #Store vales
      thetaAllstar[iiH,]=thetaA_star
      theta_initAllstar[iiH,]=theta_init_star
  
      # choose selection region so results vector only stores from original STARTDATE onwards - so all results are the same length
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
      }else{
        start.of.output1 <- length(output1$C_trace) - length(cTraceStar[iiH,]) + 1  #max(which(date.vals <= startdate))+1 #inclusive so add 1
        length.of.output1 <- length(output1$C_trace)
          cTraceStar[iiH,]=output1$C_trace[(start.of.output1:length.of.output1)]
          sTraceStar[iiH,]=output1$S_trace[(start.of.output1:length.of.output1)]
          rTraceStar[iiH,]=output1$R_trace[(start.of.output1:length.of.output1)]
          xTraceStar[iiH,]=output1$X_trace[(start.of.output1:length.of.output1)]
      }
      
      # Calculate prior density for current and proposed theta set
      prior.theta <- ComputePrior(iiH, c(thetatab_current,thetaAlltab_current[iiH,]), c(theta_star,thetaA_star), covartheta = cov_matrix_thetaA)
      prior.star <- log(prior.theta$prior.star*prior.star)
      prior.current <- log(prior.theta$prior*prior.current)
    } # end loop over regions
    
    # Calculate probability function - MH algorithm
    if(cov_matrix_thetaA["chi","chi"]>0){
    q_theta_given_theta_star = as.numeric(log(thetaAllstar[iiH,'chi'])) + as.numeric(log(thetaAllstar[iiH,'psi'])) + as.numeric(log(thetaAllstar[iiH,'iota'])) + as.numeric(log(thetaAllstar[iiH,'rep']))
    q_theta_star_given_theta = as.numeric(log(thetaAlltab_current[1, 'chi'])) + as.numeric(log(thetaAlltab_current[1, 'psi'])) + as.numeric(log(thetaAlltab_current[1, 'iota'])) + as.numeric(log(thetaAlltab_current[1, 'rep'])) 
    #q_theta_given_theta_star = 1
    #q_theta_star_given_theta = 1
    }else{
      q_theta_given_theta_star = 1 
      q_theta_star_given_theta = 1
    }
    if(m == 1){
      prior.star = 1 
    }
    
    val = exp((sim_marg_lik_star-sim_liktab_current) + (prior.star - prior.current) + (q_theta_given_theta_star - q_theta_star_given_theta)) 
    
    if(is.na(val)){
      output_prob=0}else if(is.nan(val)){
        output_prob=0}else if(is.null(val)){
          output_prob=0}else if(length(val)==0){
            output_prob=0}else{
              output_prob = min(val, 1)}
    
    # Update parameter values every k step
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
      }
      accept_rate=sum(accepttab[1:j])/j
      j <- j+1
    }
      # Update current values of parameter if{MH_algorithm_val > runif(1,0,1)}
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
    
    # Empirical covariance matrix
    if (adapting.shape >= 0) {
          covmat_thetaAll <- cov_matrix_thetaA[names(thetaA_star), names(thetaA_star)]
          theta.mean2 <- theta.mean[names(thetaA_star)]
          residual <- as.vector(thetaAlltab[j,1,] - theta.mean)
          covmat_thetaAll <- (covmat_thetaAll * (i - 1) + (i - 1)/i * residual %*% t(residual))/i
          theta.mean <- theta.mean + residual/i
          tmp <- (list(covmat_thetaAll = covmat_thetaAll, theta.mean = theta.mean))
                  covmat.empirical <- tmp$covmat_thetaAll
                  theta.mean <- tmp$theta.mean
    }
    
    if(m %% 1000==0){
      print(c(m, accept_rate, sim_liktab_current, thetaAlltab_current[1,'chi'],
            thetaAlltab_current[1,'beta_h'],
            thetaAlltab_current[1,'iota'],
            thetaAlltab_current[1,'epsilon'],
            thetaAlltab_current[1,'rho'],
            thetaAlltab_current[1,'rep'],
            thetaAlltab_current[1,'psi']))
    }  
    
  } # End MCMC loop
  
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
} #end function
