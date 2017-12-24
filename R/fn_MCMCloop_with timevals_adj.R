##' MCMC function for a single iteration using MH algorithm
##' 
##' MCMC function for a single iteration using MH algorithm
#' @param agestructure Binary indicator variable if model is age structured (1) or not (0) between children and adults. Defaults to NULL in which case no sampling happens 
#' @param sample.start.point Whether or not to sample start time 
#' @param all.Priors Whether to estimate all priors - only if running multistart!
#' @param include.count True or False - whether to include count data in likelihood. Defaults to True
#' @param startdate original start date in the data series
##' @export 

MCMCloop_withtimevalsadj <- function(agestructure, sample.start.point=T, all.Priors=F, include.count=T, startdate=as.Date("2013-10-28")){
  for (m in 1:MCMC.runs){
    # Scale COV matrices for resampling using error term epsilon
    if(m==1){
      epsilon0 = 0.001
      cov_matrix_theta=epsilon0*cov_matrix_theta0
      cov_matrix_thetaA=epsilon0*cov_matrix_thetaAll
      cov_matrix_theta_init=epsilon0*cov_matrix_theta_initAll
    }else{
      epsilon0 = max(min(0.1,exp(log(epsilon0)+(accept_rate-0.234)*0.999^m)),1e-6) # Stop epsilon getting too big or small
      cov_matrix_theta=epsilon0*cov_matrix_theta0
      cov_matrix_thetaA=epsilon0*cov_matrix_thetaAll
      cov_matrix_theta_init=epsilon0*cov_matrix_theta_initAll
    }
    
    ## Resample global theta every X step
    if(m %% 2==0){
      output_theta = SampleTheta(thetatab[m,],theta_initAlltab[1,1,],cov_matrix_theta,cov_matrix_theta_init,agestructure,global=1)
      theta_star=output_theta$thetaS
    }else{
      theta_star=thetatab[m,]
    }
    
    ## Resample local parameters every step
    prior.star=1
    sim_marg_lik_star=0
    thetaAllstar=0*thetaAlltab[m,,]
    theta_initAllstar=0*theta_initAlltab[m,,]
    cTraceStar=0*c_trace_tab[m,,]
    sTraceStar=0*s_trace_tab[m,,]
    rTraceStar=0*r_trace_tab[m,,]
    xTraceStar=0*x_trace_tab[m,,]
    sTraceCStar=0*c_trace_tabC[m,,]; cTraceAStar=0*c_trace_tabA[m,,]
    sTraceCStar=0*s_trace_tabC[m,,]; sTraceAStar=0*s_trace_tabA[m,,]
    rTraceCStar=0*r_trace_tabC[m,,]; rTraceAStar=0*r_trace_tabA[m,,]
    xTraceCStar=0*x_trace_tabC[m,,]; xTraceAStar=0*x_trace_tabA[m,,]
    #kk=1
    for(kk in itertab){ 
      iiH=kk
      data <- load.data.multistart(agestructure, add.nulls = 0, startdate, virusTab[iiH], dataTab[iiH], serology.excel, init.conditions.excel)
      list2env(data,globalenv())
      
      if(m==1){ # Don't resample on 1st step - check the zeroes!
        output_H = SampleTheta(thetaAlltab[m,iiH,],theta_initAlltab[m,iiH,],0*cov_matrix_thetaA,0*cov_matrix_theta_init,agestructure,global=0)
      }else{
        output_H = SampleTheta(thetaAlltab[m,iiH,],theta_initAlltab[m,iiH,],cov_matrix_thetaA,cov_matrix_theta_init,agestructure,global=0)
        #output_H$thetaS[['t0']]; log(output_H$thetaS[['t0']])
      } 
      thetaA_star=output_H$thetaS
      theta_init_star=output_H$theta_initS
      
      # Adjust time and date series if start point is flexible
      #data <- fn_adjust_start_time(sample.start.point = sample.start.point, t0=thetaA_star[['t0']])
      #  list2env(data,globalenv())
      #fn_adjust_start_time <- function(sample.start.point=F, t0){
      if(sample.start.point==T){
        Sample.StartTime <- (log(thetaA_star[['t0']]))
        #Sample.StartTime <- (log(1.5))
        if(Sample.StartTime>0){
          #Sample.StartTime <- Sample.StartTime*-1
          new.start.time <- startdate + (Sample.StartTime*365)
          data <- load.data.multistart(agestructure, add.nulls=0, new.start.time, virusTab[iiH], dataTab[iiH], serology.excel, init.conditions.excel)
        }else if(Sample.StartTime<0){
          new.start.time <- startdate + (Sample.StartTime*365)
          #message(paste0("sample start <0  = ",  ((thetaA_star[['t0']]))))  
          #message(paste0("sample start <0  = ",  (log(thetaA_star[['t0']]))))  
          #message(paste0("sample start <0  = ", new.start.time))  
          if(log(thetaA_star[['t0']])==-Inf){new.start.time=startdate}
          data <- load.data.multistart(agestructure, add.nulls=0, new.start.time, virusTab[iiH], dataTab[iiH], serology.excel, init.conditions.excel)
        }
      }else{
        #message(paste0("no new start time = ", startdate))  
        data <- load.data.multistart(agestructure,add.nulls=0, startdate, virusTab[iiH], dataTab[iiH], serology.excel, init.conditions.excel)
      }
      list2env(data, globalenv())
      
      # Run model simulation
      output1 = Deterministic_modelR_final(agestructure,c(theta_star,thetaA_star), theta_init_star, locationI=locationtab[iiH], seroposdates=seroposdates, episeason=episeason, include.count=include.count)
      sim_marg_lik_star=sim_marg_lik_star+output1$lik
      
      #Store vales
      thetaAllstar[iiH,]=thetaA_star
      theta_initAllstar[iiH,]=theta_init_star
      
      # choose selection region so results vector only stores from original STARTDATE onwards - so all results are the same length
      #if(sum(round(date.vals) <= startdate)==0){
      if(length(cTraceStar[iiH,]) >= length(output1$C_trace)){
        start.of.output1 <- 1
        length.of.output1 <- length(output1$C_trace)
        extra.vals <- length(cTraceStar[iiH,]) - length(output1$C_trace)
        if(extra.vals<0){extra.vals=0}
        
        cTraceStar[iiH,]=c(output1$C_trace[(start.of.output1:length.of.output1)],rep(output1$C_trace[length.of.output1], extra.vals))
        sTraceStar[iiH,]=c(output1$S_trace[(start.of.output1:length.of.output1)],rep(output1$S_trace[length.of.output1], extra.vals))
        rTraceStar[iiH,]=c(output1$R_trace[(start.of.output1:length.of.output1)],rep(output1$R_trace[length.of.output1], extra.vals))
        xTraceStar[iiH,]=c(output1$X_trace[(start.of.output1:length.of.output1)],rep(output1$X_trace[length.of.output1], extra.vals))
        if(agestructure==1){
          sTraceCStar[iiH,(1:length.of.results.storage)]=output1$S_traceC[(start.of.output1:length.of.output1)]
          sTraceAStar[iiH,(1:length.of.results.storage)]=output1$S_traceA[(start.of.output1:length.of.output1)]
          cTraceStar[iiH,(1:length.of.results.storage)]=output1$C_trace[(start.of.output1:length.of.output1)]
          rTraceCStar[iiH,(1:length.of.results.storage)]=output1$R_traceC[(start.of.output1:length.of.output1)]
          rTraceAStar[iiH,(1:length.of.results.storage)]=output1$R_traceA[(start.of.output1:length.of.output1)]
          cTraceCStar[iiH,(1:length.of.results.storage)]=output1$C_traceC[(start.of.output1:length.of.output1)]
          cTraceAStar[iiH,(1:length.of.results.storage)]=output1$C_traceA[(start.of.output1:length.of.output1)]
          xTraceCStar[iiH,(1:length.of.results.storage)]=output1$X_traceC[(start.of.output1:length.of.output1)]
          xTraceAStar[iiH,(1:length.of.results.storage)]=output1$X_traceA[(start.of.output1:length.of.output1)]
        }
      }else{
        start.of.output1 <- length(output1$C_trace) - length(cTraceStar[iiH,]) + 1  #max(which(date.vals <= startdate))+1 #inclusive so add 1
        length.of.output1 <- length(output1$C_trace)
        cTraceStar[iiH,]=output1$C_trace[(start.of.output1:length.of.output1)]
        sTraceStar[iiH,]=output1$S_trace[(start.of.output1:length.of.output1)]
        rTraceStar[iiH,]=output1$R_trace[(start.of.output1:length.of.output1)]
        xTraceStar[iiH,]=output1$X_trace[(start.of.output1:length.of.output1)]
        if(agestructure==1){
          sTraceCStar[iiH,(1:length.of.results.storage)]=output1$S_traceC[(start.of.output1:length.of.output1)]
          sTraceAStar[iiH,(1:length.of.results.storage)]=output1$S_traceA[(start.of.output1:length.of.output1)]
          cTraceStar[iiH,(1:length.of.results.storage)]=output1$C_trace[(start.of.output1:length.of.output1)]
          rTraceCStar[iiH,(1:length.of.results.storage)]=output1$R_traceC[(start.of.output1:length.of.output1)]
          rTraceAStar[iiH,(1:length.of.results.storage)]=output1$R_traceA[(start.of.output1:length.of.output1)]
          cTraceCStar[iiH,(1:length.of.results.storage)]=output1$C_traceC[(start.of.output1:length.of.output1)]
          cTraceAStar[iiH,(1:length.of.results.storage)]=output1$C_traceA[(start.of.output1:length.of.output1)]
          xTraceCStar[iiH,(1:length.of.results.storage)]=output1$X_traceC[(start.of.output1:length.of.output1)]
          xTraceAStar[iiH,(1:length.of.results.storage)]=output1$X_traceA[(start.of.output1:length.of.output1)]
        }
      }
      
      #message(paste0("T/F=",sum(date.vals <= startdate)==0," ln dataframe=",length(cTraceStar[iiH,])," // ln y.vals=",length(y.vals)," // start out=",start.of.output1, " len= ", length.of.output1, " = ", length.of.output1-start.of.output1))
      
      if(all.Priors==F){
        prior.theta <- ComputePrior(iiH, thetaAlltab[m,,],thetaAllstar)
      }else{
        prior.theta <- ComputePrior_special(iiH, thetaAlltab[m,,],thetaAllstar, "parameters_est")
      }
      prior.star <- prior.theta$prior.star*prior.star
    } # end loop over regions
    
    # Calculate probability function - MH algorithm
    q_theta_given_theta_star = 1
    q_theta_star_given_theta = 1
    
    val = exp((sim_marg_lik_star-sim_liktab[m]))*(prior.star/prior[m])*(q_theta_given_theta_star/q_theta_star_given_theta) 
    
    if(is.na(val)){
      output_prob=0
    }else{
      output_prob = min(val, 1)
    }
    
    # Update parameter values
    if(runif(1,0,1) < output_prob){
      thetatab[m+1,] = theta_star
      thetaAlltab[m+1,,] = thetaAllstar
      theta_initAlltab[m+1,,] = theta_initAllstar
      c_trace_tab[m+1,,]=cTraceStar
      s_trace_tab[m+1,,]=sTraceStar
      r_trace_tab[m+1,,]=rTraceStar
      x_trace_tab[m+1,,]=xTraceStar
      sim_liktab[m+1] = sim_marg_lik_star
      accepttab[m]=1
      prior[m+1] = prior.star
      if(agestructure==1){
        s_trace_tabC[m+1,,]=sTraceCStar; s_trace_tabA[m+1,,]=sTraceAStar
        r_trace_tabC[m+1,,]=rTraceCStar; r_trace_tabA[m+1,,]=rTraceAStar
        c_trace_tabC[m+1,,]=cTraceCStar; c_trace_tabA[m+1,,]=cTraceAStar
        x_trace_tabC[m+1,,]=xTraceCStar; x_trace_tabA[m+1,,]=xTraceAStar
      }
    }else{
      thetatab[m+1,] = thetatab[m,]
      thetaAlltab[m+1,,] = thetaAlltab[m,,]
      theta_initAlltab[m+1,,] = theta_initAlltab[m,,]
      c_trace_tab[m+1,,]=c_trace_tab[m,,]
      s_trace_tab[m+1,,]=s_trace_tab[m,,]
      r_trace_tab[m+1,,]=r_trace_tab[m,,]
      x_trace_tab[m+1,,]=x_trace_tab[m,,]
      sim_liktab[m+1] = sim_liktab[m]
      accepttab[m]=0
      prior[m+1] = prior[m] 
      if(agestructure==1){
        s_trace_tabC[m+1,,]=s_trace_tabC[m,,]; s_trace_tabA[m+1,,]=s_trace_tabA[m,,]
        r_trace_tabC[m+1,,]=r_trace_tabC[m,,]; r_trace_tabA[m+1,,]=r_trace_tabA[m,,]
        c_trace_tabC[m+1,,]=c_trace_tabC[m,,]; c_trace_tabA[m+1,,]=c_trace_tabA[m,,]
        x_trace_tabC[m+1,,]=x_trace_tabC[m,,]; x_trace_tabA[m+1,,]=x_trace_tabA[m,,]
      }
    }
    
    if(m<20){
      accept_rate=0.234
    }else{
      accept_rate=sum(accepttab[1:m])/m
    }
    
    #message(paste0("lik=",sim_liktab[m+1]," start out=",start.of.output1, "len=", length.of.output1, "=", length.of.output1-length.of.output1))
    if(m %% min(MCMC.runs,100) == 0){
      message(paste0(m,accept_rate,sim_liktab[m],as.Date(startdate+log(thetaAlltab[m+1,1,'t0']),origin="1970-01-01"),r_trace_tab[m+1,1,length(r_trace_tab[1,1,])]/thetatab[1,'npop']))
      #  save(sim_liktab,prior,accepttab,s_trace_tab,c_trace_tab,r_trace_tab,x_trace_tab,thetatab,thetaAlltab,theta_initAlltab,file=paste("posterior_outputZ/outputR_",run.name,".RData",sep=""))
    }
    
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
