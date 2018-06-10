##' MCMC function for a single iteration using MH algorithm
##' 
##' MCMC function for a single iteration using MH algorithm
#' @param agestructure Binary indicator variable if model is age structured (1) or not (0) between children and adults. Defaults to NULL in which case no sampling happens 
#' @param include.count True or False - whether to include count data in likelihood. Defaults to True
##' @export 

MCMCloop <- function(agestructure, include.count=T){
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
      data <- load.data(agestructure, virusTab[iiH], dataTab[iiH], serology.excel, init.conditions.excel)
        list2env(data,globalenv())
      
      if(m==1){ # Don't resample on 1st step - check the zeroes!
        output_H = SampleTheta(thetaAlltab[m,iiH,],theta_initAlltab[m,iiH,],0*cov_matrix_thetaA,0*cov_matrix_theta_init,agestructure,global=0)
      }else{
        output_H = SampleTheta(thetaAlltab[m,iiH,],theta_initAlltab[m,iiH,],cov_matrix_thetaA,cov_matrix_theta_init,agestructure,global=0)
      } 
      thetaA_star=output_H$thetaS
      theta_init_star=output_H$theta_initS
      
      # Run model simulation
      output1 = Deterministic_modelR_final(agestructure,c(theta_star,thetaA_star), theta_init_star, locationI=locationtab[iiH], seroposdates=seroposdates, episeason=episeason, include.count=include.count)
      sim_marg_lik_star=sim_marg_lik_star+output1$lik
      #Store vales
      thetaAllstar[iiH,]=thetaA_star
      theta_initAllstar[iiH,]=theta_init_star
        cTraceStar[iiH,(1:length(output1$C_trace))]=output1$C_trace
        sTraceStar[iiH,(1:length(output1$S_trace))]=output1$S_trace
        rTraceStar[iiH,(1:length(output1$R_trace))]=output1$R_trace
        xTraceStar[iiH,(1:length(output1$X_trace))]=output1$X_trace
      if(agestructure==1){
        sTraceCStar[iiH,(1:length(output1$S_traceC))]=output1$S_traceC
        sTraceAStar[iiH,(1:length(output1$S_traceA))]=output1$S_traceA
        cTraceStar[iiH,(1:length(output1$C_trace))]=output1$C_trace
        rTraceCStar[iiH,(1:length(output1$R_traceC))]=output1$R_traceC
        rTraceAStar[iiH,(1:length(output1$R_traceA))]=output1$R_traceA
        cTraceCStar[iiH,(1:length(output1$C_traceC))]=output1$C_traceC
        cTraceAStar[iiH,(1:length(output1$C_traceA))]=output1$C_traceA
        xTraceCStar[iiH,(1:length(output1$X_traceC))]=output1$X_traceC
        xTraceAStar[iiH,(1:length(output1$X_traceA))]=output1$X_traceA
      }
      
      prior.theta <- ComputePrior(iiH, thetaAlltab[m,,],thetaAllstar)
      prior.star <- prior.theta$prior.star*prior.star
    } # end loop over regions

  # Calculate probability function - MH algorithm
    #include prior density of initial recovered? 
    #p_thetaI_starC = priorRec(nELISA_C[3]*as.numeric(theta_Istar[1,"r_initC"])/as.numeric(thetaAllstar[1,"npopC"]))
    #p_thetaIC = priorRec(nELISA_C[3]*as.numeric(thetaItab[1,"r_initC"])/as.numeric(thetaAllstar[1,"npopC"]))
    
    # Calculate probability of correct children/adult reporting 
    #p_ratio_star = 1 #priorRatio(c_case_ratioStar*c_a_N[2])
    #p_ratio_base = 1 #priorRatio(c_case_ratioTab*c_a_N[2])
    
    # P(theta | theta_star) #symmetric sampling dist so ratio = 1 
    q_theta_given_theta_star = 1
    q_theta_star_given_theta = 1
    
    val = exp((sim_marg_lik_star-sim_liktab[m]))*(prior.star/prior[m])*(q_theta_given_theta_star/q_theta_star_given_theta) 
    if(is.na(val)){
        output_prob=0
    }else{
        output_prob = min(val, 1)}
    
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
} # End funciton

