##' Simulate function for SIR
##' 
##' Simulate function for SIR - called as part of denv_MCMC
#' @param dt time step
#' @param theta vector of parameter values. Must include npop1, beta, gamma and rep
#' @param theta_init vector of initial valeus for S I R and C
#' @param time.vals Time series to estimate over
##' @export 

Run_simulation<-function(dt, theta, theta_init,time.vals){
  # Define simulation model ODEs
  simulate_deterministic <- function(theta, init.state, times) {
    SIR_ode <- function(time, state, theta) {
      ## define variables
      S <- state[["s_init"]]  
      E <- state[["e_init"]] 
      I <- state[["i1_init"]] 
      R <- state[["r_init"]] 
      C <-state[["c_init"]]   
      N1 <- 342000
      ## extract parameters
      beta <- theta[["beta"]] 
      alpha <- theta[["alpha"]]
      gamma <- theta[["gamma"]]
      repR <- theta[["rep"]]   
      
      foi1 <- beta*(I)/N1 #Force of infection
      
      # Define transition equations
      dS <- -S*foi1         #change in suceptibles
      dE <- S*foi1-alpha*E  #change in pre-infectious
      dI <- alpha*E - gamma*I  #change in infected
      dR <- gamma*I         #change in recovered
      dC <- alpha*E     #prop cases reported*(tranistions to I compartment)
      
      return(list(c(dS,dE,dI,dR,dC)))
    }
    traj <- as.data.frame(ode(init.state, times, SIR_ode, theta, method = "ode45"))
    return(traj)
  }
  
  # Read in initial conditions
  init1 <- c(s_init=theta_init[["s_init"]],
             e_init=theta_init[["e_init"]],
             i1_init=theta_init[["i1_init"]],
             r_init=theta_init[["r_init"]],
             c_init=0)
  
  # Output simulation data
  output <- simulate_deterministic(theta,init1,seq(0,max(time.vals),dt)) 
  S_traj <- output[match(time.vals,output$time),"s_init"]               
  I_traj <- output[match(time.vals,output$time),"i1_init"]              
  
  # Calculate incidence in two populations
  cases1 <- output[match(time.vals,output$time),"c_init"]
  casecount <- cases1-c(0,cases1[1:(length(time.vals)-1)])
  
  # Return selected outputs
  return(list(C1_trace=casecount,S_trace=S_traj,I_trace=I_traj))
}

##' R Code for SEIR DENV epidemic in a single population
##' 
##' Single SEIR epidemic
#' @param y.vals data series
#' @param time.vals time series
#' @param MCMC.runs number of iterations for MCMC loop
##' @export 

DENV_mcmc <- function(y.vals, time.vals, MCMC.runs){
# Set initial conditions
popsize <- 342000
theta <- c(beta=0.27, # tranmission rate
           gamma=1/7, # mean duration of infectiousness
           alpha=1/5.9, # mean incubation period
           rep=0.03, # proportion of cases reported
           repvol=0.1 
)

# Set initial conditions
initial_inf <- 50 # initial number of infected people at start date (Estimated in MCMC loop)
theta_init <- c(s_init=NA,e_init=0,i1_init=initial_inf,r_init=0)
theta_init[["s_init"]] <- popsize-theta_init[["i1_init"]]

## Setup covatiance matrices and results storage
# cov matrix theta (global)
nparam=length(theta) 
npc=c(1,1,1,1,1)
cov_matrix_thetaAll = diag(npc)
  colnames(cov_matrix_thetaAll)=names(theta)
  rownames(cov_matrix_thetaAll)=names(theta)

## thetaAlltab theta (local)
thetaAlltab=array(NA, dim=c(MCMC.runs+1,locnn,length(theta)),dimnames=list(NULL,NULL,names(theta)))
thetaAlltab[1,1,]=theta

# cov matrix theta_init
npc=c(0,0,1,0)
cov_matrix_thetainit = diag(npc)
  colnames(cov_matrix_thetainit)=names(theta_init)
  rownames(cov_matrix_thetainit)=names(theta_init)

# theta_init_Alltab
thetainitAlltab=array(NA, dim=c(MCMC.runs+1,locnn,length(theta_init)),dimnames=list(NULL,NULL,names(theta_init)))
thetainitAlltab[1,1,]=theta_init

# Arrays for storing results from MCMC loop
accepttab=rep(NA,(MCMC.runs))
sim_liktab=rep(-Inf,(MCMC.runs+1))
max.length = length(time.vals)
c_trace_tab=array(NA, dim=c(MCMC.runs+1,locnn,max.length)) 

## setup prior distributions
var.prior <- 0.2
priorGamma<-function(x){dgamma(x,shape=5/(var.prior), scale=var.prior)} # Prior - infectious period
priorAlpha<-function(x){dgamma(x,shape=5.9/(var.prior), scale=var.prior)} # Prior - extrinsic incubation period

for (m in 1:MCMC.runs){
    # Scale COV matrices for resampling using error term epsilon
    if(m==1){
      epsilon0 = 0.001
      cov_matrix_thetaA=epsilon0*cov_matrix_thetaAll
    }else{
      epsilon0 = max(min(0.1,exp(log(epsilon0)+(accept_rate-0.234)*0.999^m)),1e-6) # Stop epsilon getting too big or small
      cov_matrix_thetaA=epsilon0*cov_matrix_thetaAll
    }
    
    prior.star=1
    sim_marg_lik_star=0
    thetaAllstar=0*thetaAlltab[m,,]
    thetaInitstar=0*thetainitAlltab[m,,]
    cTraceStar=0*c_trace_tab[m,,]
    #kk=1
    for(kk in itertab){ 
      iiH=kk
      data <- load.data.multistart(add.nulls = 0, startdate=as.Date("2013-10-27"), virusTab[iiH], dataTab[iiH], serology.excel, init.conditions.excel)
        list2env(data,globalenv())
      
      if(m==1){ # Don't sample theta on 1st step 
        theta_cand = theta
        theta_init_cand = theta_init
      }else{
        theta_cand = as.numeric(exp(rmvnorm(1, log(thetaAlltab[m,iiH,]), cov_matrix_thetaA)))
        names(theta_cand) = names(theta)
        if(sum(names(theta_cand)=="rep")>0){ # check theta contains this vector
          theta_cand[["rep"]]=min(theta_cand[["rep"]],2-theta_cand[["rep"]]) # Ensure reporting between zero and 1
        }
        theta_init_cand = as.numeric(exp(rmvnorm(1, log(thetainitAlltab[m,iiH,]), cov_matrix_thetainit)))
        names(theta_init_cand) = names(theta_init)
      } 
      thetaA_star=theta_cand
      thetaI_star=theta_init_cand
      
      # Run model simulation
      output1 <- Run_simulation(dt=7,thetaA_star,thetaI_star,time.vals)
        casecount <- output1$C1_trace
        likelihood <- sum(log(dnbinom(y.vals,mu=thetaA_star[["rep"]]*(casecount),size=1/thetaA_star[["repvol"]])))
        sim_marg_lik_star=likelihood
      
      #Store vales
      thetaAllstar[iiH,]=thetaA_star
      thetaInitstar[iiH,]=thetaI_star
      cTraceStar[iiH,]=casecount 
    } 
    # Compute prior for current and proposal theta
    prior.theta=priorAlpha(1/thetaAlltab[1,1,'alpha'])*priorAlpha(1/thetaAlltab[1,1,'gamma'])
    prior.star=priorAlpha(1/thetaAllstar[1,'alpha'])*priorAlpha(1/thetaAllstar[1,'gamma'])
    
    # Calculate probability function 
    q_theta_given_theta_star = 1
    q_theta_star_given_theta = 1
    
    val = exp((sim_marg_lik_star-sim_liktab[m]))*(prior.star/prior.theta)*(q_theta_given_theta_star/q_theta_star_given_theta) 
    
    if(is.na(val)){
      output_prob=0
    }else{
      output_prob = min(val, 1)
    }
    
    # Update parameter values
    if(runif(1,0,1) < output_prob){
      thetaAlltab[m+1,,] = thetaAllstar
      thetainitAlltab[m+1,,] = thetaInitstar
      c_trace_tab[m+1,,]=cTraceStar
      sim_liktab[m+1] = sim_marg_lik_star
      accepttab[m]=1
    }else{
      thetaAlltab[m+1,,] = thetaAlltab[m,,]
      thetainitAlltab[m+1,,] = thetainitAlltab[m,,]
      c_trace_tab[m+1,,]= c_trace_tab[m,,]
      sim_liktab[m+1] = sim_liktab[m]
      accepttab[m]=0
    }
    
    if(m<20){
      accept_rate=0.234
    }else{
      accept_rate=sum(accepttab[1:m])/m
    }
    
  } # end of MCMC loop
    
  return(list(sim_liktab=sim_liktab,
                accepttab=accepttab,
                c_trace_tab=c_trace_tab,
                thetaAlltab=thetaAlltab,
                thetainitAlltab=thetainitAlltab
              ))
} # end of denv_MCMC function
