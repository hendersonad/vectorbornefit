#' ODEs for SEIR-SEI 
#' 
#' This function solves an SEIR-SEI vector borne disease model.
#' @param theta Vector of parameters for model
#' @param init.state List of initial values
#' @param time.vals.sim List of time values to simulate over
#' @export

simulate_deterministic <- function(theta, init.state, time.vals.sim) {
  SIR_ode <- function(time, state, theta) {
    ## extract parameters from theta
    beta_h1 <-  theta[["beta_h"]] *  seasonal_f(time, date0=0,amp=theta[["beta_v_amp"]],mid=theta[["beta_v_mid"]]) * decline_f(time,date0=theta[["shift_date"]],mask=0,base=1,grad=0.5,mid=0) 
    beta_h3 <-  theta[["beta_h_3"]] * beta_h1
    beta_h2 <-  theta[["beta_h_2"]] * beta_h3 
    beta_v1 <-  theta[["beta_v"]] * beta_h1 
    beta_v2 <-  theta[["beta_v"]] * beta_h2 
    beta_v3 <-  theta[["beta_v"]] * beta_h3 
    Nsize <-   theta[["npop"]]
    NsizeC <-   theta[["npopC"]]
    NsizeA <-   theta[["npopA"]]
    delta_v  <- theta[["MuV"]]
    alpha_v <-  theta[["Vex"]]
    alpha_h <-  theta[["Exp"]]
    gamma <-    theta[["Inf"]]
    
    ## extract initial states from theta_init
    SC <- state[["s_initC"]]
    EC <- state[["e_initC"]]
    IC <- state[["i_initC"]]
    RC <- state[["r_initC"]]
    CC <- state[["c_initC"]] 
    SMC <- state[["sm_initC"]]
    EMC <- state[["em_initC"]]
    IMC <- state[["im_initC"]]
    
    SA <- state[["s_initA"]]
    EA <- state[["e_initA"]]
    IA <- state[["i_initA"]]
    RA <- state[["r_initA"]]
    CA <- state[["c_initA"]] 
    SMA <- state[["sm_initA"]]
    EMA <- state[["em_initA"]]
    IMA <- state[["im_initA"]]
    
    # Human child population
    dSC  =  - SC*(beta_h1*IMC) 
    dEC  =  SC*(beta_h1*IMC) - alpha_h*EC  
    dIC  = alpha_h*EC  - gamma*IC
    dRC  = gamma*IC
    dCC  = alpha_h*EC 
    
    # Human adult population
    dSA  = - SA*(beta_h3*IMC)      
    dEA  =  SA*(beta_h3*IMC) - alpha_h*EA
    dIA  = alpha_h*EA  - gamma*IA
    dRA  = gamma*IA
    dCA  = alpha_h*EA
    
    # Mosquito populations
    dSMC = delta_v - SMC*(beta_v1*(IC+IA)/Nsize ) - delta_v*SMC   
    dEMC = SMC*(beta_v1*(IC+IA)/Nsize) - (delta_v+alpha_v)*EMC  
    dIMC = alpha_v*EMC-delta_v*IMC
    
    dSMA = delta_v - SMA*(beta_v2*IC/NsizeC+beta_v3*IA/NsizeA) - delta_v*SMA  
    dEMA = SMA*(beta_v2*IC/NsizeC+beta_v3*IA/NsizeA) - (delta_v+alpha_v)*EMA
    dIMA = alpha_v*EMA-delta_v*IMA
    
    
    return(list(c(dSC,dEC,dIC,dRC,dCC,dSMC,dEMC,dIMC,
                  dSA,dEA,dIA,dRA,dCA,dSMA,dEMA,dIMA)))
  }
  #Solve ODEs
  traj <- as.data.frame(ode(init.state, time.vals.sim, SIR_ode, theta, method = "ode45"))
  return(traj)
}

