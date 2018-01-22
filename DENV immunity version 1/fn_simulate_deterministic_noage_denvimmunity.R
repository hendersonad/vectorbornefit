#' ODEs for SEIR-SEI with noage structure
#' 
#' This function solves an SEIR-SEI vector borne disease model.
#' @param theta Vector of parameters for model
#' @param init.state List of initial values
#' @param time.vals.sim List of time values to simulate over
#' @export

simulate_deterministic_noage_DENVimm <- function(theta, init.state, time.vals.sim) {
  SIR_ode <- function(time, state, theta) {
    ## extract parameters from theta
    beta_h1 <-  theta[["beta_h"]] * seasonal_f(time, date0=theta[["shift_date"]],amp=theta[["beta_v_amp"]],mid=theta[["beta_v_mid"]]) * decline_f(time,date0=theta[["shift_date"]],mask=theta[['beta_mask']],base=theta[['beta_base']],grad=theta[['beta_grad']],mid=theta[['beta_mid']]) 
    beta_v1 <-  theta[["beta_v"]] * beta_h1 
    Nsize <-   theta[["npop"]]
    delta_v  <- theta[["MuV"]] 
    alpha_v <-  theta[["Vex"]]
    alpha_h <-  theta[["Exp"]]
    gamma <-    theta[["Inf"]]
    CrossImmunity <- CrossImmune(time,date0=0,mask=theta[["wane_mask"]],base=0.2,grad=2,mid=0.3)
    immunity <- CrossImmunity$imm
    waning <- CrossImmunity$wane
    #print(paste0(immunity, "  ", waning))
    ## extract initial states from theta_init
    S <- state[["s_init"]]
    E <- state[["e_init"]]
    I <- state[["i_init"]]
    R <- state[["r_init"]]
    C <- state[["c_init"]] 
    SM <- state[["sm_init"]]
    EM <- state[["em_init"]]
    IM <- state[["im_init"]]
    immune <- state[["immune"]]
    
    # Human population
    dS  =  - S*(beta_h1*IM) - immunity*S + waning*immune
    dimmune = immunity*S - waning*immune
    dE  =  S*(beta_h1*IM) - alpha_h*E  
    dI  = alpha_h*E  - gamma*I
    dR  = gamma*I
    dC  = alpha_h*E 
    
    # Mosquito population
    dSM = delta_v - SM*(beta_v1*I/Nsize) - delta_v*SM   
    dEM = SM*(beta_v1*I/Nsize) - (delta_v+alpha_v)*EM  
    dIM = alpha_v*EM-delta_v*IM
    return(list(c(dS,dimmune,dE,dI,dR,dC,dSM,dEM,dIM)))
  }
  #Solve ODEs
  traj <- as.data.frame(ode(init.state, time.vals.sim, SIR_ode, theta, method = "ode45"))
  return(traj)
}

 