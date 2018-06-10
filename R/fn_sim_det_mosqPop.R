#' ODEs for SEIR-SEI with noage structure
#' 
#' This function solves an SEIR-SEI vector borne disease model.
#' @param theta Vector of parameters for model
#' @param init.state List of initial values
#' @param time.vals.sim List of time values to simulate over
#' @export

simulate_deterministic_mosqPop <- function(theta, init.state, time.vals.sim) {
  SIR_ode <- function(time, state, theta) {
    ## extract parameters from theta
    beta_h1 <-  theta[["beta_h"]] * seasonal_f(time, date0=theta[["shift_date"]],amp=theta[["beta_v_amp"]],mid=theta[["beta_v_mid"]]) 
    #beta_h1 <-  theta[["beta_h"]] 
    beta_v1 <-  theta[["beta_v"]] * beta_h1 
    Nsize   <-  theta[["npop"]]
    nu_v    <-  theta[["MuV"]] 
    delta_v <-  0 #theta[["MuV"]] # * death_f(time,base=theta[['beta_base']])
    iota_v  <-  theta[["MuV"]] * death_f(time,base=theta[['beta_base']]) # the iota (i) stands for intervention
    alpha_v <-  theta[["Vex"]] 
    alpha_h <-  theta[["Exp"]]
    gamma   <-  theta[["Inf"]]
    
    ## extract initial states from theta_init
    S <- state[["s_init"]]
    E <- state[["e_init"]]
    I <- state[["i_init"]]
    R <- state[["r_init"]]
    C <- state[["c_init"]] 
    SM <- state[["sm_init"]]
    EM <- state[["em_init"]]
    IM <- state[["im_init"]]
    
    ## mosquito pop size
    #MNsize <-   SM+EM+IM
    MNsize <-   Nsize*theta[["mosq_ratio"]]

    # Human population
    dS  =  - S*(beta_h1*IM/MNsize) 
    dE  =  S*(beta_h1*IM/MNsize) - alpha_h*E  
    dI  = alpha_h*E  - gamma*I
    dR  = gamma*I
    dC  = alpha_h*E 
    
    # Mosquito population
    dSM = (nu_v*(MNsize)) - SM*(beta_v1*I/Nsize) - delta_v*SM - iota_v*SM
    dEM = SM*(beta_v1*I/Nsize) - (delta_v+alpha_v+iota_v)*EM
    dIM = alpha_v*EM - delta_v*IM - iota_v*IM
    
    return(list(c(dS,dE,dI,dR,dC,dSM,dEM,dIM)))
  }
  #Solve ODEs
  traj <- as.data.frame(ode(init.state, time.vals.sim, SIR_ode, theta, method = "ode45"))

  return(traj)
}
