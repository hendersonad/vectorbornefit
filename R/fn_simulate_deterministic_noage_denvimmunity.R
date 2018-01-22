#' ODEs for SEIR-SEI with noage structure
#' 
#' This function solves an SEIR-SEI vector borne disease model.
#' @param theta Vector of parameters for model
#' @param init.state List of initial values
#' @param time.vals.sim List of time values to simulate over
#' @export
#theta=final_thetaAll; init.state=theta_init_star; time.vals.sim=time.vals

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
    chi <-    theta[["chi"]]
    
    beta_d <- theta[['beta_d']]
    alpha_d <- theta[['alpha_d']]
    gamma_d <- theta[['gamma_d']]
    omega_d <- theta[['omega_d']]
    
    ## extract initial states from theta_init
    S <- state[["s_init"]]
    E <- state[["e_init"]]
    I <- state[["i_init"]]
    R <- state[["r_init"]]
    C <- state[["c_init"]] 
    Sd <- state[["sd_init"]]
    Ed <- state[["ed_init"]]
    Id <- state[["id_init"]]
    Rd <- state[["rd_init"]]
    Cd <- state[["cd_init"]]
    SM <- state[["sm_init"]]
    EM <- state[["em_init"]]
    IM <- state[["im_init"]]
    
    # Human population
    dS  =  - S*(beta_h1*IM) - chi*Sd*(beta_d*Id/Nsize) + chi*omega_d*Rd
    dE  =  S*(beta_h1*IM) - alpha_h*E  
    dI  = alpha_h*E  - gamma*I
    dR  = gamma*I
    dC  = alpha_h*E 
    
    # Denv infection and immunity
    dSd = -Sd*(beta_d*Id/Nsize)
    dEd = Sd*(beta_d*Id/Nsize) - alpha_d*Ed 
    dId = alpha_d*Ed - gamma_d*Id  
    dRd = gamma_d*Id - omega_d*Rd
    dCd = alpha_d*Ed
    
    # Mosquito population
    dSM = delta_v - SM*(beta_v1*I/Nsize) - delta_v*SM   
    dEM = SM*(beta_v1*I/Nsize) - (delta_v+alpha_v)*EM  
    dIM = alpha_v*EM-delta_v*IM
    
    return(list(c(dS,dE,dI,dR,dC,dSd,dEd,dId,dRd,dCd,dSM,dEM,dIM)))
  }
  #Solve ODEs
  traj <- as.data.frame(ode(init.state, time.vals.sim, SIR_ode, theta, method = "ode45"))
  return(traj)
}
