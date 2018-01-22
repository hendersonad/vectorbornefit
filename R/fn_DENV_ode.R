#' ODEs for SEIR-SEI with noage structure
#' 
#' This function solves an SEIR-SEI vector borne disease model.
#' @param theta Vector of parameters for model
#' @param init.state List of initial values
#' @param time.vals.sim List of time values to simulate over
#' @export

DENV_ode <- function(theta, init.state, time.vals.sim) {
  SIR_ode <- function(time, state, theta) {
    ## extract parameters from theta
    beta_h1 <-  theta[["beta_h"]]
    Nsize <-   theta[["npop"]]
    gamma <-    theta[["Inf"]]
    sigma <-    theta[["sigma"]]
    
    ## extract initial states from theta_init
    S <- state[["s_init"]]
    I <- state[["i_init"]]
    R <- state[["r_init"]]
    C <- state[["c_init"]] 
    
    # Human population
    dS  = -S*(beta_h1*I/Nsize) + sigma*R
    dI  = S*(beta_h1*I/Nsize) - gamma*I
    dR  = gamma*I - sigma*R
    dC  = S*(beta_h1*I/Nsize)
    
    return(list(c(dS,dI,dR,dC)))
  }
  #Solve ODEs
  traj <- as.data.frame(ode(init.state, time.vals.sim, SIR_ode, theta, method = "ode45"))
  return(traj)
}

