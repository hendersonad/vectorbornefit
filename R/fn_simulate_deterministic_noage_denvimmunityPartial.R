#' ODEs for SEIR-SEI with cross immunity from denv infection
#' 
#' This function solves an SEIR-SEI vector borne disease model.
#' @param theta Vector of parameters for model
#' @param init.state List of initial values
#' @param time.vals.sim List of time values to simulate over
#' @export
#theta=theta; init.state=init1; time.vals.sim=time.vals

simulate_deterministic_noage_DENVimm_partial <- function(theta, init.state, time.vals.sim) {
  #time=1;state=init1;
  SIR_ode <- function(time, state, theta) {
    ## extract parameters from theta
    Nsize <-   theta[["npop"]]
    chi <-    theta[["chi"]]
    psi <- theta[["psi"]]
    
    # No DENV outbreak until denv_start parameter reached in time.vals
    if(time<theta[["denv_end"]]){ 
      chi <- theta[["chi"]]
    }else{
      chi <- 0
    }
    
    # OR no ZIKV outbreak until zikv_start reached in time.vals
    if(time<theta[["zika_start"]]){ 
      beta_h1 <- 0
      beta_v1 <- 0
      delta_v <- 0
      alpha_v <- 0
      alpha_h <- 0
      gamma <- 0
      rho <- 0
    }else{
      beta_h1 <- theta[['beta_h']] * seasonal_f(time, date0=theta[["shift_date"]],amp=theta[["beta_v_amp"]],mid=theta[["beta_v_mid"]]) * decline_f(time,date0=theta[["shift_date"]],mask=theta[['beta_mask']],base=theta[['beta_base']],grad=theta[['beta_grad']],mid=theta[['beta_mid']]) 
      beta_v1 <- theta[['beta_v']] * beta_h1 
      delta_v  <- theta[["MuV"]] 
      alpha_v <-  theta[["Vex"]]
      alpha_h <-  theta[["Exp"]]
      gamma <-    theta[["Inf"]]
      rho <-      1/theta[["rho"]]
    }

    ## extract initial states from theta_init
    S <- state[["s_init"]]
    E <- state[["e_init"]]
    I <- state[["i_init"]]
    R <- state[["r_init"]]
    S2 <- state[["s2_init"]]
    E2 <- state[["e2_init"]]
    I2 <- state[["i2_init"]]
    R2 <- state[["r2_init"]]
    C <- state[["c_init"]] 
    SM <- state[["sm_init"]]
    EM <- state[["em_init"]]
    IM <- state[["im_init"]]
    D3 <- state[["D3_init"]]
    FP <- state[["FP_init"]]
      
    ## extinction if not at least 1 infected
    Ipos = extinct(I,1) # Need at least one infective
    
    # Human population
    dS  = - S*(beta_h1*IM)*Ipos - (D3*(1/7))
    dE  =  S*(beta_h1*IM)*Ipos - alpha_h*E
    dI  = alpha_h*E  - gamma*I + (psi*FP)
    dR  = gamma*I - rho*R
    
    # Human population - that has had DENV
    dS2  = (D3*(1/7)) #- S2*((1-chi)*beta_h1*IM)*Ipos 
    dE2  = S2*((1-chi)*beta_h1*IM)*Ipos - alpha_h*E2
    dI2  = alpha_h*E2  - gamma*I2
    dR2  = gamma*I2 - rho*R2
    
    dC  = alpha_h*(E+E2) 
    
    # Mosquito population
    dSM = delta_v - SM*(beta_v1*(I+I2)/Nsize)*Ipos - delta_v*SM   
    dEM = SM*(beta_v1*(I+I2)/Nsize)*Ipos - (delta_v+alpha_v)*EM  
    dIM = alpha_v*EM-delta_v*IM
    
    # underlying infections
    dD3 = D3reg(time) - D3
    dFP = Ctreg(time) - FP
    
    return(list(c(dS,dE,dI,dR,dS2,dE2,dI2,dR2,dC,dSM,dEM,dIM,dD3,dFP)))
  }
  #init.state=c(init1,"D3_init"=D3reg(1),"FP_init"=Ctreg(1))
  #time.vals.sim=time.vals
  traj <- as.data.frame(ode(init.state, time.vals.sim, SIR_ode, theta, method = "ode45"))
  return(traj)
}

