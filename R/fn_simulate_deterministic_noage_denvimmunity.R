#' ODEs for SEIR-SEI with cross immunity from denv infection
#' 
#' This function solves an SEIR-SEI vector borne disease model.
#' @param theta Vector of parameters for model
#' @param init.state List of initial values
#' @param time.vals.sim List of time values to simulate over
#' @export
#theta=theta; init.state=init1; time.vals.sim=time.vals; state=init.state

simulate_deterministic_noage_DENVimm <- function(theta, init.state, time.vals.sim) {
  #time=1;state=init1;
  SIR_ode <- function(time, state, theta) {
    ## extract parameters from theta
    Nsize <-   theta[["npop"]]
    if(!is.na(theta[['omega_d']])){
      omega_d <- 1/theta[['omega_d']]}else{
      omega_d <- 0}
    if(!is.na(theta[['chi']])){
      chi <- theta[['chi']]}else{
      chi <- 0}
    if(!is.na(theta[['psi']])){
      psi <- theta[['psi']]}else{
      psi <- 0}
    
    # No of ZIKA FP introductions
    if(time<=(18*7)){
      Ct <- Ctreg(time)}else{ #no of zika infections
      Ct=0}
    #psi <- 2e-4
    #plot(psi*sapply(time.vals.sim, function(x){Ctreg(x)}), type='l')
    #sum(sapply(time.vals.sim, function(x){Ctreg(x)}))
    # No DENV outbreak until denv_start parameter reached in time.vals
    if(time<theta[["denv_start"]]){ 
      beta_d <- 0
      alpha_d <- 0 
      gamma_d <- 0
    }else{
      beta_d <- theta[['beta_d']]
      alpha_d <- theta[['alpha_d']]
      gamma_d <- theta[['gamma_d']]
    }
    # And no ZIKV outbreak until zikv_start reached in time.vals
    if(time<theta[["zika_start"]]){ 
      beta_h1 <- 0
      beta_v1 <- 0
      delta_v <- 0
      alpha_v <- 0
      alpha_h <- 0
      gamma <- 0
      rho <- 0
    }else{
      beta_h1 <- theta[['beta_h']] * 
                  seasonal_f(time, date0=theta[["shift_date"]],amp=theta[["beta_v_amp"]],mid=theta[["beta_v_mid"]]) * 
                  decline_f(time, mid=theta[["beta_mid"]], width=theta[["beta_grad"]], base=theta[["beta_base"]])
      beta_v1 <- theta[['beta_v']] * beta_h1 
      delta_v  <- theta[["MuV"]] 
      alpha_v <-  theta[["Vex"]]
      alpha_h <-  theta[["Exp"]]
      gamma <-    theta[["Inf"]]
      if(!is.na(theta[['rho']])){
        rho <- 1/theta[['rho']]}else{
        rho <- 0}    
    }

    ## extract initial states from theta_init
    S <- state[["s_init"]]
    E <- state[["e_init"]]
    I <- state[["i_init"]]
    R <- state[["r_init"]]
    C <- state[["c_init"]] 
    Sd <- state[["sd_init"]]
    Ed <- state[["ed_init"]]
    Id <- state[["id_init"]]
    T1d <- state[["t1d_init"]]
    T2d <- state[["t2d_init"]]
    Cd <- state[["cd_init"]]
    SM <- state[["sm_init"]]
    EM <- state[["em_init"]]
    IM <- state[["im_init"]]
    
    ## extinction if not at least 1 infected
    Ipos = extinct(I,1) # Need at least one infective
    
    # Human population
    dS  =  - S*(beta_h1*IM)*Ipos - chi*Sd*(beta_d*Id/Nsize) + chi*(2*omega_d*T2d) 
    dE  =  S*(beta_h1*IM)*Ipos - alpha_h*E  
    dI  = alpha_h*E  - gamma*I + (psi*Ct)
    dR  = gamma*I - rho*R
    dC  = alpha_h*E 
    
    # Denv infection and immunity
    dSd = -Sd*(beta_d*Id/Nsize)
    dEd = Sd*(beta_d*Id/Nsize) - alpha_d*Ed 
    dId = alpha_d*Ed - gamma_d*Id  
    dT1d = gamma_d*Id - 2*omega_d*T1d
    dT2d = 2*omega_d*T1d - 2*omega_d*T2d
    dCd = alpha_d*Ed
    
    # Mosquito population
    dSM = delta_v - SM*(beta_v1*I/Nsize)*Ipos - delta_v*SM   
    dEM = SM*(beta_v1*I/Nsize)*Ipos - (delta_v+alpha_v)*EM  
    dIM = alpha_v*EM-delta_v*IM
    
    return(list(c(dS,dE,dI,dR,dC,dSd,dEd,dId,dT1d,dT2d,dCd,dSM,dEM,dIM)))
  }
  traj <- as.data.frame(ode(init.state, time.vals.sim, SIR_ode, theta, method = "ode45"))
  return(traj)
}

