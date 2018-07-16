##' Compute prior
##' 
##' @param iiH Numeric indicator of location given by location vector "locationtab[]"
##' @param thetaAlltab Current values of parameter set theta
##' @param thetaAllstar Proposed values of parameter set theta
##' @export
# thetaAlltab=thetaAlltab_current[iiH,];  thetaAllstar=thetaA_star; covartheta = cov_matrix_thetaA
ComputePrior <- function(iiH, thetaAlltab, thetaAllstar, covartheta){
  # Compute prior of proposed parameter set theta
  p_theta_star =  priorInf(1/thetaAllstar["Inf"])*
                  priorExp(1/thetaAllstar["Exp"])*
                  priorVEx(1/thetaAllstar["Vex"])*
                  priorMuV(1/thetaAllstar["MuV"])*
                  priorOmega(thetaAllstar["omega_d"])
  

  # If vector control is included in the model - include prior density of proposed 
  #   strength of control measure
  if(vector.control==T){
    p_theta_star <- p_theta_star*priorBeta_base(thetaAllstar["beta_base"])
  }
  # If seasonal transmission is estimated, include in prior
  if(seasonal.transmission==T){
    p_theta_star <- p_theta_star*                  
      priorBeta_amp(thetaAllstar["beta_v_amp"])*
      priorBeta_mid(thetaAllstar["beta_v_mid"])
  }
  
  # Compute prior of current parameter set theta
  p_theta = priorInf(1/thetaAlltab["Inf"])*
            priorExp(1/thetaAlltab["Exp"])*
            priorVEx(1/thetaAlltab["Vex"])*
            priorMuV(1/thetaAlltab["MuV"])*
            priorOmega(thetaAlltab["omega_d"])
  
  # If vector control is included in the model - include prior density of proposed 
  #   strength of control measure
  if(vector.control==T){
    p_theta <- p_theta*priorBeta_base(thetaAlltab["beta_base"])
  }
  # If seasonal transmission is estimated, include in prior
  if(seasonal.transmission==T){
    p_theta <- p_theta*                  
      priorBeta_amp(thetaAlltab["beta_v_amp"])*
      priorBeta_mid(thetaAlltab["beta_v_mid"])
  }
  # formatting and replace NAs with 0 
  names(p_theta_star)=NULL;names(p_theta)=NULL
  if(is.na(p_theta_star)==T){p_theta_star=0}
  if(is.na(p_theta)==T){p_theta=0}
  return(list(prior.star=p_theta_star,prior=p_theta))
}  
