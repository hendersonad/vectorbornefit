#' Sample the vector of parameters from a multivariate normal distribution
#' 
#' Sample a new vector of parameters (either global or local) and initial conditions using a multivariate normal dist
#' @param theta_in Vector of previous values of parameters
#' @param theta_init_in Vector of previous values of initial conditions
#' @param covartheta Covariance matrix for parameters. Must be same width as length of theta_in
#' @param covartheta_init Covariance matrix for initial conditions. Must be same width as length of theta_init_in
#' @param agestructure Binary indicator variable if model is age structured (1) or not (0) between children and adults. Defaults to NULL in which case no sampling happens
#' @param global Binary indicator variable if sampling global (1) or local (0) parameters. Defaults to NULL in which case no sampling happens
#' @export
#' @examples
#' SampleTheta(thetaAlltab[m,iiH,],theta_initAlltab[m,iiH,],cov_matrix_thetaA,cov_matrix_theta_init,agestructure,global=0)
#' SampleTheta(thetatab[m,],theta_initAlltab[1,1,],cov_matrix_theta,cov_matrix_theta_init,agestructure,global=1)

#theta_in=thetaAlltab[m,iiH,]; theta_init_in=theta_initAlltab[m,iiH,]; covartheta=cov_matrix_thetaA; covartheta_init=cov_matrix_theta_init

SampleThetaSIR<-function(theta_in, theta_init_in, covartheta, covartheta_init, agestructure=NULL, global=NULL){
  ## Parameters
    # sample new parameters from nearby: 
      mean_vector_theta = theta_in
      mean_vector_theta0=mean_vector_theta
      theta_star = as.numeric(exp(rmvnorm(1,log(mean_vector_theta0), covartheta)))
      names(theta_star)=names(theta_in)
      
      if(sum(names(theta_star)=="rep")>0){ # check theta contains this vector
        theta_star[["rep"]]=min(theta_star[["rep"]],2-theta_star[["rep"]]) # Ensure reporting between zero and 1
      }
      
      if(sum(names(theta_star)=="beta_v_amp")>0){
        theta_star[["beta_v_amp"]]=min(theta_star[["beta_v_amp"]],2-theta_star[["beta_v_amp"]]) # Ensure amplitude between zero and 1
      }
      
      if(sum(names(theta_star)=="beta_h_2")>0){
        theta_star[["beta_h_2"]]=min(theta_star[["beta_h_2"]],2-theta_star[["beta_h_2"]]) # Ensure beta is between zero and 1
        theta_star[["beta_h_3"]]=min(theta_star[["beta_h_3"]],2-theta_star[["beta_h_3"]]) # Ensure beta is between zero and 1
      }
  
    ## Initial conditions
    mean_vector_theta_init = theta_init_in
    
    ## Sample initial conditions from mulitvariate normal
    theta_init_star = as.numeric(exp(rmvnorm(1,log(mean_vector_theta_init), covartheta_init)))
    names(theta_init_star)=names(theta_init_in)
    
    if(global==0){
        # Human initial conditions: Sample I, set E equal to it and S is residual (Pop-E-I-R)
        theta_init_star[["r_init"]]= min(theta_init_star[["r_init"]], 2*thetatab[m,"npop"] - theta_init_star[["r_init"]]) # Bound at population size
        theta_init_star[["s_init"]] = (thetatab[m,"npop"]-theta_init_star[["i1_init"]]-theta_init_star[["r_init"]])
        theta_init_star=theta_init_star
    }else if(global==1){
      theta_init_star=theta_init_in
    }
  return(list(thetaS=theta_star,theta_initS=theta_init_star))
}
