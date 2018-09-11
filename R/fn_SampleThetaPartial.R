#' Sample the vector of parameters from a multivariate normal distribution
#' 
#' Sample a new vector of parameters (either global or local) and initial conditions using a multivariate normal dist
#' @param theta_in Vector of previous values of parameters
#' @param theta_init_in Vector of previous values of initial conditions
#' @param covartheta Covariance matrix for parameters. Must be same width as length of theta_in
#' @param covartheta_init Covariance matrix for initial conditions. Must be same width as length of theta_init_in
#' @param global Binary indicator variable if sampling global (1) or local (0) parameters. Defaults to NULL in which case no sampling happens
#' @export

SampleThetaPartial<-function(theta_in, theta_init_in, covartheta, covartheta_init, global=NULL){
  ## Parameters
    # sample new parameters from nearby using multivariate normal distribution: 
      mean_vector_theta = theta_in
      mean_vector_theta0 = mean_vector_theta
      theta_star = as.numeric(exp(rmvnorm(1,log(mean_vector_theta0), covartheta)))
      names(theta_star)=names(theta_in)
      
      if(sum(names(theta_star)=="rep")>0){ # check theta contains this vector
        theta_star[["rep"]]=min(theta_star[["rep"]],2-theta_star[["rep"]]) # Ensure reporting between zero and 1
      }
      if(sum(names(theta_star)=="iota")>0){ # check theta contains this vector
        theta_star[["iota"]]=min(theta_star[["iota"]],2-theta_star[["iota"]]) # Ensure reporting between zero and 1
      }
      #
      #if(sum(names(theta_star)=="beta_v_amp")>0){
      #  theta_star[["beta_v_amp"]]=min(theta_star[["beta_v_amp"]],2-theta_star[["beta_v_amp"]]) # Ensure amplitude between zero and 1
      #}
      #
      #if(sum(names(theta_star)=="beta_h_2")>0){
      #  theta_star[["beta_h_2"]]=min(theta_star[["beta_h_2"]],2-theta_star[["beta_h_2"]]) # Ensure beta is between zero and 1
      #  theta_star[["beta_h_3"]]=min(theta_star[["beta_h_3"]],2-theta_star[["beta_h_3"]]) # Ensure beta is between zero and 1
      #}
      #
      #if(sum(names(theta_star)=="t0")>0){ # check theta contains this vector
      #  theta_star[["t0"]]=max(0,theta_star[["t0"]]) 
      #}
      #
      #if(sum(names(theta_star)=="chi")>0){ # check theta contains this vector
      #  theta_star[["chi"]]=min(theta_star[["chi"]],2-theta_star[["chi"]]) # Ensure reporting between zero and 1
      #}
      
  ## Initial conditions
  theta_init_star = theta_init_in
    
  initial_inf=as.numeric(theta_star['inf0'])
  init_vec=as.numeric(theta_star['vec0']/2)
  
  popsizeTot=theta_init_star["s_init"]+theta_init_star["e_init"]+theta_init_star["i1_init"]+theta_init_star["r_init"]

  theta_init_star["r_init"]=0
  theta_init_star["e_init"]=initial_inf; theta_init_star["i1_init"]=initial_inf
  theta_init_star["em_init"]=init_vec; theta_init_star["im_init"]=init_vec
  
  theta_init_star["s_init"]=popsizeTot-theta_init_star["i1_init"]-theta_init_star["e_init"]-theta_init_star["r_init"]
  theta_init_star["sm_init"]=1-theta_init_star["em_init"]-theta_init_star["im_init"]
  
  theta_init_star["ed_init"]=0; theta_init_star["id_init"]=thetainit_denv[["i1_init"]]; theta_init_star["t1d_init"]=0; theta_init_star["t1d_init"]=0
  theta_init_star["sd_init"]=popsizeTot-theta_init_star["id_init"]-theta_init_star["ed_init"]-theta_init_star["t1d_init"]-theta_init_star["t2d_init"]
  
  theta_init_star["s2_init"]=0
  theta_init_star["e2_init"]=0
  theta_init_star["i1_2_init"]=0
  theta_init_star["r2_init"]=0
  
  return(list(thetaS=theta_star,theta_initS=theta_init_star))
}
