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

#theta_in=thetaAlltab[m,iiH,]; theta_init_in=theta_initAlltab[m,iiH,]; covartheta=cov_matrix_thetaA; covartheta_init=cov_matrix_theta_init; global=1
#theta_in=thetaAlltab[1,1,]; theta_init_in=theta_initAlltab[1,1,]; covartheta=cov_matrix_thetaAll; covartheta_init=cov_matrix_theta_initAll; global=0
 
SampleTheta<-function(theta_in, theta_init_in, covartheta, covartheta_init, agestructure=NULL, global=NULL){
  ## Parameters
    # sample new parameters from nearby: 
      mean_vector_theta = theta_in
      mean_vector_theta0 = mean_vector_theta
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
      
      if(sum(names(theta_star)=="t0")>0){ # check theta contains this vector
        theta_star[["t0"]]=min(theta_star[["t0"]],2-theta_star[["t0"]]) # Ensure start time between zero and 1 (therefore log is negative)
      }
      
    ## Initial conditions
    theta_init_star = theta_init_in
    
  initial_inf=as.numeric(theta_star['inf0']) #*(popsizeTot/2))
  init_vec=as.numeric(theta_star['vec0']/2)
  
  popsizeTot=theta_init_star["s_init"]+theta_init_star["e_init"]+theta_init_star["i1_init"]+theta_init_star["r_init"]

  # initial recovered 
  if(agestructure==1){
    theta_init_star["r_initC"]=0
    theta_init_star["r_initA"]=0
      theta_init_star["e_initC"]=initial_inf; theta_init_star["i1_initC"]=initial_inf
      theta_init_star["em_initC"]=init_vec; theta_init_star["im_initC"]=init_vec
      theta_init_star["e_initA"]=initial_inf; theta_init_star["i1_initA"]=initial_inf
      theta_init_star["em_initA"]=init_vec; theta_init_star["im_initA"]=init_vec
      
      theta_init_star["s_initC"]=popsizeC-theta_init_star["i1_initC"]-theta_init_star["e_initC"]-theta_init_star["r_initC"]
      theta_init_star["sm_initC"]=1-theta_init_star["em_initC"]-theta_init_star["im_initC"]
      
      theta_init_star["s_initA"]=popsizeA-theta_init_star["i1_initA"]-theta_init_star["e_initA"]-theta_init_star["r_initA"]
      theta_init_star["sm_initA"]=1-theta_init_star["em_initA"]-theta_init_star["im_initA"]    
  }else{
    if(baselineSero==T){
      theta_init_star["r_init"]=(nLUM[1]/nPOP[1])*popsizeTot
    }else{theta_init_star["r_init"]=0}
      
      theta_init_star["e_init"]=initial_inf; theta_init_star["i1_init"]=initial_inf
      theta_init_star["em_init"]=init_vec; theta_init_star["im_init"]=init_vec
      
      theta_init_star["s_init"]=popsizeTot-theta_init_star["i1_init"]-theta_init_star["e_init"]-theta_init_star["r_init"]
      theta_init_star["sm_init"]=1-theta_init_star["em_init"]-theta_init_star["im_init"]
      
      theta_init_star["ed_init"]=0; theta_init_star["id_init"]=theta_init_star[["i1_init"]]; theta_init_star["rd_init"]=0
      theta_init_star["sd_init"]=popsizeTot-theta_init_star["id_init"]-theta_init_star["ed_init"]-theta_init_star["rd_init"]
  }


#    ## Sample initial conditions from mulitvariate normal
#    theta_init_star = as.numeric(exp(rmvnorm(1,log(mean_vector_theta_init), covartheta_init)))
#    names(theta_init_star)=names(theta_init_in)
#    
#    if(global==0){
#      if(agestructure==1){
#      # Mosquito init conditions: Sample I, fix between 0 and 1, set E=I, S is residual (1-E-I)
#      theta_init_star[["im_initC"]] = min(theta_init_star[["im_initC"]],1-theta_init_star[["im_initC"]])
#        theta_init_star[["em_initC"]] = theta_init_star[["im_initC"]] 
#        theta_init_star[["sm_initC"]] = 1-theta_init_star[["im_initC"]]-theta_init_star[["em_initC"]]
#      theta_init_star[["im_initA"]] = min(theta_init_star[["im_initA"]],1-theta_init_star[["im_initA"]])
#        theta_init_star[["em_initA"]] = theta_init_star[["im_initA"]]
#        theta_init_star[["sm_initA"]] = 1-theta_init_star[["im_initA"]]-theta_init_star[["em_initA"]]
#      
#      # Human initial conditions: Sample I, set E equal to it and S is residual (Pop-E-I-R)
#      theta_init_star[["r_initC"]]= min(theta_init_star[["r_initC"]], 2*thetatab[m,"npopC"] - theta_init_star[["r_initC"]]) # Bound at population size
#      theta_init_star[["r_initA"]]= min(theta_init_star[["r_initA"]], 2*thetatab[m,"npopA"] - theta_init_star[["r_initA"]]) # Bound at population size
#      
#      theta_init_star[["e_initC"]] = theta_init_star[["i1_initC"]]
#      theta_init_star[["e_initA"]] = theta_init_star[["i1_initA"]]
#      
#      theta_init_star[["s_initC"]] = (thetatab[m,"npopC"]-theta_init_star[["i1_initC"]]-theta_init_star[["e_initC"]]-theta_init_star[["r_initC"]])
#      theta_init_star[["s_initA"]] = (thetatab[m,"npopA"]-theta_init_star[["i1_initA"]]-theta_init_star[["e_initA"]]-theta_init_star[["r_initA"]])
#  
#      theta_init_star=theta_init_star
#      }else if(agestructure==0){
#        # Mosquito init conditions: Sample I, fix between 0 and 1, set E=I, S is residual (1-E-I)
#        if(mosquitoes=="density"){
#          mosqpop=1
#        }else{
#          mosqpop=thetatab[m,"npop"]*thetatab[m,"mosq_ratio"]
#        }
#        theta_init_star[["im_init"]] = min(theta_init_star[["im_init"]],mosqpop-theta_init_star[["im_init"]])
#        theta_init_star[["em_init"]] = theta_init_star[["im_init"]] 
#        theta_init_star[["sm_init"]] = mosqpop-theta_init_star[["im_init"]]-theta_init_star[["em_init"]]
#        # Human initial conditions: Sample I, set E equal to it and S is residual (Pop-E-I-R)
#        theta_init_star[["r_init"]]= min(theta_init_star[["r_init"]], 2*thetatab[m,"npop"] - theta_init_star[["r_init"]]) # Bound at population size
#        theta_init_star[["e_init"]] = theta_init_star[["i1_init"]]
#        theta_init_star[["s_init"]] = (thetatab[m,"npop"]-theta_init_star[["i1_init"]]-theta_init_star[["e_init"]]-theta_init_star[["r_init"]])
#        theta_init_star=theta_init_star
#      }
#    }else if(global==1){
#      theta_init_star=theta_init_in
#    }
  return(list(thetaS=theta_star,theta_initS=theta_init_star))
}
