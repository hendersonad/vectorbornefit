##' Compute prior
##' 
##' @param iiH Location 
##' @param thetaAlltab Current values of theta
##' @param thetaAllstar Proposed values of theta
##' @param parameters_est_file Name of file with list of parameters to estimate
##' @export

ComputePrior_special <- function(iiH, thetaAlltab, thetaAllstar, parameter_est_file){
  parameters_est <- read.csv(paste0("data_sets/",parameter_est_file,".csv"), stringsAsFactors = F) 
  parms_to_est <- parameters_est$parameters_est
  parms_to_est <- c(parms_to_est,"Inf.")
  
  p_theta_star=0; p_theta=0
  for(k in parms_to_est){
    if(sum(names(thetaAllstar[iiH,])==k)>0){
      proposal <- thetaAllstar[iiH,k]
      original <- thetaAlltab[iiH,k]
      if(k=='Vex'| k=='Exp' | k=='MuV' | k=='Inf'){
        proposal = 1/proposal
        original = 1/original
      }
      p_theta_star = p_theta_star + dnorm(proposal, mean=param1[k,1], sd=param1[k,1])
      p_theta = p_theta + dnorm(original, mean=param1[k,1], sd=param1[k,1])
    }else{
      p_theta_star = p_theta_star
      p_theta = p_theta
    }
  }
  
  names(p_theta_star)=NULL;names(p_theta)=NULL
  
  if(is.na(p_theta_star)==T){p_theta_star=0}
  if(is.na(p_theta)==T){p_theta=0}
  
  if(is.nan(p_theta_star)==T){p_theta_star=0}
  if(is.nan(p_theta)==T){p_theta=0}
  
  return(list(prior.star=p_theta_star,prior=p_theta))
}  
