##' Compute prior
##' 
##' @param iiH Location 
##' @param thetaAlltab Current values of theta
##' @param thetaAllstar Proposed values of theta
##' @export

#thetaAlltab=thetaAlltab[m,,]
#thetaAllstar=thetaAllstar
#thetaAlltab=c(thetatab[m,],thetaAlltab[m,iiH,])
#thetaAllstar=c(theta_star,thetaA_star)
ComputePrior <- function(iiH, thetaAlltab, thetaAllstar){

p_theta_star = priorInf(1/thetaAllstar["Inf"])*priorExp(1/thetaAllstar["Exp"])*
                priorVEx(1/thetaAllstar["Exp"])*priorMuV(1/thetaAllstar["MuV"])*
                priorBeta_amp(thetaAllstar["beta_v_amp"])*priorBeta_mid(thetaAllstar["beta_v_mid"])*
                priorchi(thetaAllstar["chi"])*
                priort0(thetaAllstar["t0"])
if(vector.control==T){
  p_theta_star <- p_theta_star*priorBeta_base(thetaAllstar["beta_base"])
}
p_theta = priorInf(1/thetaAlltab["Inf"])*priorExp(1/thetaAlltab["Exp"])*
            priorVEx(1/thetaAlltab["Exp"])*priorMuV(1/thetaAlltab["MuV"])*
            priorBeta_amp(thetaAlltab["beta_v_amp"])*priorBeta_mid(thetaAlltab["beta_v_mid"])*
            priorchi(thetaAlltab["chi"])*
            priort0(thetaAlltab["t0"])
if(vector.control==T){
  p_theta <- p_theta*priorBeta_base(thetaAlltab["beta_base"])
}

names(p_theta_star)=NULL;names(p_theta)=NULL

if(is.na(p_theta_star)==T){p_theta_star=0}
if(is.na(p_theta)==T){p_theta=0}

return(list(prior.star=p_theta_star,prior=p_theta))
}  
