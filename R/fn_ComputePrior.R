##' Compute prior
##' 
##' @param iiH Location 
##' @param thetaAlltab Current values of theta
##' @param thetaAllstar Proposed values of theta
##' @export

#thetaAlltab=thetaAlltab[m,,]
#thetaAllstar=thetaAllstar

ComputePrior <- function(iiH, thetaAlltab, thetaAllstar){

p_theta_star = priorInf(1/thetaAllstar[iiH,"Inf"])*priorExp(1/thetaAllstar[iiH,"Exp"])*
                priorVEx(1/thetaAllstar[iiH,"Exp"])*priorMuV(1/thetaAllstar[iiH,"MuV"])*
                priorBeta_amp(thetaAllstar[iiH,"beta_v_amp"])*priorBeta_mid(thetaAllstar[iiH,"beta_v_mid"])*
                priorchi(thetaAllstar[iiH,"chi"])*
                priort0(thetaAllstar[iiH,"t0"])
p_theta = priorInf(1/thetaAlltab[iiH,"Inf"])*priorExp(1/thetaAlltab[iiH,"Exp"])*
            priorVEx(1/thetaAlltab[iiH,"Exp"])*priorMuV(1/thetaAlltab[iiH,"MuV"])*
            priorBeta_amp(thetaAlltab[iiH,"beta_v_amp"])*priorBeta_mid(thetaAlltab[iiH,"beta_v_mid"])*
            priorchi(thetaAlltab[iiH,"chi"])*
            priort0(thetaAlltab[iiH,"t0"])

names(p_theta_star)=NULL;names(p_theta)=NULL

if(is.na(p_theta_star)==T){p_theta_star=0}
if(is.na(p_theta)==T){p_theta=0}

return(list(prior.star=p_theta_star,prior=p_theta))
}  
