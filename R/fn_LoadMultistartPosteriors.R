##' load.posteriors
##' 
##' Load and combine posterior chains following MCMC run
#' @param agestructure Binary indicator variable if model is age structured (1) or not (0) between children and adults. Defaults to NULL in which case no sampling happens 
#' @param iiH Locationtab index
#' @param mcmc.burn Between 0 and 1. Proportion of iterations to discard.
#' @param month.start Month the series started
##' @export 
##' 

load.multistart.posteriors <- function(agestructure, iiH, mcmc.burn, month.start){
  #m.tot=length(list.files(path = paste("post_output_multistart/",month.start,"/",sep="")))
  m.tot=2
  thetatabA=NULL
  theta_inittabA=NULL
  c_trace_tab0=NULL
  s_trace_tab0=NULL
  r_trace_tab0=NULL
  x_trace_tab0=NULL
  s_trace_tab0_C=NULL
  s_trace_tab0_A=NULL
  r_trace_tab0_C=NULL
  r_trace_tab0_A=NULL
  x_trace_tab0_C=NULL
  x_trace_tab0_A=NULL
  
  for(iiM in 1:m.tot){
    load(paste("post_output_multistart/",month.start,"/outputR_",month.start,"_",iiM,".RData",sep=""))
    
    thetatab=cbind(data.frame(thetatab),data.frame(thetaAlltab[,iiH,]))
    theta_inittab=data.frame(theta_initAlltab[,iiH,])
    
    mcmc_samples=length(sim_liktab)
    maxB=sum(sim_liktab!=-Inf)/mcmc_samples
    minB=mcmc.burn*maxB
    picks=c(round(minB*mcmc_samples):round(maxB*mcmc_samples))
    
    thetatabA=rbind(thetatabA,thetatab[picks,])
    theta_inittabA=rbind(theta_inittabA,theta_inittab[picks,])
    
    c_trace_tab0 = rbind(c_trace_tab0,c_trace_tab[picks,iiH,])
    s_trace_tab0 = rbind(s_trace_tab0,s_trace_tab[picks,iiH,])
    r_trace_tab0 = rbind(r_trace_tab0,r_trace_tab[picks,iiH,])
    x_trace_tab0 = rbind(x_trace_tab0,x_trace_tab[picks,iiH,])
    
    if(agestructure==1){
      s_trace_tab0_C=rbind(s_trace_tab0_C,s_trace_tabC[picks,iiH,])
      s_trace_tab0_A=rbind(s_trace_tab0_A,s_trace_tabC[picks,iiH,])
      r_trace_tab0_C=rbind(r_trace_tab0_C,r_trace_tabC[picks,iiH,])
      r_trace_tab0_A=rbind(r_trace_tab0_A,r_trace_tabA[picks,iiH,])
      x_trace_tab0_C=rbind(x_trace_tab0_C,x_trace_tabC[picks,iiH,])
      x_trace_tab0_A=rbind(x_trace_tab0_A,x_trace_tabA[picks,iiH,])
    }
  }
  
  picks=c(1:length(thetatabA[,1]))
  thetatab=thetatabA
  theta_inittab=theta_inittabA
  c_trace_tab=c_trace_tab0
  s_trace_tab=s_trace_tab0
  r_trace_tab=r_trace_tab0
  x_trace_tab=x_trace_tab0
  s_trace_tabC=s_trace_tab0_C
  s_trace_tabA=s_trace_tab0_A
  r_trace_tabC=r_trace_tab0_C
  r_trace_tabA=r_trace_tab0_A
  x_trace_tabC=x_trace_tab0_C
  x_trace_tabA=x_trace_tab0_A
  
  return(list(picks=picks,
              thetatab=thetatab,
              theta_inittab=theta_inittab,
              c_trace_tab=c_trace_tab,
              s_trace_tab=s_trace_tab,
              r_trace_tab=r_trace_tab,
              x_trace_tab=x_trace_tab,
              s_trace_tabC=s_trace_tabC,
              s_trace_tabA=s_trace_tabA,
              r_trace_tabC=r_trace_tabC,
              r_trace_tabA=r_trace_tabA,
              x_trace_tabC=x_trace_tabC,
              x_trace_tabA=x_trace_tabA,
              sim_liktab=sim_liktab,
              accepttab=accepttab
  ))
}
