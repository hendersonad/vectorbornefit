#' Setup results vectors and dataframes and load initial values
#' 
#' Setup results vectors and dataframes and load initial values
#' @param agestructure Binary indicator variable if model is age structured (1) or not (0) between children and adults. Defaults to NULL in which case no sampling happens 
#' @param iiH PLace in locationtab of location being processed
#' @param parameter_est_file Name of csv file with parameters to be estimated included
#' @export
#parameter_est_file="parameters_est_denv"
results.set.up_SIR <- function(agestructure, iiH, parameter_est_file){
# Set up vectors for storing their values during MCMC loop
  thetaAll=data.frame(rep(NA,locnn))
  for (i in 1:(length(thetaR_IC_local$param)-1)){
    thetaAll <- cbind(thetaAll, rep(NA,locnn))  
  }
  names(thetaAll) <- thetaR_IC_local$param

  compartments=c('s_init','i1_init','r_init')
    
  theta_initAll=data.frame(rep(NA,locnn))
  
  for (i in 1:(length(compartments)-1)){
    theta_initAll <- cbind(theta_initAll, rep(NA,locnn))  
  }
  names(theta_initAll) <- compartments


# Loop to load data for specific location & set initial conditions
for(iiH in itertab){
  ##  Global parameters
  c1=(names(thetaR_IC_local)==locationtab[iiH])
  
  theta <- NULL
  for (i in thetaR_IC_global$param){
    theta <- c(theta, thetaR_IC_global[thetaR_IC_global$param==i, c1])  
  }
  names(theta) <- thetaR_IC_global$param
  
  # Scale age groups by total size
  popsize=theta["npop"]
  popsizeC=theta["npopC"]
  popsizeA=theta["npopA"]
  popsizeTot=(popsizeC+popsizeA); names(popsizeTot)='npop'
  popsizeC=round(popsize* popsizeC/popsizeTot); popsizeA=round(popsize* popsizeA/popsizeTot)
  
  ## Local parameters
  for (i in thetaR_IC_local$param){
    thetaAll[iiH, names(thetaAll)==i] <- thetaR_IC_local[thetaR_IC_local$param==i, c1] 
  }
  
  if(sum(names(thetaAll)=="Exp")>0){
    thetaAll[iiH,"Exp"] <- 1/thetaAll[iiH,"Exp"]
  }
  if(sum(names(thetaAll)=="Inf")>0){
    thetaAll[iiH,"Inf"] <- 1/thetaAll[iiH,"Inf"]
  }
  if(sum(names(thetaAll)=="Vex")>0){
    thetaAll[iiH,"Vex"] <- 1/thetaAll[iiH,"Vex"]
  }
  if(sum(names(thetaAll)=="MuV")>0){
    thetaAll[iiH,"MuV"] <- 1/thetaAll[iiH,"MuV"]
  }
  if(sum(names(thetaAll)=="sigma")>0){
    thetaAll[iiH,"sigma"] <- 1/thetaAll[iiH,"sigma"]
  }
  
  ## Initial compartment conditions
  initial_inf=as.numeric(thetaAll[iiH,'inf0']) #*(popsizeTot/2))
  init_vec=as.numeric(thetaAll[iiH,'vec0'])
  
  # initial recovered 
  theta_initAll[iiH,"r_init"]=0
    theta_initAll[iiH,"i1_init"]=initial_inf
    theta_initAll[iiH,"s_init"]=popsizeTot-theta_initAll[iiH,"i1_init"]-theta_initAll[iiH,"r_init"]
  }

## Covariance matrices 
parameters_est <- read.csv(paste0("data_sets/",parameter_est_file,".csv"), stringsAsFactors = F) 
parms_to_est <- parameters_est$parameters_est
compartments_to_est <- parameters_est$compartments_est

#theta - global
nparam=length(theta) 
npc=rep(0,nparam)
npc[match(parms_to_est,names(theta))]=1
cov_matrix_theta0 = diag(npc)
  colnames(cov_matrix_theta0)=names(theta)
  rownames(cov_matrix_theta0)=names(theta)

#thetaAll - local
nparamA=length(thetaAll[1,]) 
npcA=rep(0,nparamA)
npcA[match(parms_to_est,names(thetaAll[1,]))]=1
cov_matrix_thetaAll = diag(npcA)
  colnames(cov_matrix_thetaAll)=names(thetaAll[1,])
  rownames(cov_matrix_thetaAll)=names(thetaAll[1,])

#theta_initAll 
lmv=length(theta_initAll[1,]) 
npcInit=rep(0,lmv)
npcInit[match(compartments_to_est,names(theta_initAll[1,]))]=1
cov_matrix_theta_initAll = diag(npcInit)
  colnames(cov_matrix_theta_initAll)=names(theta_initAll[1,])
  rownames(cov_matrix_theta_initAll)=names(theta_initAll[1,])

## Set up empty frames to store results
thetatab=matrix(NA,nrow=(MCMC.runs+1),ncol=length(theta))
colnames(thetatab)=names(theta)
thetatab[1,]=theta

thetaAlltab=array(NA, dim=c(MCMC.runs+1,locnn,length(thetaAll[1,])),dimnames=list(NULL,NULL,names(thetaAll)))
thetaAlltab[1,,]=as.matrix(thetaAll)

theta_initAlltab=array(NA, dim=c(MCMC.runs+1,locnn,length(theta_initAll[1,])),dimnames=list(NULL,NULL,names(theta_initAll)))
theta_initAlltab[1,,]=as.matrix(theta_initAll)

prior=rep(1,(MCMC.runs+1))
sim_liktab=rep(-Inf,(MCMC.runs+1))
accepttab=rep(NA,(MCMC.runs))
max.length = length(time.vals)

#total pop
#c
c_trace_tab=array(NA, dim=c(MCMC.runs+1,locnn,max.length)) 
#s
s_trace_tab=array(NA, dim=c(MCMC.runs+1,locnn,max.length)) 
#r
r_trace_tab=array(NA, dim=c(MCMC.runs+1,locnn,max.length)) 
#x
x_trace_tab=array(NA, dim=c(MCMC.runs+1,locnn,max.length)) 
#age split
#c
c_trace_tabC=array(NA, dim=c(MCMC.runs+1,locnn,max.length))
c_trace_tabA=array(NA, dim=c(MCMC.runs+1,locnn,max.length))
#s
s_trace_tabC=array(NA, dim=c(MCMC.runs+1,locnn,max.length))
s_trace_tabA=array(NA, dim=c(MCMC.runs+1,locnn,max.length))
#r
r_trace_tabC=array(NA, dim=c(MCMC.runs+1,locnn,max.length))
r_trace_tabA=array(NA, dim=c(MCMC.runs+1,locnn,max.length))
#x
x_trace_tabC=array(NA, dim=c(MCMC.runs+1,locnn,max.length))
x_trace_tabA=array(NA, dim=c(MCMC.runs+1,locnn,max.length))

return(list(prior=prior,
            sim_liktab=sim_liktab,
            accepttab=accepttab,
            max.length=max.length,
            c_trace_tab=c_trace_tab,
            s_trace_tab=s_trace_tab,
            r_trace_tab=r_trace_tab,
            x_trace_tab=x_trace_tab,
            c_trace_tabC=c_trace_tabC,
            c_trace_tabA=c_trace_tabA,
            s_trace_tabC=s_trace_tabC,
            s_trace_tabA=s_trace_tabA,
            r_trace_tabC=r_trace_tabC,
            r_trace_tabA=r_trace_tabA,
            x_trace_tabC=x_trace_tabC,
            x_trace_tabA=x_trace_tabA,
            thetatab=thetatab,
            thetaAlltab=thetaAlltab,
            theta_initAlltab=theta_initAlltab,
            cov_matrix_theta0=cov_matrix_theta0,
            cov_matrix_thetaAll=cov_matrix_thetaAll,
            cov_matrix_theta_initAll=cov_matrix_theta_initAll))
}
