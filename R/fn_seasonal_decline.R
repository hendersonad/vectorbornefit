#' Seasonal transmission rate function
#' 
#' Calculates the seasonal tranmsission rate
#' @param time t in time.vals of calc
#' @param date0 reference index value for start of transmission. Defaults to 0
#' @param amp amplitude of seasonal wave. Defaults to 0
#' @param mid midpoint of seasonal wave. Defaults to -0.048 to fix at Feb
#' @keywords Seasonality
#' @export
#' @examples 
#' seasonal_f(date0=0, amp=theta[["beta_v_amp"]], mid=theta[["beta_v_mid"]])

seasonal_f <- function(time, date0=0, amp=0, mid=pi*(3/4)){
  (1 + amp *sin(((time-date0)/365 + mid)*2*pi))
}

#' Flexible decline function transmission rate function
#' 
#' Sigmoid function to allow drop in transmission rate
#' @param time t in time.vals of calc
#' @param date0 reference index value for start of transmission. Defaults to 0
#' @param mask Strength of reduction in transmission rate. Defaults to 0 (i.e. no effect)
#' @param base Height of starting point of function. Defaults to 1
#' @param grad Strength of decline in transmission once begun. Defaults to 0.5
#' @param mid Point at which decline can begin. Defaults to 0
#' @keywords Sigmoid Decline
#' @export
#' @examples 
#' decline_f(mask=1, grad=0.75)

decline_f <- function(time,date0=0,mask=0,base=1,grad=0.5,mid=0){
  (1 - mask*base/(1+exp(-10*grad*((time+date0)/365-mid))))  
}

#'@export
death_f <- function(time,base=1){
  if(time>max(time.vals)){
    0
  }else{
  nearest.week <- max(1,min(which(time.vals>=time)))
  (1 + tanh(time.interventions[nearest.week]*base))
  }
}

#' @export
death_f_linear <- function(time,base=0.5){
#  max(0,1+sum((apply(as.matrix(index.interventions), 1, function(x){((mask*(base*time.interventions[x]))/(1 + exp(-grad*(time-time.vals[x])))) - ((mask*(base*time.interventions[x]))/(1 + exp(-grad*(time-time.vals[x-1]))))}))))
  nearest.week <- max(1,min(which(time.vals>=time)))
  (1 - (time.interventions[nearest.week]*base))  
}

#' Reporting cases function
#' 
#' Given a number of infected hosts, this function estimates the resulting number of cases
#' @param cases List of cases from model
#' @param rep Reporting rate used to generate those cases
#' @param repvol Overdispersion parameter used in Likelihood to estimate those cases
#' @export

ReportC<-function(cases, rep, repvol){
  mu00=cases
  mu01=sapply(mu00,function(x){max(x,0)})
  sapply(mu01,function(x){rnbinom(1, mu=rep*x,size=1/repvol)})
}

#' Cross immunity
#' 
#' Given a level of dengue infection, susceptible Zika hosts are depressed
#' @param time
#' @param immune 
#' @param pop 
#' @export

CrossImmune <- function(time,date0=0,mask=1,base=0.2,grad=2,mid=0.3,wane=1){
  imm <- mask*base/(1+exp(-10*grad*((time+date0)/365-mid))) 
  wane <- (mask*base/(1+exp(-10*grad*((time+date0)/365-wane)))) 
  return(list(imm=imm,wane=wane))
}

#' Flexible decline function transmission rate function
#' 
#' Sigmoid function to allow drop in transmission rate
#' @param time t in time.vals of calc
#' @param date0 reference index value for start of transmission. Defaults to 0
#' @param mask Strength of reduction in transmission rate. Defaults to 0 (i.e. no effect)
#' @param base Height of starting point of function. Defaults to 1
#' @param grad Strength of decline in transmission once begun. Defaults to 0.5
#' @param mid Point at which decline can begin. Defaults to 0
#' @param date.interventions Vector of dates when vector control interventions took place
#' @keywords Sigmoid Decline
#' @export

#deprecated:
#death_f <- function(time,date0=0,mask=0,base=1,grad=0.5,mid=0,date.interventions=NA){
#  no.interventions <- length(date.interventions)
#  death.rate.sum <- 0
#  for(int in 1:no.interventions){
#    time.interventions <- min(which(date.vals>date.interventions[int]))
#    death.rate <- mask*base/(1+exp(-10*grad*((time+365.25-time.vals[time.interventions])/365-mid))) - mask*base/(1+exp(-10*grad*((time+365.25-time.vals[time.interventions+1])/365-mid))) 
#    death.rate.sum <- death.rate.sum+death.rate
#  }
#  death.rate.sum
#}

death_f_log <- function(time,date0=0,mask=1,base=0.5,grad=1e6,mid=1){
  1+sum((apply(as.matrix(index.interventions), 1, function(x){((mask*(base/(1+exp(-time.interventions[x]))))/(1 + exp(-grad*(time-time.vals[x])))) - ((mask*(base/(1+exp(-time.interventions[x]))))/(1 + exp(-grad*(time-time.vals[x-1]))))})))
}

#' @export
death_f_simple <- function(time,date0=0,mask=1,base=0.5,grad=1e6,mid=1){
  1+sum((apply(as.matrix(index.interventions), 1, function(x){((mask*base)/(1 + exp(-grad*(time-time.vals[x])))) - ((mask*base)/(1 + exp(-grad*(time-time.vals[x-1]))))})))
}


#' @export
death_f_basic <- function(time, base){
  nearest.week <- max(1,which(time.vals<time+3.5 & time.vals>time-3.5))
  1-time.interventions[nearest.week]*base
}

#' @export
death_f_logistic <- function(time, base){
  nearest.week <- max(1,which(time.vals<time+3.5 & time.vals>time-3.5))
  1-(base/(1+exp(-time.interventions[nearest.week])))
}
# time.interventions <- v.c.vals.df[,locationtab[iiH]]
# (1-apply(as.matrix(time.vals), 1, function(x){death_f(x)}))
#1 - death_f(time.vals)

#(t1 <- apply(as.matrix(time.vals),1,function(x){death_f(x, base=0.5)}))
#plot(t1, type='l')
#par(new=T)
#plot(time.interventions,col=2, xaxt='n', yaxt='n')
#
#x=1:10
#plot(0.9/(1+exp(-x)), type='l')
