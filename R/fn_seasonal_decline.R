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