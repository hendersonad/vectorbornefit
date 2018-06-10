#' calculate_r0
#'
#' @param th_in Binary indicator variable if model is age structured (1) or not (0) between children and adults. Defaults to NULL in which case no sampling happens 
#' @param sus_c Suceptible child human hosts
#' @param sus_a Suceptible adult human hosts
#' @param sm_c Suceptible child moisquitoes
#' @param sm_a Suceptible adult moisquitoes
#' @param b_vary Defaults to 1
#' @export

calculate_r0_v2 <- function(th_in,sus_c=1,sus_a=1,sm_c=1,sm_a=1,b_vary=1){
  
  # Rate humans get infected -- FORMULATION WITH ONE MOSQUITO POP AND DIFFERENT BITING RATES
  b_hh = 0;
  b_hv = th_in$beta_h
  
  b_vh = th_in$beta_v; 
  b_vv = 0
  
  r0_hh = b_vary*(b_hh/th_in$MuV)*(th_in$Exp/(th_in$Exp+th_in$MuV)) 
  r0_hv = b_vary*(b_hv/th_in$MuV)*(th_in$Exp/(th_in$Exp+th_in$MuV)) 
  r0_vh = b_vary*(b_vh/th_in$Inf.)
  r0_vv = b_vary*(b_vv/th_in$Inf.)
  
  r0_post = NULL
  for(ii in 1:length(r0_hv)){
    r0_post=c(r0_post,( max(Re(eigen(matrix(c(r0_hh[ii],r0_hv[ii],r0_vh[ii],r0_vv[ii]),nrow=2))$values)))  )
  }
  
  rr_cc = (sus_c*b_cc_h/th_in$Inf.)*(sm_c*b_cc_v/th_in$MuV)*(th_in$Exp/(th_in$Exp+th_in$MuV)) + (sus_c*b_ca_h/th_in$Inf.)*(sm_a*b_ac_v/th_in$MuV)*(th_in$Exp/(th_in$Exp+th_in$MuV))
  rr_ac = (sus_a*b_ac_h/th_in$Inf.)*(sm_c*b_cc_v/th_in$MuV)*(th_in$Exp/(th_in$Exp+th_in$MuV)) + (sus_a*b_aa_h/th_in$Inf.)*(sm_a*b_ac_v/th_in$MuV)*(th_in$Exp/(th_in$Exp+th_in$MuV))
  rr_ca = (sus_c*b_ca_h/th_in$Inf.)*(sm_a*b_aa_v/th_in$MuV)*(th_in$Exp/(th_in$Exp+th_in$MuV)) + (sus_c*b_cc_h/th_in$Inf.)*(sm_c*b_ca_v/th_in$MuV)*(th_in$Exp/(th_in$Exp+th_in$MuV))
  rr_aa = (sus_a*b_aa_h/th_in$Inf.)*(sm_a*b_aa_v/th_in$MuV)*(th_in$Exp/(th_in$Exp+th_in$MuV)) + (sus_a*b_ac_h/th_in$Inf.)*(sm_c*b_ca_v/th_in$MuV)*(th_in$Exp/(th_in$Exp+th_in$MuV))
  
  rr_post = NULL
  for(ii in 1:length(r0_cc)){
    rr_post=c(rr_post,( max(Re(eigen(matrix(c(rr_cc[ii],rr_ac[ii],rr_ca[ii],rr_aa[ii]),nrow=2))$values)))  )
  }
  # NEED sqrt ?
  
  #jj=(rr_post==max(rr_post[round(length(rr_post)/2):length(rr_post)]))
  #print( matrix(c(rr_cc[jj],rr_ac[jj],rr_ca[jj],rr_aa[jj]),nrow=2) )
  
  return( list(r0_out=r0_post,rr_out=rr_post) )
  
}
