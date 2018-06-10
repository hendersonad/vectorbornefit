#' calculate_r0
#'
#' @param th_in Binary indicator variable if model is age structured (1) or not (0) between children and adults. Defaults to NULL in which case no sampling happens 
#' @param sus_c Suceptible child human hosts
#' @param sus_a Suceptible adult human hosts
#' @param sm_c Suceptible child moisquitoes
#' @param sm_a Suceptible adult moisquitoes
#' @param b_vary Defaults to 1
#' @export

calculate_r0 <- function(th_in,sus_c=1,sus_a=1,sm_c=1,sm_a=1,b_vary=1){
  
  # Rate humans get infected
  # b_cc_h = b_vary * th_in$beta ; b_aa_h = b_cc_h * th_in$beta3
  # b_ca_h = th_in$beta2* b_aa_h ;  b_ac_h=b_ca_h 
  
  # Rate vectors get infected
  # b_cc_v=b_cc_h * th_in$beta_v    ; b_aa_v=b_aa_h * th_in$beta_v
  # b_ca_v=b_ca_h * th_in$beta_v    ; b_ac_v=b_ca_v
  
  # Rate humans get infected -- FORMULATION WITH ONE MOSQUITO POP AND DIFFERENT BITING RATES
  b_cc_h = b_vary * th_in$beta_h    ;   b_aa_h = 0
  b_ca_h = b_cc_h * th_in$beta_h  ;   b_ac_h=0
  
  # Rate vectors get infected
  b_cc_v=b_cc_h * th_in$beta_v    ; b_aa_v=0
  b_ca_v=b_ca_h * th_in$beta_v    ; b_ac_v=0
  
  r0_cc = (b_vary*b_cc_h/th_in$Inf.)*(b_cc_v/th_in$MuV)*(th_in$Exp/(th_in$Exp+th_in$MuV)) + (b_vary*b_ca_h/th_in$Inf.)*(b_ac_v/th_in$MuV)*(th_in$Exp/(th_in$Exp+th_in$MuV))
  r0_ac = (b_vary*b_ac_h/th_in$Inf.)*(b_cc_v/th_in$MuV)*(th_in$Exp/(th_in$Exp+th_in$MuV)) + (b_vary*b_aa_h/th_in$Inf.)*(b_ac_v/th_in$MuV)*(th_in$Exp/(th_in$Exp+th_in$MuV))
  r0_ca = (b_vary*b_ca_h/th_in$Inf.)*(b_aa_v/th_in$MuV)*(th_in$Exp/(th_in$Exp+th_in$MuV)) + (b_vary*b_cc_h/th_in$Inf.)*(b_ca_v/th_in$MuV)*(th_in$Exp/(th_in$Exp+th_in$MuV))
  r0_aa = (b_vary*b_aa_h/th_in$Inf.)*(b_aa_v/th_in$MuV)*(th_in$Exp/(th_in$Exp+th_in$MuV)) + (b_vary*b_ac_h/th_in$Inf.)*(b_ca_v/th_in$MuV)*(th_in$Exp/(th_in$Exp+th_in$MuV))
  
  r0_post = NULL
  for(ii in 1:length(r0_cc)){
    r0_post=c(r0_post,( max(Re(eigen(matrix(c(r0_cc[ii],r0_ac[ii],r0_ca[ii],r0_aa[ii]),nrow=2))$values)))  )
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
