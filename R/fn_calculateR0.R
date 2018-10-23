#' calculate_r0
#'
#' @param th_in Binary indicator variable if model is age structured (1) or not (0) between children and adults. Defaults to NULL in which case no sampling happens 
#' @param sus_c Suceptible child human hosts
#' @param sus_a Suceptible adult human hosts
#' @param sm_c Suceptible child moisquitoes
#' @param sm_a Suceptible adult moisquitoes
#' @param b_vary Defaults to 1
#' @export

calculate_r0_adam <- function(th_in,sus_c=1,sus_a=1,sm_c=1,sm_a=1,b_vary=1){
  
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
  
  # NEED TO DEBUG THIS
  
  # th_in=thetatab[b_ii,];sus_c=s_pickC[cosPick];sus_a=s_pickA[cosPick];sm_c=x_pickC[cosPick];sm_a=x_pickA[cosPick];b_vary=beta_ii*beta_control
  
  # Rate humans get infected
  # b_cc_h = b_vary * th_in$beta ; b_aa_h = b_cc_h * th_in$beta3
  # b_ca_h = th_in$beta2* b_aa_h ;  b_ac_h=b_ca_h 
  
  # Rate vectors get infected
  # b_cc_v=b_cc_h * th_in$beta_v    ; b_aa_v=b_aa_h * th_in$beta_v
  # b_ca_v=b_ca_h * th_in$beta_v    ; b_ac_v=b_ca_v
  
  # Rate humans get infected -- FORMULATION WITH ONE MOSQUITO POP AND DIFFERENT BITING RATES
  b_hv = b_vary * th_in$beta_h ; b_hh = 0
  
  # Rate vectors get infected
  b_vh = b_hv * th_in$beta_v 
  b_vv = 0
  
  rr_hh=rep(0,length(b_vary)); rr_vv = rr_hh;
  
  rr_hv = (sus_c*b_hv/th_in$MuV)*(th_in$Exp/(th_in$Exp+th_in$MuV)) 
  rr_vc = (sm_c*b_vh/th_in$Inf.)
  
  rr_post = NULL
  for(ii in 1:length(b_vary)){
    rr_post=c(rr_post,(max(Re(eigen(matrix(c(rr_hh[ii],rr_hv[ii],rr_vc[ii],rr_vv[ii]),nrow=2))$values)))  )
  }
  # NEED sqrt ?
  
  #jj=(rr_post==max(rr_post[round(length(rr_post)/2):length(rr_post)]))
  #print( matrix(c(rr_cc[jj],rr_ac[jj],rr_ca[jj],rr_aa[jj]),nrow=2) )
  
  return( list(r0_out=r0_post, rr_out=rr_post,rr_mat = matrix(c(rr_hh[1],rr_hv[1],rr_vc[1],rr_vv[1]),nrow=2,byrow=T)) )

}
