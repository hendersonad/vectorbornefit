#' Plot_figure2
#'
#' @param agestructure Binary indicator variable if model is age structured (1) or not (0) between children and adults. Defaults to NULL in which case no sampling happens 
#' @param btsp Number of iterations from MCMC loop to include
#' @param location Number of locations to loop over. Defaults to 1
#' @param seroposdates Should be same as in MCMC loop
#' @param nplots Number of plots - 2 or 3? Defaults to 2
#' @export

plot_figure2 <- function(agestructure, btsp, location=1, seroposdates, nplots=2){
  dataP=NULL
  locnnLoop=2
  par(mfrow=c(nplots,length(location)))
  par(mgp=c(2,0.7,0))
  par(mar = c(3,4,1,4))
  labelN=1
  for(iiH in location){ # NOTE CURRENTLY JUST 2013/14
    # Load time series dataset - need to initial timeser as global (from main model.R)
    data <- load.data(agestructure, virusTab[iiH], dataTab[iiH], serology.excel, init.conditions.excel)
    list2env(data,globalenv())
    
    # Import multiple chains if available
    load_posterior_1 <- load.posteriors(agestructure=0, iiH, mcmc.burn=0.2)
    list2env(load_posterior_1,globalenv())
    
    tMax <- length(y.vals) #length(y.vals)
    cvector=matrix(NA,nrow=btsp,ncol=tMax)
    ivector=matrix(NA,nrow=btsp,ncol=tMax)
    svector=matrix(NA,nrow=btsp,ncol=tMax)
    rvector=matrix(NA,nrow=btsp,ncol=tMax)
    svectorC=matrix(NA,nrow=btsp,ncol=tMax)
    svectorA=matrix(NA,nrow=btsp,ncol=tMax)
    rvectorC=matrix(NA,nrow=btsp,ncol=tMax)
    rvectorA=matrix(NA,nrow=btsp,ncol=tMax)
    xvector=matrix(NA,nrow=btsp,ncol=tMax)
    
    for(ii in 1:btsp){
      pick=sample(picks,1)
      cvector[ii,]=ReportC(c_trace_tab[pick,1:tMax],thetatab[pick,'rep'],thetatab[pick,'repvol'])
      ivector[ii,]=c_trace_tab[pick,1:tMax]
      svector[ii,]=s_trace_tab[pick,1:tMax]
      rvector[ii,]=r_trace_tab[pick,1:tMax]/thetatab[pick,]$npop
      xvector[ii,]=x_trace_tab[pick,1:tMax]
      if(agestructure==1){
        svectorC[ii,]= s_trace_tabC[pick,1:tMax]/thetatab[pick,]$npopC 
        svectorA[ii,]= s_trace_tabA[pick,1:tMax]/thetatab[pick,]$npopA 
        svector=list(svectorC=svectorC, svectorA=svectorA)
        rvectorC[ii,]= r_trace_tabC[pick,1:tMax]/thetatab[pick,]$npopC 
        rvectorA[ii,]= r_trace_tabA[pick,1:tMax]/thetatab[pick,]$npopA 
        rvector=list(rvectorC=svectorC, rvectorA=svectorA)
      }
      svector[ii,]= s_trace_tab[pick,1:tMax]
    }
    
    # Estimated number of cases 
    medP=apply(cvector,2,function(x){median(x)})
    ciP1=apply(cvector,2,function(x){quantile(x,0.025)})
    ciP2=apply(cvector,2,function(x){quantile(x,0.975)})
    ciP150=apply(cvector,2,function(x){quantile(x,0.25)})
    ciP250=apply(cvector,2,function(x){quantile(x,0.75)})
    
    # Proportion recovered
    medP_R=apply(rvector,2,function(x){median(x)}); 
    medP_R=c(medP_R,rep(tail(medP_R,1),150)) 
    ciP1_R=apply(rvector,2,function(x){quantile(x,0.025)}); 
    ciP1_R=c(ciP1_R,rep(tail(ciP1_R,1),150)) 
    ciP2_R=apply(rvector,2,function(x){quantile(x,0.975)}); 
    ciP2_R=c(ciP2_R,rep(tail(ciP2_R,1),150))
    
    # Number of sucetpible mosquitoes
    med_sm=apply(xvector,2,function(x){median(x)})
    ci_sm1=apply(xvector,2,function(x){quantile(x,0.025)})
    ci_sm2=apply(xvector,2,function(x){quantile(x,0.975)})
    
    # Number of infected
    med_Inf=apply(ivector,2,function(x){median(x)})
    ci_inf1=apply(ivector,2,function(x){quantile(x,0.025)})
    ci_inf2=apply(ivector,2,function(x){quantile(x,0.975)})
    ci_inf150=apply(ivector,2,function(x){quantile(x,0.25)})
    ci_inf250=apply(ivector,2,function(x){quantile(x,0.75)})
    
    # - - - - - - - 
    # PLOT RESULTS
    date_list_all = seq(min(date.vals),max(date.vals)+1000,7) # This expands time.vals to full list
    date_plot=min(date.vals)+ 0 - 7 
    xRange=c(date_plot, max(date.vals))
    date.vals.plotsero=date.vals[1:length(y.vals)]
    
    xSelect=(1:length(date.vals))
    xSelect2= 1:ceiling(as.numeric(max(xRange)-min(xRange))/7)
    
    datacol=rgb(0.4,0.4,0.4)
    col2=rgb(1,0.5,0,1)
    col2a=rgb(1,0.5,0,0.2)
    
    ## PLOT TIMESERIES
    ylimmax=1.05*max(y.vals,ci_inf2)
    plot(date_list_all[xSelect],med_Inf[xSelect],ylim=c(0,ylimmax),xlab="",ylab=ifelse(iiH==1,"confirmed cases","suspected cases"),pch=19,cex=1,col='white',yaxs="i",xlim=xRange)
    
    polygon(c(date_list_all[xSelect],rev(date_list_all[xSelect])),c(ci_inf150[xSelect],rev(ci_inf250[xSelect])),lty=0,col=rgb(0,0.3,1,0.3))
    polygon(c(date_list_all[xSelect],rev(date_list_all[xSelect])),c(ci_inf1[xSelect],rev(ci_inf2[xSelect])),lty=0,col=rgb(0,0.3,1,0.2))
    lines(date_list_all[xSelect],med_Inf[xSelect],type="l",col=rgb(0,0.3,1),xaxt="n",yaxt="n",xlab="",ylab="")
    
    points(date.vals.plotsero[y.vals!=0],y.vals[y.vals != 0],ylim=c(0,ylimmax),pch=19,cex=1,col='black')
    
    title(adj=0,main=LETTERS[labelN])
    labelN=labelN+1
    
    # PLOT SEROCONVERSION
    ylimmax=1.05*max(ci_sm2)
    ylimmin=0.95*min(ci_sm1)
    plot(date_list_all[xSelect],med_Inf[xSelect],ylim=c(ylimmin,ylimmax),xlab="Date",ylab="Rel. Susceptible mosquitoes",pch=19,cex=1,col='white',yaxs="i",xlim=xRange)
    
    polygon(c(date_list_all[xSelect2],rev(date_list_all[xSelect2])),c(ci_sm1[xSelect2],rev(ci_sm2[xSelect2])),lty=0,col=col2a)
    lines(date_list_all[xSelect2],med_sm[xSelect2],lty=2,col=col2,xaxt="n",yaxt="n",xlab="",ylab="")
    
    title(adj=0,main=LETTERS[labelN])
  }#end iiH loop
  dev.copy(pdf,paste("post_plotsZ/Figure_infections_sm_",iiH,".pdf",sep=""),width=6,height=6)#,height=8)
  dev.off()
}# end function


