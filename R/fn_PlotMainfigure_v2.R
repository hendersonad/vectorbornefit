#' Plot_figure1
#'
#' @param agestructure Binary indicator variable if model is age structured (1) or not (0) between children and adults. Defaults to NULL in which case no sampling happens 
#' @param btsp Number of iterations from MCMC loop to include
#' @param location Number of locations to loop over. Defaults to 1
#' @param seroposdates Should be same as in MCMC loop
#' @param nplots Number of plots - 2 or 3? Defaults to 2
#' @export

#plot_figure1(agestructure=0,btsp=4000,location=1,seroposdates, nplots=2)

#v2 - adjusts for estimate of t0
plot_figure1_v2 <- function(agestructure, btsp, location=1, seroposdates, nplots=2, sample.start.point=T, startdate=as.Date("2013-10-28")){
  dataP=NULL
  locnnLoop=2
  par(mfrow=c(nplots,length(location)))
  par(mgp=c(2,0.7,0))
  par(mar = c(3,4,1,4))
  labelN=1
  
  # Import multiple chains if available
  load_posterior_1 <- load.posteriors(agestructure=0, iiH, mcmc.burn=0.2)
    list2env(load_posterior_1,globalenv())
  
  for(iiH in location){ # NOTE CURRENTLY JUST 2013/14
    # Load time series dataset - need to initial timeser as global (from main model.R)
      data <- load.data.multistart(agestructure, startdate, virusTab[iiH], dataTab[iiH], serology.excel, init.conditions.excel)
        list2env(data,globalenv())
        tMax <- length(y.vals) 
    
    cvector=matrix(NA,nrow=btsp,ncol=tMax)
    ivector=matrix(NA,nrow=btsp,ncol=tMax)
    svector=matrix(NA,nrow=btsp,ncol=tMax)
    rvector=matrix(NA,nrow=btsp,ncol=tMax)
    xvector=matrix(NA,nrow=btsp,ncol=tMax)
    svectorC=matrix(NA,nrow=btsp,ncol=tMax)
    svectorA=matrix(NA,nrow=btsp,ncol=tMax)
    rvectorC=matrix(NA,nrow=btsp,ncol=tMax)
    rvectorA=matrix(NA,nrow=btsp,ncol=tMax)
    
    for(ii in 1:btsp){
      pick=sample(picks,1)
      cvector[ii,]=ReportC(c_trace_tab[pick,1:tMax],thetatab[pick,'rep'],thetatab[pick,'repvol'])
      ivector[ii,]=c_trace_tab[pick,1:tMax]
      svector[ii,]=s_trace_tab[pick,1:tMax]
      xvector[ii,]=x_trace_tab[pick,1:tMax]
      rvector[ii,]=r_trace_tab[pick,1:tMax]/thetatab[pick,]$npop
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
    ciP1=apply(cvector,2,function(x){quantile(x,0.025, na.rm=T)})
    ciP2=apply(cvector,2,function(x){quantile(x,0.975, na.rm=T)})
    ciP150=apply(cvector,2,function(x){quantile(x,0.25, na.rm=T)})
    ciP250=apply(cvector,2,function(x){quantile(x,0.75, na.rm=T)})
    
    # Proportion recovered
    medP_R=apply(rvector,2,function(x){median(x, na.rm=T)}); 
    medP_R=c(medP_R,rep(tail(medP_R,1),150)) 
    ciP1_R=apply(rvector,2,function(x){quantile(x,0.025, na.rm=T)}); 
      ciP1_R=c(ciP1_R,rep(tail(ciP1_R,1),150)) 
    ciP2_R=apply(rvector,2,function(x){quantile(x,0.975, na.rm=T)}); 
      ciP2_R=c(ciP2_R,rep(tail(ciP2_R,1),150))
    
    # Number of sucetpible mosquitoes
    med_sm=apply(xvector,2,function(x){median(x, na.rm=T)})
    ci_sm1=apply(xvector,2,function(x){quantile(x,0.025, na.rm=T)})
    ci_sm2=apply(xvector,2,function(x){quantile(x,0.975, na.rm=T)})
    
    # Number of infected
    med_Inf=apply(ivector,2,function(x){median(x, na.rm=T)})
    ci_inf1=apply(ivector,2,function(x){quantile(x,0.025, na.rm=T)})
    ci_inf2=apply(ivector,2,function(x){quantile(x,0.975, na.rm=T)})
    ci_inf150=apply(ivector,2,function(x){quantile(x,0.25, na.rm=T)})
    ci_inf250=apply(ivector,2,function(x){quantile(x,0.75, na.rm=T)})
    
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
    ylimmax=1.05*max(y.vals,ciP2)
    plot(date_list_all[xSelect],medP[xSelect],ylim=c(0,ylimmax),xlab="",ylab=ifelse(iiH==1,"confirmed cases","suspected cases"),pch=19,cex=1,col='white',yaxs="i",xlim=xRange)
  
    polygon(c(date_list_all[xSelect],rev(date_list_all[xSelect])),c(ciP150[xSelect],rev(ciP250[xSelect])),lty=0,col=rgb(0,0.3,1,0.3))
    polygon(c(date_list_all[xSelect],rev(date_list_all[xSelect])),c(ciP1[xSelect],rev(ciP2[xSelect])),lty=0,col=rgb(0,0.3,1,0.2))
    lines(date_list_all[xSelect],medP[xSelect],type="l",col=rgb(0,0.3,1),xaxt="n",yaxt="n",xlab="",ylab="")
    
    points(date.vals.plotsero[y.vals!=0],y.vals[y.vals != 0],ylim=c(0,ylimmax),pch=19,cex=1,col='black')
    
    title(adj=0,main=LETTERS[labelN])
    labelN=labelN+1
    
    # PLOT SEROCONVERSION
    par(new=TRUE)
    plot(date_list_all[xSelect2],medP_R[xSelect2],type="l",col=NULL,xaxt="n",yaxt="n",xlab="",ylab="",ylim=c(0,1),xlim=xRange)
    
    polygon(c(date_list_all[xSelect2],rev(date_list_all[xSelect2])),c(ciP1_R[xSelect2],rev(ciP2_R[xSelect2])),lty=0,col=col2a)
    lines(date_list_all[xSelect2],medP_R[xSelect2],lty=2,col=col2,xaxt="n",yaxt="n",xlab="",ylab="")
    
    axis(4,col=col2,col.axis=col2)
    mtext("proportion seropositive", side=4, line=2,col=col2,cex=0.7) # Label for 2nd axis
    
    shiftA=-1
    i=1; lci=NULL;uci=NULL;
    if(agestructure==0){
      for(date in seroposdates){
        binomtest <- binom.test(x=nLUM[i],n=nPOP[i])
        lci[i] <- binomtest$conf.int[1]
        uci[i] <- binomtest$conf.int[2]
  
        points(date-shiftA,nLUM[i]/nPOP[i],pch=1,col="red")
        lines(c(date,date)-shiftA,c(lci[i],uci[i]),pch=19,col="red")
        i <- i+1
      }
      }else if(agestructure==1){
        for(date in seroposdates){
        binomtestC <- binom.test(x=nLUM[i]+nLUM,n=nPOP[i]) 
        binomtestA <- binom.test(x=nLUM[(length(nLUM)/2)+1], n=nPOP[(length(nLUM)/2)+1])
        lci_C[i] <- binomtestC$conf.int[1]
        uci_C[i] <- binomtestC$conf.int[2]
        lci_A[i] <- binomtestA$conf.int[1]
        uci_A[i] <- binomtestA$conf.int[2]
        
        points(date-shiftA,nLUM[i]/nPOP[i],pch=1,col="red")
        lines(c(date,date)-shiftA,c(lci_C[i],uci_C[i]),pch=19,col="red")
        
        points(date-shiftA,nLUM[(length(nLUM)/2)+1]/nPOP[(length(nLUM)/2)+1],pch=1,col="red")
        lines(c(date,date)-shiftA,c(lci_A[i],uci_A[i]),pch=19,col="red")
        i <- i+1
        }
      }
    
    # PLOT seasonality function
    if(locationtab[iiH]=="Central2014" | substr(locationtab[iiH],1,4)=="Zika" ){
      
      extra_date = (tMax+1):max(xSelect2)  #max(date_list_all)  + seq(7,300,7)
      date_listSeason = date_list_all[xSelect2] #c(date_list_all, extra_date )
      
      plotCosR0 = NULL; plotCosRR = NULL
      btstrap = sample(picks,100,replace=T)
      
      for(ii in 1:length(btstrap)){
        b_ii = btstrap[ii]
        
        t.start = 0
        time.V = (1:length(date_listSeason))*7 #time.vals + t.start
        date0 = (as.Date("2013-10-28")-date_listSeason[1]) %>% as.numeric() # Shift back as simulation starts from 2013-10-28
        beta_ii = seasonal_f(time.V, date0, amp=thetatab[b_ii,'beta_v_amp'], mid=thetatab[b_ii,'beta_v_mid'])
        
          s_pick = c(s_trace_tab[b_ii,1:tMax],rep(s_trace_tab[b_ii,tMax],length(extra_date)) )/thetatab$npop[b_ii] 
          x_pick = c(x_trace_tab[b_ii,1:tMax],rep(x_trace_tab[b_ii,tMax],length(extra_date)) ) 
          c_pick = c(c_trace_tab[b_ii,1:tMax],rep(c_trace_tab[b_ii,tMax],length(extra_date)) )/thetatab$npop[b_ii] 
          r_pick = c(r_trace_tab[b_ii,1:tMax],rep(r_trace_tab[b_ii,tMax],length(extra_date)) )/thetatab$npop[b_ii] 
      
        output_rr = calculate_r0(th_in=thetatab[b_ii,],sus_c=s_pick,sus_a=0,sm_c=x_pick,sm_a=0,b_vary=beta_ii)
        
        r0_post = output_rr$r0_out
        rr_post = output_rr$rr_out
        
        # DEBUG CHECK R0
        plotCosR0=rbind(plotCosR0,  r0_post)
        plotCosRR=rbind(plotCosRR,  rr_post)
      }
        c.nume<-function(x){
          bp1=c(median(x),quantile(x,0.025),quantile(x,0.975))
          as.numeric(bp1)}
        plotCosMR0 = apply(plotCosR0,2,c.nume) ;         plotCosMRR = apply(plotCosRR,2,c.nume) 
      #R0
      plot(date_listSeason,plotCosMRR[3,],type="l",col="white",ylim=c(0,4),xlab="",ylab = "reproduction number",yaxs="i",xlim=xRange)
      polygon(c(date_listSeason,rev(date_listSeason)),c(plotCosMR0[2,],rev(plotCosMR0[3,])),lty=0,col=rgb(0,0.3,1,0.2))
      lines(date_listSeason,plotCosMR0[1,],type="l",col=rgb(0,0.3,1),xaxt="n",yaxt="n",xlab="",ylab="")
      
      polygon(c(date_listSeason,rev(date_listSeason)),c(plotCosMRR[2,],rev(plotCosMRR[3,])),lty=0,col=rgb(0,0.6,0.3,0.2))
      lines(date_listSeason,plotCosMRR[1,],type="l",col=rgb(0,0.6,0.3,1),xaxt="n",yaxt="n",xlab="",ylab="")
      
      lines(c(min(date_listSeason),max(date_listSeason)),c(1,1),col="black",lty=2)
      
      title(adj=0,main=LETTERS[labelN]) #paste(locationtabF[iiH],sep=""))
      labelN=labelN+1
  }
  dev.copy(pdf,paste("post_plotsZ/fig1a_modelfit_cases",iiH,".pdf",sep=""),width=6,height=6)#,height=8)
  dev.off()
  
  
  ## PLOT FIGURE 2
  dataP=NULL
  locnnLoop=2
  par(mfrow=c(nplots,length(location)))
  par(mgp=c(2,0.7,0))
  par(mar = c(3,4,1,4))
  labelN=1
 
    ## PLOT INFECTIONS
    ylimmax=1.05*max(y.vals,ci_inf2)
    plot(date_list_all[xSelect],med_Inf[xSelect],ylim=c(0,ylimmax),xlab="",ylab="suspected infections",pch=19,cex=1,col='white',yaxs="i",xlim=xRange)
    
    polygon(c(date_list_all[xSelect],rev(date_list_all[xSelect])),c(ci_inf150[xSelect],rev(ci_inf250[xSelect])),lty=0,col=rgb(0,0.3,1,0.3))
    polygon(c(date_list_all[xSelect],rev(date_list_all[xSelect])),c(ci_inf1[xSelect],rev(ci_inf2[xSelect])),lty=0,col=rgb(0,0.3,1,0.2))
    lines(date_list_all[xSelect],med_Inf[xSelect],type="l",col=rgb(0,0.3,1),xaxt="n",yaxt="n",xlab="",ylab="")
    
    points(date.vals.plotsero[y.vals!=0],y.vals[y.vals != 0],ylim=c(0,ylimmax),pch=19,cex=1,col='black')
    
    title(adj=0,main=LETTERS[labelN])
    labelN=labelN+1
    
    # PLOT SUSEPTIBLE MOSQUITOES
    ylimmax=1.05*max(ci_sm2)
    ylimmin=0.95*min(ci_sm1)
    plot(date_list_all[xSelect],med_Inf[xSelect],ylim=c(ylimmin,ylimmax),xlab="Date",ylab="Rel. Susceptible mosquitoes",pch=19,cex=1,col='white',yaxs="i",xlim=xRange)
    
    polygon(c(date_list_all[xSelect2],rev(date_list_all[xSelect2])),c(ci_sm1[xSelect2],rev(ci_sm2[xSelect2])),lty=0,col=col2a)
    lines(date_list_all[xSelect2],med_sm[xSelect2],lty=2,col=col2,xaxt="n",yaxt="n",xlab="",ylab="")
    
    title(adj=0,main=LETTERS[labelN])
  }#end iiH loop
  dev.copy(pdf,paste("post_plotsZ/fig1b_infectionsplot",iiH,".pdf",sep=""),width=6,height=6)#,height=8)
  dev.off()
}# end function
  

