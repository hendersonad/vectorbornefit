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
plot_figure1_v4 <- function(agestructure, btsp, location=1, seroposdates, nplots=2, sample.start.point=T, startdate=as.Date("2013-10-28")){
  require(gridExtra)
  
  # Import multiple chains if available
  labelN=1
  load_posterior_1 <- load.posteriors(agestructure=0, iiH, mcmc.burn=0.2)
    list2env(load_posterior_1,globalenv())
  
  for(iiH in itertab){ # NOTE CURRENTLY JUST 2013/14
    # Load time series dataset - need to initial timeser as global (from main model.R)
      data <- load.data.multiloc.multistart(agestructure, add.nulls = 0, startdate, virusTab[iiH], dataTab[iiH], serology.excel, init.conditions.excel)
        list2env(data,globalenv())
      y.vals <- y.vals.df[,iiH]
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
    ciP1_R=apply(rvector,2,function(x){quantile(x,0.025, na.rm=T)}); 
    ciP2_R=apply(rvector,2,function(x){quantile(x,0.975, na.rm=T)}); 
    
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
    
    y.vals.plot <- y.vals
    y.vals.plot[y.vals.plot==0] <- NA
    
    dataframe.p1 <- data.frame(date.vals=date.vals[1:length(y.vals.plot)], y.vals.plot, medP, ciP1, ciP2, ciP150, ciP250, 
                              medP_R, ciP1_R, ciP2_R)
    
    i=1; lci=NULL;uci=NULL;lci_A=NULL;uci_A=NULL;points=NULL;points_A=NULL; date=NULL
    if(agestructure==0){
      for(date in seroposdates){
        binomtest <- binom.test(x=nLUM[i],n=nPOP[i])
        lci[i] <- binomtest$conf.int[1]
        uci[i] <- binomtest$conf.int[2]
        points[i] <- nLUM[i]/nPOP[i]
        i <- i+1
      }
    dataframe.sero <- data.frame(points,lci,uci, seroposdates)
    }else if(agestructure==1){
      for(date in seroposdates){
        binomtestC <- binom.test(x=nLUM[i]+nLUM,n=nPOP[i]) 
        binomtestA <- binom.test(x=nLUM[(length(nLUM)/2)+1], n=nPOP[(length(nLUM)/2)+1])
        lci[i] <- binomtestC$conf.int[1]
        uci[i] <- binomtestC$conf.int[2]
        lci_A[i] <- binomtestA$conf.int[1]
        uci_A[i] <- binomtestA$conf.int[2]
        points[i] <- (nLUM[i]/nPOP[i])
        points_A[i] <- (nLUM[(length(nLUM)/2)+1]/nPOP[(length(nLUM)/2)+1])
        i <- i+1
      }
    dataframe.sero <- data.frame(points,lci,uci, lci_A, uci_A, points_A, seroposdates)
    }

    dataframe.p2 <- data.frame(date.vals=date.vals[1:length(y.vals.plot)], y.vals.plot, med_Inf, ci_inf1, ci_inf2, ci_inf150, ci_inf250)
    
    (p2 <- ggplot(dataframe.p2) + 
      geom_ribbon(aes(x=date.vals, ymin=ci_inf1, ymax=ci_inf2),col=col1a, fill=col1a, alpha=0.2) +
      geom_ribbon(aes(x=date.vals, ymin=ci_inf150, ymax=ci_inf250),col=col1a, fill=col1a, alpha=0.3) +
      geom_line(aes(x=date.vals, y=med_Inf), col=col1, alpha=1, linetype=1) +
      geom_point(aes(x=date.vals, y=y.vals.plot), col=col2) +
      labs(x = "Year", y = "Suspected infections", title = (LETTERS[labelN])) +
      scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
      theme_zika_fiji())
    
    labelN=labelN+1
    
    (p1 <- ggplot(dataframe.p1) + 
      geom_ribbon(aes(x=date.vals, ymin=ciP1_R, ymax=ciP2_R),col=col1a, fill=col1a) +
      geom_line(aes(x=date.vals, y=medP_R), col=col1, linetype=2) +
      geom_point(data=dataframe.sero, aes(x=seroposdates, y=points), col=col2, pch=1, size=3) +
      geom_linerange(data=dataframe.sero, aes(x=seroposdates, ymax=uci, ymin=lci),col=col2, size=0.6) +
      labs(x = "Year", y = "Proportion seropositive", title = (LETTERS[labelN])) +
      scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
      scale_y_continuous(limits = c(0,0.75), breaks=seq(0,1,0.25)) +
      theme_zika_fiji())
    
    labelN=labelN+1
    
    grid.arrange(p2, p1, ncol=1)
    
    dev.copy(pdf,paste("post_plotsZ/fig1a_modelfit",iiH,"_",run.name,"_",locationtab[iiH],".pdf",sep=""),width=6,height=6)#,height=8)
    dev.off()
    
    # PLOT seasonality function - forecasting forward
    date_listSeason = date.vals
    plotCosR0 = NULL; plotCosRR = NULL
    btstrap = sample(picks,400,replace=T)
      
    for(ii in 1:length(btstrap)){
      b_ii = btstrap[ii]
      
      t.start = 0
      time.V = (1:length(date_listSeason))*7 #time.vals + t.start
      date0 = (startdate-date_listSeason[1]) %>% as.numeric() # Shift back as simulation starts from 2013-10-28
      beta_ii = seasonal_f(time.V, date0, amp=thetatab[b_ii,'beta_v_amp'], mid=thetatab[b_ii,'beta_v_mid'])
      
        s_pick = s_trace_tab[b_ii,1:tMax]/thetatab$npop[b_ii] 
        x_pick = x_trace_tab[b_ii,1:tMax] 
        c_pick = c_trace_tab[b_ii,1:tMax]/thetatab$npop[b_ii] 
        r_pick = r_trace_tab[b_ii,1:tMax]/thetatab$npop[b_ii] 
    
      output_rr = calculate_r0(th_in=thetatab[b_ii,],sus_c=s_pick,sus_a=0,sm_c=x_pick,sm_a=0,b_vary=beta_ii)
      
      r0_post = output_rr$r0_out
      rr_post = output_rr$rr_out
      
      # DEBUG CHECK R0
      plotCosR0=rbind(plotCosR0,  r0_post)
      plotCosRR=rbind(plotCosRR,  rr_post)
    }
      c.nume<-function(x){
        bp1=c(mean(x),quantile(x,0.25),quantile(x,0.75))
        as.numeric(bp1)}
      plotCosMR0 = apply(plotCosR0,2,c.nume) ;         plotCosMRR = apply(plotCosRR,2,c.nume) 
    
    #R0 plot
    dataframe.p3 <- data.frame(date.vals, 
                               medR0=plotCosMR0[1,], lciR0=plotCosMR0[2,], uciR0=plotCosMR0[3,],
                               medRR=plotCosMRR[1,], lciRR=plotCosMRR[2,], uciRR=plotCosMRR[3,])
    (p3 <- ggplot(dataframe.p3) + 
        geom_ribbon(aes(x=date.vals, ymin=lciR0, ymax=uciR0),col=col1a, fill=col1a) +
        geom_line(aes(x=date.vals, y=medR0), col=col1, linetype=1) +
        geom_line(aes(x=date.vals, y=1), col=1, linetype=2) +
        labs(x = "Year", y = "R0", title = (LETTERS[labelN])) +
        scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
        theme_zika_fiji())
    
    labelN=labelN+1
    
    (p4 <- ggplot(dataframe.p3) + 
        geom_ribbon(aes(x=date.vals, ymin=lciRR, ymax=uciRR),col=col2a, fill=col2a) +
        geom_line(aes(x=date.vals, y=medRR), col=col2, linetype=1) +
        geom_line(aes(x=date.vals, y=1), col=1, linetype=2) +
        labs(x = "Year", y = "RR", title = (LETTERS[labelN])) +
        scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
        theme_zika_fiji())
      
    labelN=labelN+1

  grid.arrange(p3, p4, ncol=1)
      
  dev.copy(pdf,paste("post_plotsZ/fig1b_repNo",run.name,"_",locationtab[iiH],".pdf",sep=""),width=6,height=6)#,height=8)
  dev.off()

#Quad plot
grid.arrange(p2, p3, p1, p4, ncol=2)
dev.copy(pdf,paste("post_plotsZ/fig1_",run.name,"_",locationtab[iiH],".pdf",sep=""),width=12,height=8)
dev.off()

  }# end location loop
}# end function
  

