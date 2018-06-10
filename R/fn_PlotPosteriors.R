#' Plot_Sposteriorss
#' 
#' Plots of prior and posterior densities, ESS and correlation between parameter posteriors
#' @param locnnloop number of locations to loop over. Defaults to 1
#' @param max_season Defaults to 2
#' @param plot.split Max number of params in one correlation plots
#' @param param.names vector of names of parameters to plot
#' @param param.labels vector of labels to put on corr plot
#' @export

plot_Sposteriors<-function(locnnLoop=1, max_season = 2, plot.split=6, param.names, param.labels){
  paramA=NULL
  theta_final <- NA
  
  c.text <- function(x,sigF=3){
    bp1=signif(c(median(x),quantile(x,0.025),quantile(x,0.975)),sigF)
    paste(bp1[1]," (",bp1[2],"-",bp1[3],")",sep="")
  }
  
  for(iiH in 1:locnnLoop){
    
    data <- load.data(agestructure, virusTab[iiH], dataTab[iiH], serology.excel, init.conditions.excel)
      list2env(data,globalenv())
    
    # Import multiple chains if available
    load_posterior_1 <- load.posteriors(agestructure=0, iiH, mcmc.burn=0.2)
      list2env(load_posterior_1,globalenv())
      
    beta_h1 <-  thetatab$beta_h 
    beta_h3 <-  thetatab$beta_h_3 * beta_h1
    beta_h2 <-  thetatab$beta_h_2 * beta_h3 
    beta_v1 <-  thetatab$beta_v * beta_h1 
    beta_v2 <-  thetatab$beta_v * beta_h2 
    beta_v3 <-  thetatab$beta_v * beta_h3 
    Nsize <-    thetatab$npop
    NsizeC <-   thetatab$npopC
    NsizeA <-   thetatab$npopA
    delta_v  <- thetatab$MuV
    alpha_v <-  thetatab$Vex
    alpha_h <-  thetatab$Exp
    gamma <-    thetatab$Inf.
    imm0 <-     thetatab$imm0
    rep <-      thetatab$rep 
    repvol <-   thetatab$repvol 
    t0 <-       thetatab$t0
    
    amp <-      thetatab$beta_v_amp
    mid <-      thetatab$beta_v_mid
    
    inf0 <-     theta_inittab$i1_init
    
    r0_vtab = (beta_h1[picks]/gamma[picks])*
                (beta_v1[picks]/delta_v[picks])*
                (alpha_v[picks]/(alpha_v[picks]+delta_v[picks]))*max_season
    
    par(mfrow=c(3,4))
    par(mar = c(5,5,1,1))
    
    # if(iiH==1){
    colW='grey' 
    colM=rgb(0.2,0.4,1)
    colB='white'
    baW=0.3
    breaks0=seq(0,50,1)
    brekN=15
    
    hist(1/alpha_v[picks],xlab=expression("latent period (V)"),main=NULL,prob=TRUE,border=colB,col=colW,xlim=c(0,30))
    curve(priorVEx(x), col="red", lwd=2, add=TRUE, yaxt="n")
      theta_final[["Vex"]] <- 1/median(alpha_v[picks])
    
    hist(1/delta_v[picks],xlab=expression("lifespan (V)"),main=NULL,prob=TRUE,border=colB,col=colW,xlim=c(0,20))
    curve(priorMuV(x), col="red", lwd=2, add=TRUE, yaxt="n")
      theta_final[["MuV"]] <- 1/median(delta_v[picks])
    
    hist(1/alpha_h[picks],xlab=expression("latent period (H)"),main=NULL,prob=TRUE,border=colB,col=colW,xlim=c(0,20))
    curve(priorExp(x), col="red", lwd=2, add=TRUE, yaxt="n")
      theta_final[["Exp"]] <- 1/median(alpha_h[picks])
    
    hist(1/gamma[picks],xlab=expression("infectious period (H)"),main=NULL,prob=TRUE,border=colB,col=colW,xlim=c(0,20))
    curve(priorInf(x), col="red", lwd=2, add=TRUE, yaxt="n")
      theta_final[["Inf"]] <- 1/median(gamma[picks])

    hist(beta_h1[picks],xlab=expression("beta"),prob=TRUE,main=NULL,border=colB,col=colW)#,xlim=c(0,30))
      theta_final[["beta_h"]] <- median(beta_h1[picks])
      
    hist(amp[picks],xlab=expression("Seasonal amplitude: beta"[amp]),main=NULL,prob=TRUE,border=colB,col=colW, xlim=c(0.2,0.7))
    curve(priorBeta_amp(x), col="red", lwd=2, add=TRUE, yaxt="n")
      theta_final[["beta_v_amp"]] <- median(amp[picks])
      
    hist(mid[picks],breaks=breaks0,xlab=expression("Seasonal midpoint: beta"[mid]),main=NULL,prob=TRUE,border=colB,col=colW,xlim=c(0,6))
    curve(priorBeta_mid(x), col="red", lwd=2, add=TRUE, yaxt="n")
      theta_final[["beta_v_mid"]] <- 1/median(mid[picks])
      
    hist(beta_v1[picks],xlab=expression("beta_v"),prob=TRUE,main=NULL,border=colB,col=colW)#,xlim=c(0,30))
      theta_final[["beta_v"]] <- median(beta_v1[picks])
    
    hist(rep[picks],xlab="proportion of cases reported",main=NULL,border=colB,col=colW,prob=TRUE)
      theta_final[["rep"]] <- median(rep[picks])
    
    hist(repvol[picks],xlab="reporting dispersion",main=NULL,border=colB,col=colW,prob=TRUE)
      theta_final[["repvol"]] <- median(repvol[picks])
    
    #hist(inf0[picks],xlab=expression('I'[0]^h_C),main=NULL,border=colB,col=colW,prob=TRUE)
    #  theta_final[["i1_init"]] <- median(inf0[picks])
    
    hist(imm0[picks],xlab=expression('Initial Immune'),main=NULL,border=colB,col=colW,prob=TRUE)
      theta_final[["imm0"]] <- median(imm0[picks])
    
    hist(log(t0)[picks],xlab=expression('Start point (weeks)'),main=NULL,border=colB,col=colW,prob=TRUE)
      theta_final[["t0"]] <- median(t0[picks])
      
    dev.copy(pdf,paste("post_plotsZ/PlotPosteriors_",iiH,".pdf",sep=""),width=15,height=8) #,locationtab[iiH],
    dev.off()
    
    # PLOT Parameter correlation
    par(mfcol=c(min(plot.split,length(param.names)),min(plot.split,length(param.names))))
    par(mar = c(2,2,1,1),mgp=c(1.8,0.5,0))
    
    thetatab0 = thetatab %>% data.frame()
    thinner.theta=thetatab0[sample(length(thetatab0$beta_h),1000,replace=T),]
    
    for(ii in 1:min(plot.split,length(param.names))){
      for(jj in 1:min(plot.split,length(param.names))){
        if(ii<=jj){
          if(ii == jj){
            hist(thetatab0[[param.names[ii]]],xlab=param.labels[ii],main=paste("ESS=",round(effectiveSize(thetatab0[[param.names[ii]]]))))
          }else{
            plot(x=1,y=1,pch=19,cex=0.2, xlab="", ylab="",col="white",xaxt="n",yaxt="n",axes=F)
          }
        }else{
          plot(thinner.theta[[param.names[ii]]],thinner.theta[[param.names[jj]]],pch=19,cex=0.3, xlab=param.labels[ii], ylab=param.labels[jj])
        }
      }
    }
    dev.copy(png,"post_plotsZ/CorrelationPlot_yr_1.png",units="cm",width=25,height=25,res=200)
    dev.off()
    
    if(length(param.names)>plot.split){
      
    par(mfcol=c(length(param.names)-plot.split,length(param.names)-plot.split))
    par(mar = c(3,3,1,1),mgp=c(1.8,0.5,0))
    
    thetatab0 = thetatab %>% data.frame()
    thinner.theta=thetatab0[sample(length(thetatab0$beta_h),1000,replace=T),]
    
    for(ii in (plot.split+1):length(param.names)){
      for(jj in (plot.split+1):length(param.names)){
        if(ii<=jj){
          if(ii == jj){
            hist(thetatab0[[param.names[ii]]],xlab=param.labels[ii],main=paste("ESS=",round(effectiveSize(thetatab0[[param.names[ii]]]))))
          }else{
            plot(x=1,y=1,pch=19,cex=0.2, xlab="", ylab="",col="white",xaxt="n",yaxt="n",axes=F)
          }
        }else{
          plot(thinner.theta[[param.names[ii]]],thinner.theta[[param.names[jj]]],pch=19,cex=0.3, xlab=param.labels[ii], ylab=param.labels[jj])
        }
        
      }
    }
    
    dev.copy(png,"post_plotsZ/CorrelationPlot_yr_2.png",units="cm",width=25,height=25,res=200)
    dev.off()
    }
    
    # Include bootstrap reporting uncertainty - denominator is total people infected
    wks=length(time.vals)
    
    repBTS=sapply(sample(picks,1000,replace=T),
                  function(x){rr.bt=(Nsize[x]-s_trace_tab[x,wks]); 
                  rnbinom(1,mu=(rr.bt*rep[x]),size=1/repvol[x])/rr.bt })
    
    param1=cbind(
      c.text(r0_vtab,2),
      c.text(100*repBTS,2),
      c.text(100*(1-s_trace_tab[picks,wks]/Nsize[picks]),2),
      c.text(100*(r_trace_tab[picks,1]/Nsize[picks]),2),
      c.text(beta_h1[picks],2),
      c.text(beta_v1[picks],2),
      c.text(1/alpha_h[picks],2),
      c.text(1/gamma[picks],2),
      c.text(1/alpha_v[picks],2),
      c.text(1/delta_v[picks],2),
      c.text(rep[picks],2),
      c.text(repvol[picks],2),
      c.text(inf0[picks],2),
      c.text(amp[picks],2),
      c.text(mid[picks],2),
      c.text(t0[picks],2)
    )
    
    rownames(param1)=c(locationtab[iiH])
    colnames(param1)=c("R0","propn reported (%)","final size","prop_imm","beta_h","beta_v","alpha_h","gamma","alpha_v","delta","r","phi","I_H(0)","beta_amp","beta_mid","t0")
    
    paramA=rbind(paramA,param1)
    
  }
  
  write.csv(t(paramA),paste("post_plotsZ/param1_",iiH,".csv",sep=""))
  
}
