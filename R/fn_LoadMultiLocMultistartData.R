#' Load data for model fitting 
#' 
#' This function loads surveillance cases, serology results and initial parameter values for vector-borne model fitting
#' @param agestructure Binary indicator variable if model is age structured (1) or not (0) between children and adults. Defaults to NULL in which case no sampling happens
#' @param add.nulls Number of time series values to add to end (default is 0)
#' @param startdate Start month of time series
#' @param Virus Name of virus for serology results. e.g. 'Zika'
#' @param cases.file.name Name of csv file with surveillance data
#' @param serology.file.name Name of csv file with serological data 
#' @param init.values.file.name Name of csv file with initial parameter values
#' @export

#agestructure = 0
#startdate  = new.start.time
#Virus  = virusTab[iiH] 
#cases.file.name  = dataTab[iiH]
#serology.file.name  = serology.excel 
#init.values.file.name = init.conditions.excel 
#add.nulls=0
#load.data(agestructure=1, 'Zika', "Fiji_2016Z_timeseries", "Fiji_serology", "thetaR_IC_zika_14")
load.data.multiloc.multistart <- function(agestructure=NULL, add.nulls=0, startdate, Virus, cases.file.name, serology.file.name, init.values.file.name){
  ## Load case data. First column must be 'date'. Column with data must be under same name as in locationtab
  y.vals.df=NULL
  timeser <- read.csv(paste0("data_sets/",cases.file.name,".csv"), stringsAsFactors = F)  # Load ZIKA  data
  for(iiH in itertab){
    colID=(names(timeser)==locationtab[iiH])
    endwk=length(timeser[,1])
    startdate=min(as.Date(timeser$date,origin="1970-01-01"), startdate)
    enddate=max(as.Date(timeser$date,origin="1970-01-01")+add.nulls)
    pickdate=(as.Date(timeser$date,origin="1970-01-01")>=startdate & as.Date(timeser$date,origin="1970-01-01")<=enddate)
    y.vals.null <- rep(0, round((min(as.Date(timeser$date,origin="1970-01-01")) - startdate)/7))
    y.vals.add <- rep(0, (add.nulls/7))
    y.vals = c(y.vals.null, as.numeric(timeser[,colID])[pickdate], y.vals.add)
    endwk=endwk+length(y.vals.null)
    firstentry=min(c(1:endwk)[])
    #pickdate=c(firstentry:endwk)
    y.vals=y.vals[firstentry:endwk]
    if(length(y.vals.null)==0){
      date.vals.null=NULL
    }else{
      date.vals.null=seq(startdate, startdate+((length(y.vals.null)-1)*7), 7)
    }
    if(add.nulls==0){
      date.vals.add=NULL
    }else{
      date.vals.add =seq(round(max(as.Date(timeser$date, origin="1970-01-01"))), round(enddate), 7)
    }
    ## adjust date vals by difference between new.start.time and original specified
    date.vals <- as.Date(timeser$date)
    adj <- length(y.vals.null) - (min(as.Date(timeser$date,origin="1970-01-01")) - startdate)/7
    date.vals <- date.vals + (adj*7)
    date_list=c(date.vals.null, date.vals[pickdate], date.vals.add)
    time.vals = as.numeric(date_list-min(date_list)) + (firstentry-1) * 7 # Shift by start date so both time series line up -- IMPORTANT FOR SEASONALITY
    time.vals = time.vals+time.vals[2]-time.vals[1]
    #time.vals = c(time.vals,seq(max(time.vals)+7,max(time.vals)+7*add.null.dates,7)) # Add extra onto end if needed
    date.vals = as.Date(time.vals, origin = startdate) 
    
    y.vals.df <- cbind(y.vals.df, y.vals)
  }
  colnames(y.vals.df) <- locationtab
  rownames(y.vals.df) <- time.vals
  
  ## Load serology data
  serology <- read.csv(paste0("data_sets/",serology.file.name,".csv"), stringsAsFactors = F)  # Load ZIKA  data
  nLUM=NULL; nPOP=NULL
  for(loc in locationtab){
  if(agestructure==1){
    nLUM <- c(serology$sero2013[serology$virus==Virus & serology$age == 'c' & serology$serology == 'pos'],
              serology$sero2015[serology$virus==Virus & serology$age == 'c' & serology$serology == 'pos'],
              serology$sero2017[serology$virus==Virus & serology$age == 'c' & serology$serology == 'pos'],
                 serology$sero2013[serology$virus==Virus & serology$age == 'a' & serology$serology == 'pos'],
                 serology$sero2015[serology$virus==Virus & serology$age == 'a' & serology$serology == 'pos'],
                 serology$sero2017[serology$virus==Virus & serology$age == 'a' & serology$serology == 'pos']
                )
    nPOP <- c(serology$sero2013[serology$virus==Virus & serology$age == 'c' & serology$serology == 'total'],
              serology$sero2015[serology$virus==Virus & serology$age == 'c' & serology$serology == 'total'],
              serology$sero2017[serology$virus==Virus & serology$age == 'c' & serology$serology == 'total'],
                 serology$sero2013[serology$virus==Virus & serology$age == 'a' & serology$serology == 'total'],
                 serology$sero2015[serology$virus==Virus & serology$age == 'a' & serology$serology == 'total'],
                 serology$sero2017[serology$virus==Virus & serology$age == 'a' & serology$serology == 'total']
                )
  }else{
    nLUM.loc <- c(serology$sero2013pos[serology$virus==Virus & serology$age == 'c' & serology$location == loc]+
                  serology$sero2013pos[serology$virus==Virus & serology$age == 'a' & serology$location == loc],
                  serology$sero2015pos[serology$virus==Virus & serology$age == 'c' & serology$location == loc]+
                  serology$sero2015pos[serology$virus==Virus & serology$age == 'a' & serology$location == loc],
                  serology$sero2017pos[serology$virus==Virus & serology$age == 'c' & serology$location == loc]+
                  serology$sero2017pos[serology$virus==Virus & serology$age == 'a' & serology$location == loc]  
                )
    nLUM <- rbind(nLUM, nLUM.loc)
    nPOP.loc <- c(serology$sero2013total[serology$virus==Virus & serology$age == 'c' & serology$location == loc]+
              serology$sero2013total[serology$virus==Virus & serology$age == 'a' & serology$location == loc],
                  serology$sero2015total[serology$virus==Virus & serology$age == 'c' & serology$location == loc]+
                  serology$sero2015total[serology$virus==Virus & serology$age == 'a' & serology$location == loc],
                  serology$sero2017total[serology$virus==Virus & serology$age == 'c' & serology$location == loc]+ 
                  serology$sero2017total[serology$virus==Virus & serology$age == 'a' & serology$location == loc]  
                )
    nPOP <- rbind(nPOP, nPOP.loc)
    
  }
  }
  row.names(nLUM)=locationtab
  row.names(nPOP)=locationtab
  
  ## Load initial parameter values
  thetaR_IC_global <- read.csv(paste0("data_sets/",init.values.file.name,"_global.csv"),stringsAsFactors=FALSE)
  thetaR_IC_local <- read.csv(paste0("data_sets/",init.values.file.name,"_local.csv"),stringsAsFactors=FALSE)
  
  ## Load vector control dates
if(vector.control==T){ 
  vector_control_dates <- read.csv(paste0("data_sets/", vec.control.excel, ".csv"), stringsAsFactors = F)
  v.c.vals.df <- matrix(time.vals, nrow=length(time.vals), ncol=1)
  for(iiH in itertab){
    colID=(names(vector_control_dates)==locationtab[iiH])
    endwk=length(vector_control_dates[,1])
    startdate=min(as.Date(vector_control_dates$date,origin="1970-01-01"), startdate)
    enddate=max(as.Date(vector_control_dates$date,origin="1970-01-01")+add.nulls)
    pickdate=(as.Date(vector_control_dates$date,origin="1970-01-01")>=startdate & as.Date(vector_control_dates$date,origin="1970-01-01")<=enddate)
    v.c.vals.null <- rep(0, round((min(as.Date(vector_control_dates$date,origin="1970-01-01")) - startdate)/7))
    v.c.vals.add <- rep(0, (add.nulls/7))
    v.c.vals = c(v.c.vals.null, as.numeric(vector_control_dates[,colID])[pickdate], v.c.vals.add)
    
    v.c.vals.df <- cbind(v.c.vals.df, v.c.vals)
  }
  colnames(v.c.vals.df) <- c('time.vals', locationtab)
  }else{
  v.c.vals.df <- data.frame(NULL)
}
  ## 
  return(list(y.vals.df=y.vals.df, v.c.vals.df=v.c.vals.df, time.vals=time.vals, date.vals=date.vals, nLUM=nLUM, nPOP=nPOP, 
              thetaR_IC_global=thetaR_IC_global, thetaR_IC_local=thetaR_IC_local))
}



