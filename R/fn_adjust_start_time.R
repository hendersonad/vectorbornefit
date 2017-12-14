##' Function to adjusst start time and data when sampling start times
##' 
##' Funciton that reloads data and fills in zeroes for time series and date series if sampling start time
##' @param sample.start.point True/False whether to sample start time. Defaults to F 
##' @param t0 New sampled value of t0
##' @export 

fn_adjust_start_time <- function(sample.start.point=F, t0){
  if(sample.start.point==T){
    Sample.StartTime <- (log(t0))
    if(Sample.StartTime>0){
      #Sample.StartTime <- Sample.StartTime*-1
      new.start.time <- startdate + (Sample.StartTime*30)
      data <- load.data.multistart(agestructure, add.nulls=0,new.start.time, virusTab[iiH], dataTab[iiH], serology.excel, init.conditions.excel)
    }else if(Sample.StartTime<0){
      new.start.time <- startdate + (Sample.StartTime*30)
      data <- load.data.multistart(agestructure, add.nulls=0,new.start.time, virusTab[iiH], dataTab[iiH], serology.excel, init.conditions.excel)
    }
  }else{
    data <- load.data.multistart(agestructure,add.nulls=0, startdate, virusTab[iiH], dataTab[iiH], serology.excel, init.conditions.excel)
  }
  return(data)
}

