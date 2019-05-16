#install.packages("devtools")
library("devtools")
#devtools::install_github("klutometis/roxygen")
library(roxygen2)

setwd(paste0("/Users/",Sys.info()["user"],"/Documents/R"))
# create documentation
setwd("./vectorbornefit")
document()
# install
setwd("..")
install("vectorbornefit")
