#install.packages("devtools")
library("devtools")
#devtools::install_github("klutometis/roxygen")
library(roxygen2)

setwd("/Users/AHenderson/Documents/R")
#create("vectorbornefit")

# create documentation
setwd("./vectorbornefit")
document()

# install
setwd("..")
install("vectorbornefit")

#test


?MCMCloop_withtimevalsadj
# install_github('vectorbornefit','a-henderson91')  