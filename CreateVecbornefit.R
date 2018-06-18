#install.packages("devtools")
library("devtools")
#devtools::install_github("klutometis/roxygen")
library(roxygen2)

setwd("/Users/lsh1510922/Documents/R")
#create("vectorbornefit")

# create documentation
setwd("./vectorbornefit")
document()

# install
setwd("..")
install("vectorbornefit")

?ComputePrior
# install_github('vectorbornefit','a-henderson91')  
