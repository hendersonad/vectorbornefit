#install.packages("devtools")
library("devtools")
#devtools::install_github("klutometis/roxygen")
library(roxygen2)
getwd()
setwd("../vectorbornefit")
#create("vectorbornefit")

# create documentation
document()

# install
setwd("..")
install("vectorbornefit")

?MCMCloop_withtimevalsadj
# install_github('vectorbornefit','a-henderson91')  