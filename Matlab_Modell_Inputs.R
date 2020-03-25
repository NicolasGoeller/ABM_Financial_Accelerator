library(R.matlab)
library(data.table)
library(tidyverse)

abm_output <- readMat("ABM_Replicator.mat")

# Separate paramter values from actual outpu
keep <- c()
for (i in 1:length(abm_output)) {
  if (length(abm_output[[i]]) > 1){
    keep <- c(keep, i)
  }
}

abm_accelerator <- abm_output[keep]
parameters <- abm_output[-keep]

names(abm_accelerator)

## Transform matrix in data.table format
Rrify <- function(matrix){
  steps <- c()
  for (i in 1:nrow(matrix)){
  steps <-  c(steps, paste0("T", as.character(i)))
  }
  data <- data.table(t(matrix))
  names(data) <- steps
  return(data)
}
## D firms data
d_networth <- Rrify(abm_accelerator["Ad"][[1]])
d_production <- Rrify(abm_accelerator["Yd"][[1]])
d_labour <- Rrify(abm_accelerator["Nd"][[1]])
d_loans <- Rrify(abm_accelerator["Bd"][[1]])
d_leverage <- Rrify(abm_accelerator["Ld"][[1]])
du_interest <- Rrify(abm_accelerator["Rud"][[1]])
db_interest <- Rrify(abm_accelerator["Rbd"][[1]])
d_profit <- Rrify(abm_accelerator["PId"][[1]])

## U firms data
u_networth <- Rrify(abm_accelerator["Au"][[1]])
u_production <- Rrify(abm_accelerator["Qu"][[1]])
u_labour <- Rrify(abm_accelerator["Nu"][[1]])
u_loans <- Rrify(abm_accelerator["Bu"][[1]])
u_leverage <- Rrify(abm_accelerator["Lu"][[1]])
ub_interest <- Rrify(abm_accelerator["Rbu"][[1]])
u_profit <- Rrify(abm_accelerator["PIu"][[1]])

## B banks data
b_networth <- Rrify(abm_accelerator["Ab"][[1]])
b_profit <- Rrify(abm_accelerator["PIb"][[1]])

## Credit networks
ub_network <- Rrify(abm_accelerator["BU"][[1]])
db_network <- Rrify(abm_accelerator["BD"][[1]])
du_network <- Rrify(abm_accelerator["UD"][[1]])

hist(db_network$T1)

(1^(0.9))
