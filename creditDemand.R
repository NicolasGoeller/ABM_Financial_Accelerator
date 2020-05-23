library(R.matlab)


credit_demand <- function(phi, beta, delta, wage){

eqa <- data.frame(matrix(c(1:2000), ncol=1))
names(eqa) <- c("worth")
#eqa$prod <- (eqa$worth^beta)*phi
#eqa$lab <- eqa$prod  * delta * wage
#eqa$loan <- eqa$lab - eqa$worth
eqa$cred <- wage*delta*phi*(eqa$worth^beta) - eqa$worth
#print(all(eqa$loan <= 0))
return(eqa)
}

phi2 <- credit_demand(2,0.9,0.5,1)
phi3 <- credit_demand(3,0.9,0.5,1)
phi2.5 <- credit_demand(2.5,0.9,0.5,1)
phi4 <- credit_demand(4,0.9,0.5,1)

abm_output <- readMat("ABM_Replicator.mat")

credd <- abm_output[["Bd"]]
abm_original <- readMat("PCC100.mat")
or_p <- abm_original[["pd"]]
mod_p <- abm_output[["u"]][1000,]
mean(or_p)
mean(mod_p)
