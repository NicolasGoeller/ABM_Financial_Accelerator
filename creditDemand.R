
credit_demand <- function(phi, beta, delta, wage){

eqa <- data.frame(matrix(c(1:2000), ncol=1))
names(eqa) <- c("worth")
eqa$prod <- (eqa$worth^beta)*phi
eqa$lab <- eqa$prod  * delta * wage
eqa$loan <- eqa$lab - eqa$worth
#eqa$cred <- 1*0.5*2*(eqa$worth^0.9) - eqa$worth
print(all(eqa$loan <= 0))
return(eqa)
}

phi2 <- credit_demand(2,0.9,0.5,1)
phi3 <- credit_demand(3,0.9,0.5,1)
phi2.5 <- credit_demand(2.5,0.9,0.5,1)
