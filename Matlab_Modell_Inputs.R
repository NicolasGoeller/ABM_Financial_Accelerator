output_transform <- function() {

  # Required packages  
require(R.matlab) 
require(data.table)
require(rlist)  

# Import data
abm_output <- readMat("ABM_Replicator.mat")

# Separate parameter values from actual output
keep <- c()
for (i in 1:length(abm_output)) {
  if (length(abm_output[[i]]) > 1){
    keep <- c(keep, i)
  }
}

abm_accelerator <- abm_output[keep]
parameters <- abm_output[-keep]
saveRDS(parameters, "parameters.rds")

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

### Better put each frame as element of list
names_list <- c("d_prices", "d_networth", "d_production", "d_int_goods", "d_labour",
                "d_loans", "d_leverage", "du_interest", "db_interest", "d_profit",
                "u_networth", "u_production", "u_labour", "u_loans", "u_leverage",
                "ub_interest", "u_profit", 
                "b_networth", "b_profit", 
                "ub_network", "db_network", "du_network")

## D firms data: 9 tables
d_prices <- Rrify(abm_accelerator["u"][[1]])
d_networth <- Rrify(abm_accelerator["Ad"][[1]])
d_production <- Rrify(abm_accelerator["Yd"][[1]])
d_int_goods <- Rrify(abm_accelerator["Qd"][[1]])
d_labour <- Rrify(abm_accelerator["Nd"][[1]])
d_loans <- Rrify(abm_accelerator["Bd"][[1]])
d_leverage <- Rrify(abm_accelerator["Ld"][[1]])
du_interest <- Rrify(abm_accelerator["Rud"][[1]])
db_interest <- Rrify(abm_accelerator["Rbd"][[1]])
d_profit <- Rrify(abm_accelerator["PId"][[1]])

## U firms data: 7 tables
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

# Zip outputs back together
abm_accelerator <- list(d_prices, d_networth, d_production, d_int_goods, d_labour, d_loans, d_leverage, du_interest, db_interest, d_profit,
                        u_networth, u_production, u_labour, u_loans, u_leverage, ub_interest, u_profit,
                        b_networth, b_profit,
                        ub_network, db_network, du_network)

names(abm_accelerator) <- names_list
return(abm_accelerator)

}
