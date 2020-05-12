source("Matlab_Modell_Inputs.R")

abm_accelerator <- output_complete()

abm_aggregate <- output_aggregate()

parameters <- readRDS("parameters.rds")

str(abm_accelerator)

for (i in names(abm_accelerator)){
  print(i)
  df <- data.frame(abm_accelerator[i][[1]])
  col <- as.vector(df[,"T1"])
  print(class(col))
  print(col[1:10])
  
}



d_net <- abm_accelerator["d_networth"][[1]]
d_pro <- abm_accelerator["d_production"][[1]]
d_lab <- abm_accelerator["d_labour"][[1]]
d_lon <- abm_accelerator["d_loans"][[1]]
d_lev <- abm_accelerator["d_leverage"][[1]]
d_gew <- abm_accelerator["d_profit"][[1]]
d_pri <- abm_accelerator["d_prices"][[1]]

d_bk <-abm_accelerator["d_bankruptcy"][[1]]

test4 <- cbind(test2[,2], test3[,2])
test4$diff <- test4[,1] - test4[,2]
test3 <- d_lab$T3 - d_net$T3
test[test>0]
hist(d_pri$T1)
hist(d_pro$T3)
table(d_bk$T100)

all(d_bk == 0)
all(test == 0)

hist()

hist(d_lab$T3)
hist(d_net$T3)
##Check for complex numbers
for (i in names(abm_accelerator)){ 
  print(i) 
  df <- data.frame(abm_accelerator[i][[1]]) 
  col <- as.vector(df[,"T1"]) 
  print(class(col)) 
  print(col[1:10])
}
