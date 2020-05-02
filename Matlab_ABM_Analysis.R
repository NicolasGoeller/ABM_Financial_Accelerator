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



test <- abm_accelerator["d_networth"][[1]]
test1 <- abm_accelerator["d_production"][[1]]
test2 <- abm_accelerator["d_labour"][[1]]
test3 <- abm_accelerator["d_loans"][[1]]
test5 <- abm_accelerator["d_leverage"][[1]]
test6 <- abm_accelerator["d_profit"][[1]]

bk <-abm_accelerator["d_bankruptcy"][[1]]

test4 <- cbind(test2[,2], test3[,2])
test4$diff <- test4[,1] - test4[,2]


##Check for complex numbers
for (i in names(abm_accelerator)){ 
  print(i) 
  df <- data.frame(abm_accelerator[i][[1]]) 
  col <- as.vector(df[,"T1"]) 
  print(class(col)) 
  print(col[1:10])
}
