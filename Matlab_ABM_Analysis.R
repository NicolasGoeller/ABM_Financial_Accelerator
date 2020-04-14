source("Matlab_Modell_Inputs.R")

abm_accelerator <- output_transform()

parameters <- readRDS("parameters.rds")



for (i in names(abm_accelerator)){
  print(i)
  df <- data.frame(abm_accelerator[i][[1]])
  col <- as.vector(df[,"T1"])
  print(class(col))
  print(col[1:10])
  
}
