source("Matlab_Modell_Inputs.R")

output_complete("ABM_Replicator.mat")

dfirms <- readRDS("Data/abm_dfirm_panel.rds")
ufims <- readRDS("Data/abm_ufirm_panel.rds")
banks <- readRDS("Data/abm_bank_panel.rds")

abm_aggregate <- readRDS("Data/abm_aggregate_data.rds")


for (i in names(abm_accelerator)){
  print(i)
  df <- data.frame(abm_accelerator[i][[1]])
  col <- as.vector(df[,"T1"])
  print(class(col))
  print(col[1:10])
  
}

