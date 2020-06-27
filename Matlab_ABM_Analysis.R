source("Matlab_Modell_Inputs.R")

output_complete("ABM_Replicator.mat")

original_abm_output("PCC100.mat")

### Model outputs of replicator model
dfirms <- readRDS("Data/abm_dfirm_panel.rds")
ufims <- readRDS("Data/abm_ufirm_panel.rds")
banks <- readRDS("Data/abm_bank_panel.rds")
abm_aggregate <- readRDS("Data/abm_aggregate_data.rds")

### Model outputs of original
orig_abm <- readRDS("Data/abm_ori_data_1.rds")
analysis_abm <- readRDS("Data/abm_ori_facts_1.rds")


sum(abm_aggregate$BRu)

for (i in names(abm_accelerator)){
  print(i)
  df <- data.frame(abm_accelerator[i][[1]])
  col <- as.vector(df[,"T1"])
  print(class(col))
  print(col[1:10])
  
}

