source("Matlab_Modell_Inputs.R")

output_complete("ABM_Replicator.mat")

orig_file <- "C:/Users/Nicolas/OneDrive - Zeppelin-University gGmbH/Dokumente/Studium/Humboldtprojekt/Humbold_ABM/PCC_orig100.mat"

original_abm_output("PCC100.mat")

### Model outputs of replicator model
dfirms <- readRDS("Data/abm_dfirm_panel.rds")
ufims <- readRDS("Data/abm_ufirm_panel.rds")
banks <- readRDS("Data/abm_bank_panel.rds")
abm_aggregate <- readRDS("Data/abm_aggregate_data.rds")

### Model outputs of original
orig_abm <- readRDS("Data/abm_ori_data_1.rds")
analysis_abm <- readRDS("Data/abm_ori_facts_1.rds")


sum(abm_aggregate$BRd)
hist(abm_aggregate$BDu)

hist(orig_abm$FALLU)

hist(orig_abm$FALLD)
hist(abm_aggregate$BRd)

plot(abm_aggregate$Yd, type= "l")

?plot()

for (i in names(abm_accelerator)){
  print(i)
  df <- data.frame(abm_accelerator[i][[1]])
  col <- as.vector(df[,"T1"])
  print(class(col))
  print(col[1:10])
  
}

