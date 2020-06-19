# I don't see the point of reporting ABB, ADD, AUU, ERU, and ERD as they are not discussed in the paper
# Possible changes in Matlab
### calculation of FSD, to show the change over time, not only the snapshot of the last round
### Calculation of Network degree distribution
# But still: it is probably not necessary to go that deep for replication standards, more important is own analysis & policy

original_abm_output <- function(matfile){
  
  require(R.matlab)
  abm_original <- readMat(matfile)

  keep <- sapply(abm_original, function(i) length(i) == 1000)
  original_results <- abm_original[keep]
  original_results <- data.frame(do.call(cbind, original_results))
  names(original_results) <- names(abm_original[keep])


  leave <- sapply(abm_original, function(i) (length(i) == 250 | length(i) == 500 | length(i) == 100))
  abm_original <- abm_original[!leave & !keep]

  results <- c("mc" ,"CCC", "CORR.DU.B", "CORR.D.B", "CORR.D.U", "CORR.U.B",
               "PB", "PBB", "PBD", "PBU", "SK", "KR")
  abm_analysis <- abm_original[results]

  saveRDS(original_results, paste0("Data/abm_ori_data_", as.character(abm_analysis[["mc"]][1]), ".rds"))
  saveRDS(abm_analysis, paste0("Data/abm_ori_facts_", as.character(abm_analysis[["mc"]][1]), ".rds"))
}

original_abm_output("PCC100.mat")

orig_abm <- readRDS("Data/abm_ori_data_1.rds")
analysis_abm <- readRDS("Data/abm_ori_facts_1.rds")
