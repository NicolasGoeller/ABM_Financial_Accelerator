library(R.matlab)
library(matrixStats)

abm_output <- readMat("ABM_Replicator.mat")

# Separate parameter values from actual output
keep <- c()
for (i in 1:length(abm_output)) {
  if (length(abm_output[[i]]) > 1){
    keep <- c(keep, i)
  }
}
abm_output <- abm_output[keep]
str(abm_output)

abm_aggregate <- function(matrix){
  out <- data.table()
  means <- rowMeans(matrix)
  stds <- rowSds(matrix)
  out <- cbind(out, means)
  out <- cbind(out, stds)
  return(out)
}

abm_aggregate(abm_output[[1]])

abm_names <- names(abm_output)
#print(abm_names)
abm_agg <- data.table()

for (i in 1:length(abm_output)){
  agg <- abm_aggregate(abm_output[[i]])
  name <- abm_names[i]
  #print(name)
  names(agg) <- c(paste0(name,"_m"), paste0(name,"_sd"))
  abm_agg <- cbind(abm_agg, agg) ###seems like its not working here
}
