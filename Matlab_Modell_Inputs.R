output_complete <- function() {

# Required packages  
require(R.matlab) 
#require(data.table)
#require(rlist)  

# Import data
abm_output <- readMat("ABM_Replicator.mat")

# Separate parameter values from actual output
keep <- sapply(abm_output, function(i) length(i) > 1)
abm_accelerator <- abm_output[keep]
parameters <- abm_output[!keep]
saveRDS(parameters, "parameters.rds")

## Transform matrix in data.frame format
Rrify <- function(mat){
  agents <- sapply(c(1:ncol(mat)), function(i) paste0("A", as.character(i)))
  mat <- data.frame(mat)
  names(mat) <- agents
  return(mat)
}

# Take out network matrices
networks <- abm_accelerator[names(abm_accelerator) %in% c("BU", "BD", "UD")]
netnames <- names(networks)
abm_accelerator <- abm_accelerator[!(names(abm_accelerator) %in% c("BU", "BD", "UD"))]
abm_names <-names(abm_accelerator)

# Reshape list of matrices
abm_data <- list()
for (i in 1:length(abm_accelerator)) {
  data <- Rrify(abm_accelerator[[i]])
  abm_data[[i]] <- data
} 

names(abm_data) <- abm_names

return(abm_accelerator)

}

output_aggregate <- function(){

  # Required packages  
  require(R.matlab) 
  #require(data.table)
  require(matrixStats)
  
  # Import data
  abm_output <- readMat("ABM_Replicator.mat")
  
  # Separate parameter values from actual output
  keep <- sapply(abm_output, function(i) length(i) > 1)
  abm_accelerator <- abm_output[keep]
  parameters <- abm_output[!keep]
  saveRDS(parameters, "parameters.rds")
  
  # Take out network matrices
  networks <- abm_accelerator[names(abm_accelerator) %in% c("BU", "BD", "UD")]
  netnames <- names(networks)
  abm_accelerator <- abm_accelerator[!(names(abm_accelerator) %in% c("BU", "BD", "UD"))]
  abm_names <-names(abm_accelerator)
  
  # Data Aggregation
  abm_aggregate <- data.frame("T"=1:1000)
  abm_aggregate[,"Yd"] <- rowSums(abm_accelerator[["Yd"]])
  abm_aggregate[,"Yu"] <- rowSums(abm_accelerator[["Qu"]])
  
  abm_aggregate[,"BRd"] <- rowSums(abm_accelerator[["BRd"]])
  abm_aggregate[,"BRu"] <- rowSums(abm_accelerator[["BRu"]])
  abm_aggregate[,"BRb"] <- rowSums(abm_accelerator[["BRb"]])
  
  abm_aggregate[,"Ad"] <- rowSums(abm_accelerator[["Ad"]])
  abm_aggregate[,"Au"] <- rowSums(abm_accelerator[["Au"]])
  abm_aggregate[,"Ab"] <- rowSums(abm_accelerator[["Ab"]])
  
  abm_aggregate[,"BB"] <- rowSums(abm_accelerator[["Bd"]]) + rowSums(abm_accelerator[["Bu"]])
  abm_aggregate[,"BAD"] <- rowSums(abm_accelerator[["Badb"]]) + rowSums(abm_accelerator[["Badu"]])
  
  abm_aggregate[,"Rud"] <- rowMeans(abm_accelerator[["Rud"]])
  abm_aggregate[,"Rbd"] <- rowMeans(abm_accelerator[["Rbd"]])
  abm_aggregate[,"Rbu"] <- rowMeans(abm_accelerator[["Rbu"]])
  
  abm_aggregate[,"Prd"] <- rowMeans(abm_accelerator[["Prd"]])
  abm_aggregate[,"Pru"] <- rowMeans(abm_accelerator[["Pru"]])
  abm_aggregate[,"Prb"] <- rowMeans(abm_accelerator[["Prb"]])
  
  abm_aggregate[,"GR"] <- log(diff(abm_aggregate$Yd))
  
  return(abm_aggregate)
}

