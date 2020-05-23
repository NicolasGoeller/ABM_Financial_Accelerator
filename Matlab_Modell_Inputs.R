output_complete <- function(matfile) {

# Required packages  
require(R.matlab) 
require(ineq)
require(matrixStats)
#require(data.table)
#require(rlist)  

# Import data
abm_output <- readMat(matfile)

# Separate parameter values from actual output
keep <- sapply(abm_output, function(i) length(i) > 1)
abm_accelerator <- abm_output[keep]
parameters <- abm_output[!keep]
saveRDS(parameters, "Data/parameters.rds")

## Transform matrix in data.frame format
Rrify <- function(mat){
  if (ncol(mat) == 100){agents <- sapply(c(1:ncol(mat)), function(i) paste0("B", as.character(i)))}
  else if (ncol(mat) == 250) {agents <- sapply(c(1:ncol(mat)), function(i) paste0("U", as.character(i)))}
  else if (ncol(mat) == 500) {agents <- sapply(c(1:ncol(mat)), function(i) paste0("D", as.character(i)))}
  
  mat <- data.frame(mat)
  names(mat) <- agents
  mat[,"T"] <- c(1:1000)
  mat <- gather(mat, key= "ID", value = "test", 1:(length(mat)-1))
  return(mat)
}

# Differentiate agent types
banks_var <-  sapply(abm_accelerator, function(i) ncol(i) == 100)
ufirms_var <-  sapply(abm_accelerator, function(i) ncol(i) == 250)
dfirms_var <-  sapply(abm_accelerator, function(i) ncol(i) == 500)

banks <- abm_accelerator[banks_var]
ufirms <- abm_accelerator[ufirms_var]
dfirms <- abm_accelerator[dfirms_var]

bank_names <- names(banks)
ufirm_names <- names(ufirms)
dfirm_names <- names(dfirms)

# Reshape list of matrices for banks
abm_banks <- Rrify(banks[[1]])
names(abm_banks)[3] <- bank_names[1]
for (i in 2:length(banks)) {
  b_data <- Rrify(banks[[i]])
  names(b_data)[3] <- bank_names[i]
  abm_banks <- left_join(abm_banks, b_data, by= c("T", "ID"))
}  

# Reshape list of matrices for u firms
abm_ufirms <- Rrify(ufirms[[1]])#
names(abm_ufirms)[3] <- ufirm_names[1]
for (i in 2:length(ufirms)) {
  u_data <- Rrify(ufirms[[i]])
  names(u_data)[3] <- ufirm_names[i]
  abm_ufirms <- left_join(abm_ufirms, u_data, by= c("T", "ID"))
} 

# Reshape list of matrices for d firms
abm_dfirms <- Rrify(dfirms[[1]])#)
names(abm_dfirms)[3] <- dfirm_names[1]
for (i in 2:length(dfirms)) {
  d_data <- Rrify(dfirms[[i]])
  names(d_data)[3] <- dfirm_names[i]
  abm_dfirms <- left_join(abm_dfirms, d_data, by= c("T", "ID"))
} 

# Save of in files for panel data
saveRDS(abm_banks, "Data/abm_bank_panel.rds")
saveRDS(abm_ufirms, "Data/abm_ufirm_panel.rds")
saveRDS(abm_dfirms, "Data/abm_dfirm_panel.rds")

# Functions for aggregation
rowGinis <- function(mat){
  gini <- c()
  for (i in 1:nrow(mat)){
    gini <- c(gini, Gini(mat[i,]))
  }
  return(gini)
}
  
rowSkews <- function(mat){
  skew <- c()
  for (i in 1:nrow(mat)){
    skew <- c(skew, (mean(mat[i,]) - median(mat[i,]))/ mean(mat[i,] -  median(mat[i,])))
  }
  return(skew)
}
  

link_udis <- function(mat){
  link_mat <- matrix(nrow = 1000, ncol = 250)
  for (j in 1:nrow(mat)){
    links <- c()
    for (i in 1:250){
      links <- c(links, length(mat[j,mat[j,] == i]))
    }
    link_mat[j,] <- links
  }
  return(link_mat)
}
  
link_bdis <- function(mat){
  link_mat <- matrix(nrow = 1000, ncol = 100)
  for (j in 1:nrow(mat)){
    links <- c()
    for (i in 1:100){
      links <- c(links, length(mat[j,mat[j,] == i]))
    }
    link_mat[j,] <- links
  }
  return(link_mat)
}
  
  
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
abm_aggregate[,"BAD"] <- rowSums(abm_accelerator[["BDb"]]) + rowSums(abm_accelerator[["BDu"]])
  
abm_aggregate[,"Rud"] <- rowMeans(abm_accelerator[["Rud"]])
abm_aggregate[,"Rbd"] <- rowMeans(abm_accelerator[["Rbd"]])
abm_aggregate[,"Rbu"] <- rowMeans(abm_accelerator[["Rbu"]])
  
abm_aggregate[,"Prd"] <- rowMeans(abm_accelerator[["PId"]])
abm_aggregate[,"Pru"] <- rowMeans(abm_accelerator[["PIu"]])
abm_aggregate[,"Prb"] <- rowMeans(abm_accelerator[["PIb"]])
  
abm_aggregate[,"GR"] <- c(NA, log(diff(abm_aggregate$Yd)))
  
abm_aggregate[,"dsize_mean"] <- rowMeans(abm_accelerator[["Ad"]])
abm_aggregate[,"dsize_median"] <- rowMedians(abm_accelerator[["Ad"]])
abm_aggregate[,"dsize_sd"] <- rowSds(abm_accelerator[["Ad"]])
abm_aggregate[,"dsize_sk"] <- rowSkews(abm_accelerator[["Ad"]])
abm_aggregate[,"dsize_ineq"] <- rowGinis(abm_accelerator[["Ad"]])
  
abm_aggregate[,"usize_mean"] <- rowMeans(abm_accelerator[["Au"]])
abm_aggregate[,"usize_median"] <- rowMedians(abm_accelerator[["Au"]])
abm_aggregate[,"usize_sd"] <- rowSds(abm_accelerator[["Au"]])
abm_aggregate[,"usize_sk"] <- rowSkews(abm_accelerator[["Au"]])
abm_aggregate[,"usize_ineq"] <- rowGinis(abm_accelerator[["Au"]])
  
abm_aggregate[,"bsize_mean"] <- rowMeans(abm_accelerator[["Ab"]])
abm_aggregate[,"bsize_median"] <- rowMedians(abm_accelerator[["Ab"]])
abm_aggregate[,"bsize_sd"] <- rowSds(abm_accelerator[["Ab"]])
abm_aggregate[,"bsize_sk"] <- rowSkews(abm_accelerator[["Ab"]])
abm_aggregate[,"bsize_ineq"] <- rowGinis(abm_accelerator[["Ab"]])
  
linkdis_u <- link_udis(abm_accelerator[["UD"]])
linkdis_b <- link_bdis(abm_accelerator[["BD"]]) + link_bdis(abm_accelerator[["BU"]])
  
abm_aggregate[,"ulink_sk"] <- rowSkews(linkdis_u)
abm_aggregate[,"ulink_gini"] <- rowGinis(linkdis_u)
abm_aggregate[,"ulink_max"] <- rowMaxs(linkdis_u)
#abm_aggregate[,"ulink_mean"] <- 
#abm_aggregate[,"ulink_median"] <- 
  
abm_aggregate[,"blink_sk"] <- rowSkews(linkdis_b)
abm_aggregate[,"blink_gini"] <- rowGinis(linkdis_b)
abm_aggregate[,"blink_max"] <- rowMaxs(linkdis_b)
#abm_aggregate[,"blink_mean"] <- 
#abm_aggregate[,"blink_median"] <-

#Compute total network measures

#Save off file for aggregates  
saveRDS(abm_aggregate, "Data/abm_aggregate_data.rds")
}
