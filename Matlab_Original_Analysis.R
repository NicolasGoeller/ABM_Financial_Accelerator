library(R.matlab)

abm_original <- readMat("PCC100.mat")

keep <- c("YD", "YU", "FALLD", "FALLU", "AD", "AU", "AB", "RBD", "RBU",
          "RU", "BAD", "FALLD", "FALLU", "FALLB", "PRU", "PRD", "PRB",
          "GR", "SK", "KR")

abm_or_var <- abm_original[keep]


