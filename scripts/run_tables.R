required_packages <- c("ggplot2", "dplyr")

new_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]
if(length(new_packages)) install.packages(new_packages)

lapply(required_packages, library, character.only = TRUE)

source("scripts/read_data_KL15_XRF.R")

dir <- paste0(getwd(), "/data/Africa_NE_200/data/")
data_kl15_xrf <- read.table(paste0(dir,'data_KL15_XRF.txt'), header = TRUE, sep = "\t")
data_kl15_xrf <- rename(data_kl15_xrf, depth = User_ID)
data_kl15_agem <- read.table(paste0(dir,'data_KL15-2_smooth_spline_ages.txt'), header = TRUE, sep = "\t")
data_kl15_agem <- rename(data_kl15_agem, age = best)
results <- read_data_kl15_xrf(data_kl15_xrf, data_kl15_agem)
data_kl15 <- results$data_kl15

data_sel <- data_kl15[4:ncol(data_kl15)-1]
colnames(data_sel) <- gsub("_cts", "", colnames(data_sel))

data_acomp <- acomp(data_sel)

var_matrix <- variation(data_acomp)
#**************table 7*******************#
(norm_var_matrix <- 1/sqrt(2) * var_matrix)

mvar(data_acomp)

D <- 13

xi_matrix <- norm_var_matrix / (ncol(data_acomp) * mvar(data_acomp))
#**************table 3*******************#
xi_norm <- xi_matrix * 13 * 12