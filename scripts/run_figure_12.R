required_packages <- c("ggplot2", "gridExtra", "grid", "compositions", 
                      "robCompositions", "zCompositions", "targets", "dplyr")

new_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]
if(length(new_packages)) install.packages(new_packages)

lapply(required_packages, library, character.only = TRUE)

source("scripts/functions/fit_comp_pca.R")
source("scripts/functions/help_functions.R")
source("scripts/functions/cond_scores_function.R")
source("scripts/functions/grad_function.R")
source("scripts/functions/simulation_functions.R")
source("scripts/functions/read_data_KL15_XRF.R")

set_1 <- get_color_palette()

# #*******reproduction plot12***********#
dir <- paste0(getwd(), "/data/Africa_NE_200/data/")
data_kl15_xrf <- read.table(paste0(dir,'data_KL15_XRF.txt'), header = TRUE, sep = "\t")
data_kl15_xrf <- rename(data_kl15_xrf, depth = User_ID)
data_kl15_agem <- read.table(paste0(dir,'data_KL15-2_smooth_spline_ages.txt'), header = TRUE, sep = "\t")
data_kl15_agem <- rename(data_kl15_agem, age = best)
results <- read_data_kl15_xrf(data_kl15_xrf, data_kl15_agem)
data <- results$data_kl15

n_simulations <- 100
n_observations <- 100
scales <- c(0.01)

simulation_results <- lapply(scales, function(scale) {
  run_complex_simulation(
    n_simulations = n_simulations,
    n_observations = n_observations,
    data = data * scale,
    complex_setting
  )
})

simulation_results_auto <- lapply(scales, function(scale) {
  run_complex_simulation(
    n_simulations = n_simulations,
    n_observations = n_observations,
    data = data * scale,
    complex_setting
  )
})

iteration_data <- data.frame(
  iteration = c(
    unlist(lapply(simulation_results[[1]]$pca_results, function(x) x$iteration)),
    unlist(lapply(simulation_results_auto[[1]]$pca_results, function(x) x$iteration))
  ),
  setting = factor(
    c(
      rep("uncorrelated scores", length(unlist(lapply(simulation_results[[1]]$pca_results, function(x) x$iteration)))),
      rep("scores with temporal trend", length(unlist(lapply(simulation_results_auto[[1]]$pca_results, function(x) x$iteration))))
    ),
    levels = c("uncorrelated scores", "scores with temporal trend")
  )
)

plot1 <- ggplot(iteration_data, aes(x = iteration, fill = setting)) +
  geom_histogram(aes(y = after_stat(density)), 
                 position = "dodge", 
                 bins = 30, alpha = 0.7) +
  geom_density(aes(color = setting), fill = NA, linewidth = 0.5, adjust = 1.3) +
  scale_fill_manual(values = c("#e0ecf4", "#9ebcda")) +
  scale_color_manual(values = c("#e0ecf4", "#9ebcda")) +
  theme_minimal() +
  labs(
    title = "Number of iterations until convergence",
    x = "Iterations",
    y = "Distribution"
  )

png("./scripts/figures/figure_12_upd.png", width = 12, height = 5, units = "in", res = 300)
print(plot1)
dev.off()