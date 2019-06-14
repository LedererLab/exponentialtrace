library(rhdf5)
library(dplyr)

dataset_name <- "Demas2003_6WK_CTRL_ME2.h5"

# Config ------------------------------------------------------------------
# Time interval of interest (ms)
interval <- 40
 
# Estimation config
sampling_seed <- 1
minit <- 1
maxit <- 2000
converge_threshold <- 1E-6

# Evaluation config
num_edge <- 30

dir_path <- dirname(rstudioapi::getSourceEditorContext()$path)

filename <- file.path(dir_path, "../ExampleData", dataset_name)

# Input functions --------------------------------------------------------------
source(file.path(dir_path, "exp_trace_model.R"))
source(file.path(dir_path, "Data_type_specific", "poisson_analog.R"))

# Data prepare ------------------------------------------------------------
epos <- h5read(filename, name = "epos")
scount <- h5read(filename, name = "sCount")
spikes <- h5read(filename, name = "spikes")

# Get spikes dataset
neuron_names <- h5read(filename, name = "names")
wide_spikes <- split(spikes, f = rep(neuron_names, times = scount))

# Merge datasets
my_data <- transform_data(index_neuron = 1, time_interval = interval)

for (i in 2:nrow(epos)){
  my_data <- my_data %>%
    full_join(transform_data(index_neuron = i, time_interval = interval),
              by = "time") %>%
    replace(is.na(.), 0)
}

# Double check
if (!all(neuron_names == colnames(my_data)[-1]) |
    !all(apply(my_data[, -1], MARGIN = 2, FUN = sum, na.rm = T) == scount) |
    !ceiling(max(spikes) * 1000 / interval) == nrow(my_data)) {
  warning("Something is wrong!")
}

my_samples <- my_data[, -1]

print("distribution of num of spikes of each neurons")
print(apply(my_samples, MARGIN = 2, max))

# Analyze -----------------------------------------------------------------
if (file.exists(paste(filename, "estimator.RData", sep = "_"))){
  load(file = paste(filename, "estimator.RData", sep = "_"))
}else{
  p <- ncol(my_samples)

  # Compute sufficient statistics for solving MLE
  Tbar <- apply(apply(my_samples, 1, T_matrix), 1, FUN = mean)
  Tbar <- matrix(Tbar, p, p)

  # Solve the MLE
  lambda0 <- get_lambda0(Tbar)

  # Sampling for approximating the normalization term
  set.seed(sampling_seed)

  MCsample <- get_mcsample(lambda0)

  if (nrow(MCsample) != p){
    stop("MC sample has dimension problem!")
  }

  print("Solving MLE...")

  estimate_start_time <- proc.time()

  rlts <- gradient_descent(sample = my_samples, MCsample = MCsample,
                           converge_threshold = converge_threshold)

  estimate_time <- proc.time() - estimate_start_time

  print(paste("Estimation is done.", "Use", round(estimate_time[3], 2), "s"))

  # Save the estimation results
  save(rlts, file = paste(filename, "estimator.RData", sep = "_"))
}

# Evaluate the results
pdf(paste(filename, "plot.pdf", sep = "_"), width = 19, height = 8)
par(mfrow = c(1, 2))

# Results of ETM
dist_table_etm <- connected_dist(rlts$estimator, num_edge = num_edge)
mean_dist_etm <- mean(dist_table_etm$`Euclidean dist`)
my_qgraph(rlts$estimator, num_edge, title = paste("ETM", round(mean_dist_etm)))

# Results of GGM
est_ggm <- solve(cov(my_samples))
dist_table_ggm <- connected_dist(est_ggm, num_edge = num_edge)
mean_dist_ggm <- mean(dist_table_ggm$`Euclidean dist`)
my_qgraph(est_ggm, num_edge, title = paste("GGM", round(mean_dist_ggm)))

dev.off()
