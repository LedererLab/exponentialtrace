
# Config ------------------------------------------------------------------
# Data type
datatype <- "Poisson"
# Number of node.
p <- 20
# Diagonal entry of parameter matrix M.
diagentry <- -1
# Off-diagonal entry of parameter matrix M.
offentry <- 0.3
# Sample size of observations
samplesize <- 250
# Number of random graphs for averaging. 
num_rand_seed <- 2
# Graph type: Erdos Renyi random graph.
graph_type <- "random"
# probability of negative interaction.
negative_prob <- 1 / 2

# Sampler configuration
burnin <- 20000
thin <- 4
num_chain <- 1

# Estimation configuration
minit <- 1
maxit <- 2000
converge_threshold <- 1E-6

# Random seed
start_seed <- 1

# Directory for output
file_path <- dirname(rstudioapi::getSourceEditorContext()$path)

# Read in functions -------------------------------------------------------
require(rlist)
source(file.path(file_path, "exp_trace_model.R"))
source(file.path(file_path, "Data_type_specific", "poisson_analog.R"))

# Implement ---------------------------------------------------------------
setting_name <- paste(datatype, "_",
                      graph_type,
                      "_node", p,
                      "_diag", diagentry,
                      "_off", offentry,
                      "_samplesize", samplesize,
                      "_trials", sprintf("%03d", num_rand_seed),
                      "_startSeed", start_seed, sep = "")

results_path <- paste(file_path, "/",
                      setting_name, "_",
                      Sys.Date(), sep = "")

# Create the results_path if needed
if (!file.exists(results_path)){
  dir.create(results_path)
}

setwd(results_path)

graph_type_ETM <- graph_type_precision <- graph_type_label <- list()

set.seed(start_seed)

start_value <- sample(1:100, 1)

for (rand_seed in ( (start_value + 1):(start_value + num_rand_seed))){
  # Set random seed for parameter matrix generation.
  set.seed(rand_seed)

  lambda <- get_lambda(p = p, graph_type = graph_type,
                       negative_prob = negative_prob,
                       diagentry = diagentry,
                       offentry = offentry)

  # Set random seed for data generation.
  set.seed(rand_seed)

  # Slice sampling
  chain_length <- (burnin + (samplesize - 1) * (thin + 1) + 1) / num_chain

  my_chain <- my_slice_sampler(parameter_matrix = lambda,
                               sample_size = chain_length)

  # Remove the burnin part
  my_samples <- fitR::burnAndThin(my_chain, burn = burnin, thin = thin)

  if (!corpcor::is.positive.definite(cov(my_samples))){
    next
  }

  # Compute sufficient statistics for solving MLE
  Tbar <- apply(apply(my_samples, 1, T_matrix), 1, FUN = mean)
  Tbar <- matrix(Tbar, p, p)

  # Solve the MLE
  lambda0 <- get_lambda0(Tbar)

  # sampling for approximating the normalization term
  sampling_seed <- rand_seed

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

  MLE_hat <- rlts$estimator

  precision_hat <- solve(cov(my_samples))

  graph_type_ETM <- list.append(graph_type_ETM,
                                abs(get_lower_tri(MLE_hat)))

  graph_type_precision <- list.append(graph_type_precision,
                                      abs(get_lower_tri(precision_hat)))

  graph_type_label <- list.append(graph_type_label,
                                  abs(get_lower_tri(lambda)) > 1E-5)

}

# Calculate delta in AUC
pred_graph_type_ETM <- ROCR::prediction(graph_type_ETM, graph_type_label)
perf_graph_type_ETM <- ROCR::performance(pred_graph_type_ETM, measure = "auc")
auc_etm_list <- unlist(perf_graph_type_ETM @y.values)

pred_graph_type_precision <- ROCR::prediction(graph_type_precision,
                                              graph_type_label)
perf_graph_type_precision <- ROCR::performance(pred_graph_type_precision,
                                               measure = "auc")
auc_gm_list <- unlist(perf_graph_type_precision @y.values)

mean_delta_gm <- mean(auc_etm_list) - mean(auc_gm_list)

print(paste("The average difference in AUC between ETM and GM of",
            num_rand_seed, "random graphs is :", mean_delta_gm))

# Save the working image.
save.image(file = file.path(results_path, paste0(setting_name, ".RData")))

# Plot average ROC curve --------------------------------------------------------------
pdf(file = file.path(results_path, "plot.pdf"), width = 5, height = 4)
par(mfrow = c(1, 1))

perf_graph_type_ETM <- ROCR::performance(pred_graph_type_ETM, "tpr", "fpr")
ROCR::plot(perf_graph_type_ETM, avg = "horizontal", lwd = 2.5,
           main = paste0(graph_type, " graph : sample size = ", samplesize,
                         "\n AUC(etm) - AUC(precision) = ",
                         round(mean_delta_gm, 3)))

perf_graph_type_precision <- ROCR::performance(pred_graph_type_precision,
                                               "tpr", "fpr")
ROCR::plot(perf_graph_type_precision,
           avg = "horizontal", col = 2, add = T, lwd = 2.5)
legend("bottomright", c("ETM", "Precision"), col = 1:2, lty = 1, cex = 0.8)

dev.off()
