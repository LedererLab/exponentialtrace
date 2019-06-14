
exp_trace <- function(input_data, parameter_matrix, reference_matrix){
  # Compute a vector of exp(-(M-M0, T(x))). Each element corresponding to
  #    one sample of p variables.
  #
  # Args:
  #   input_data: large size samples from the proposal density.
  #   parameter_matrix: A p*p parameter matrix of the target density.
  #   reference_matrix: The p*p parameter matrix of the proposal density.
  #
  # Returns:
  #   Vector corresponding to the exponential trace for each sample.

  delta_matrix <- parameter_matrix - reference_matrix

  exp_trace_value <- apply(input_data, 2,
                           function(x){
                             exp(-sum(delta_matrix * T_matrix(x)))
                           })

  return(exp_trace_value)
}

gradient <- function(MCsample, exp_trace_value){
  # Compute gradient for current step given the samples from
  #   proposal distribution.
  #
  # Args:
  #   MCsample: samples from the proposal distribution.
  #   exp_trace_value: a vector of the exponential trace for each sample.
  #
  # Returns:
  #   A p*p matrix corresponding to the gradient at current step.

  vector_numerator <- apply(MCsample, 2, T_matrix) * rep(1, p ** 2) %o%
    exp_trace_value

  numerator <- matrix(apply(vector_numerator, 1, mean), nrow = p, ncol = p)
  denominator <- mean(exp_trace_value)

  gradient_value <- Tbar - numerator / denominator

  return(gradient_value)
}

objective <- function(parameter_matrix, sufficient_data, exp_trace_value){
  # Compute the objective function for current step.
  #
  # Args:
  #   parameter_matrix: A p*p parameter matrix.
  #   sufficient_data: A p*p sufficient statistics matrix T(x).
  #   exp_trace_value: A vector of length p corresponding to
  #     the exponential trace value for each MC sample.
  #
  # Returns:
  #   Value of the objective function at the current step.

  objective_value <- sum(parameter_matrix * sufficient_data) +
    log(mean(exp_trace_value))

  return(objective_value)
}

gradient_descent <- function(sample, MCsample, constraint = FALSE,
                             t = 1, alpha = 0.3, beta = 0.5,
                             converge_threshold = 1e-4){
  # Implement gradient descent for computing the
  #   approximate maximum likelihood estimator.
  #
  # Args:
  #   sample: A n*p dataset.
  #   MCsample: Samples drawn from the proposal distribution with
  #     global parameter M0.
  #   constraint: TRUE/FALSE or start of continuous variables determining
  #     applied constraint.
  #   t: Initial step size
  #   alpha: Fraction of decreasing for backtracking line search.
  #     recommended range: (0, 0.5).
  #   beta: Fraction of step size shrinking for backtracking line search.
  #     recommended range: (0, 1).
  #   converge_threshold: Criterion for convergence.
  #
  # Returns:
  #   Maxmimum likelihood estimator of the parameter matrix M.

  # Initiation
  i <- 0
  converge <- FALSE

  parameter_update <- lambda0

  exp_trace_value_update <- exp_trace(input_data = MCsample,
                                      parameter_matrix = parameter_update,
                                      reference_matrix = lambda0)

  objective_value_update <- objective(parameter_matrix = parameter_update,
                                      sufficient_data = Tbar,
                                      exp_trace_value = exp_trace_value_update)

  objective_values <- objective_value_update

  # Start estimation
  while (converge == FALSE & i < maxit){

    i <- i + 1

    parameter_current <- parameter_update

    exp_trace_value_current <- exp_trace_value_update

    gradient_value <- gradient(MCsample = MCsample,
                               exp_trace_value = exp_trace_value_current)

    objective_value_current <- objective_value_update

    # Update
    parameter_update <- parameter_current - t * gradient_value
    
    if (constraint == TRUE){
      while (det(parameter_update) < 0){
        t <- beta * t
        parameter_update <- parameter_current - t * gradient_value
      }
    }else if (is.numeric(constraint)){
      while (det(parameter_update[constraint:p, constraint:p]) < 0){
        t <- beta * t
        parameter_update <- parameter_current - t * gradient_value
      }
    }
    
    exp_trace_value_update <- exp_trace(input_data = MCsample,
                                        parameter_matrix = parameter_update,
                                        reference_matrix = lambda0)

    objective_value_update <- objective(parameter_matrix = parameter_update,
                                        sufficient_data = Tbar,
                                        exp_trace_value =
                                          exp_trace_value_update)

    while (objective_value_update >
           objective_value_current - alpha * t * sum(gradient_value ** 2)){
      t <- beta * t

      parameter_update <- parameter_current - t * gradient_value

      exp_trace_value_update <- exp_trace(MCsample,
                                          parameter_matrix = parameter_update,
                                          reference_matrix = lambda0)

      objective_value_update <- objective(parameter_matrix = parameter_update,
                                          sufficient_data = Tbar,
                                          exp_trace_value =
                                            exp_trace_value_update)
    }

    objective_values <- c(objective_values, objective_value_update)

    # Check convergence
    change_value <- objective_value_update - objective_value_current
    converge <- abs(change_value) <= converge_threshold

    # Print progress
    if (i %% 200 == 0){
      print(paste("iteration", i))
      print(paste("step size", round(t, 5)))
      print(paste("change in negative loglikelihood",
                  round(change_value, digits = 5)))
    }

  }

  return(list(estimator = parameter_update))
}

transform_data <- function(index_neuron, time_interval){
  # Prepare the neuron spike data.
  #
  # Args:
  #   index_neuron: The neuron to be processed
  #   time_interval: The length of time bin to compute count of spikes
  #
  # Returns:
  #   A n*2 dataframe with row representing the time bin and
  #     column representing the count of spike in the time bin.
  #
  # Raises:
  #   Warning if the spike counts before and after processing do not match.

  dat <- wide_spikes[index_neuron]

  # time from second to ms.
  dat <- unlist(dat) * 1000

  # allocate spikes to time intervals
  dat_df <- table(cut(dat, breaks = seq(from = 0,
                                        to = ceiling(max(dat) / time_interval) *
                                          time_interval,
                                        by = time_interval)))

  dat_df <- data.frame(dat_df) %>% mutate(Var1 = as.character(Var1))
  colnames(dat_df) <- c("time", neuron_names[index_neuron])

  # double check
  if (sum(dat_df[, -1]) != length(dat)){
    warning("Spike counts do not match!")
  }

  return(dat_df)

}

my_qgraph <- function(est, num_edge, title){
  # Plot the undirected graph obtained by thresholding
  #   the maximum likelihood estimator.
  #
  # Args:
  #   est: An estimator of parameter matrix M.
  #   num_edge: Pre-specified number of edges.
  #   title: Title of the graph.
  #
  # Returns:
  #   A graph recovered by thresholding the estimator.

  off_entries <- est[upper.tri(est, diag = FALSE)]
  off_ordered <- abs(off_entries)[order(abs(off_entries))]

  pos_cut <- length(off_ordered) - num_edge
  cut_point <- mean(off_ordered[pos_cut : (pos_cut + 1)])

  adj_matrix <- (abs(est) > cut_point) * sign(est)

  qgraph::qgraph(adj_matrix,
                 edge.labels = FALSE,
                 diag = FALSE,
                 labels = neuron_names,
                 layout = epos,
                 title = title,
                 vsize = 10,
                 esize = 6,
                 label.prop = 0.92,
                 normalize = FALSE)
}

connected_dist <- function(est, num_edge){
  # Compute the mean Eucleadian distance of directly connected neurons
  #   recovered by thresholding the estimator.
  #
  # Args:
  #   est: An estimator of parameter matrix M.
  #   num_edge: Pre-specified number of edges.
  #
  # Returns:
  #   Mean Eucleadian distance using (x,y) epos information for
  #     directly connected neurons.

  off_entries <- est[upper.tri(est, diag = FALSE)]
  off_ordered <- abs(off_entries)[order(abs(off_entries))]

  pos_cut <- length(off_ordered) - num_edge
  cut_point <- mean(off_ordered[pos_cut : (pos_cut + 1)])

  # connected ones
  adj_matrix <- (abs(est) > cut_point) * sign(est)
  diag(adj_matrix) <- 0

  connect_tbl <- data.frame(matrix(ncol = 6))
  colnames(connect_tbl) <- c("neuron1_name", "neuron2_name",
                             "neuron1_loc", "neuron2_loc",
                             "estimator", "Euclidean dist")

  connect_pos <- data.frame(which( adj_matrix != 0 & upper.tri(adj_matrix),
                                   arr.ind = T))

  for (i in 1: nrow(connect_pos)){
    idx1 <- connect_pos[i, "row"]
    idx2 <- connect_pos[i, "col"]

    connect_tbl <- connect_tbl %>% add_row(neuron1_name = neuron_names[idx1],
                                           neuron2_name = neuron_names[idx2],
                                           neuron1_loc = paste(epos[idx1, ],
                                                               collapse = ","),
                                           neuron2_loc = paste(epos[idx2, ],
                                                               collapse = ","),
                                           estimator = est[idx1, idx2],
                                           `Euclidean dist` =
                                             dist(rbind(epos[idx1, ],
                                                        epos[idx2, ])))
  }

  return(connect_tbl[-1, ])
}

get_lower_tri <- function(input_matrix){
  # Get the lower triangle of the input matrix.
  #
  # Args:
  #   input_matrix: A matrix.
  #
  # Returns:
  #   Vector of the lower triangle.

  input_matrix[lower.tri(input_matrix)]
}

make_symmetry <- function(input_matrix){
  # Mapping the upper triangle of input matrix to the lower.
  #
  # Args:
  #   input_matrix: A matrix.
  #
  # Returns:
  #   A symmetric matrix.

  input_matrix[lower.tri(input_matrix)] <-
    t(input_matrix)[lower.tri(input_matrix)]
  return(input_matrix)
}
