T_matrix <- function(input){
  # Compute the sufficient statistics for square root transformation.
  #
  # Args:
  #   input: A vector of p variables.
  #
  # Returns:
  #   A p*p sufficient statistics matrix.
  #
  # Raises:
  #   Error if the input data is not a vector.
  
  if (!is.null(dim(input))){
    stop("the input data is not a vector!")
  }
  
  T_matrix_value <- sqrt(input) %*% t(sqrt(input))
  
  return(T_matrix_value)
}

get_lambda0 <- function(Tbar){
  # Compute the proposal diagonal matrix.
  #   with parameter matrix M0.
  #
  # Args:
  #   lambda0: Parameter matrix for the proposal distribution.
  #
  # Returns:
  #   A p*p proposal matrix.

  diag(1 / diag(Tbar))
}

get_lambda <- function(p, graph_type, negative_prob){
  # Generate random parameter matrix
  #
  # Args:
  #   p: Number of nodes.
  #   graph_type: Type of graphs, e.g., "random graph".
  #   negative_prob: Probability of negative interaction.
  #
  # Returns:
  #   A p*p parameter matrix.

  edge_sign <- c(0, -1, 1)

  # Generate adjacency matrix
  edge_vec <- sample(edge_sign, p ** 2, replace = TRUE,
                     prob = c(1 - 2 / p, 1 / p, 1 / p))
  A <- matrix(edge_vec, nrow = p, ncol = p)
  diag(A) <- 0

  # Make the matrix symmetric
  A <- make_symmetry(A)

  s <- max(apply(A, MARGIN = 1, FUN = function(x){
    sum(abs(x))
  }))
  A <- A * 1 / (s + 0.1)
  diag(A) <- 1

  if (!corpcor::is.positive.definite(A)){
    stop("the parameter matrix is not positive definite!")
  }

  return(A)
}

get_mcsample <- function(lambda0){
  # Generate samples from proposal distribution
  #   with parameter matrix M0.
  #
  # Args:
  #   lambda0: Parameter matrix for the proposal distribution.
  #
  # Returns:
  #   A p*n data matrix.

  MCsample <- matrix(nrow = p, ncol = 10000)
  for (i in 1:p){
    MCsample[i, ] <- rexp(n = 10000, rate = lambda0[i, i])
  }

  return(MCsample)
}

require(coda)
suppressMessages(require(MfUSampler))

my_slice_sampler <- function(parameter_matrix, sample_size) {
  # Slice sampler
  #
  # Args:
  #   parameter matrix: A p*p matrix M.
  #   sample_size: Number of samples
  #
  # Returns:
  #   A MCMC object.

  print("Generating data using slice sampler...")

  sample_start_time <- proc.time()

  # Log likelihood function
  logTarget <- function(x) {
    Tx <- sqrt(x) %o% sqrt(x)
    exponent <- -sum(parameter_matrix * Tx)
    return(exponent)
  }

  myControl <- MfU.Control(p, slice.lower = c(rep(0, p)))

  # Initial value
  initial_sample <- sample(1:100, p)

  my_samples <- MfU.Sample.Run(x = initial_sample, f = logTarget,
                               control = myControl,
                              uni.sampler = "slice", nsmp = sample_size)

  # Convert the generated sample to mcmc objective
  my_samples <- mcmc(my_samples)

  sample_time <- proc.time() - sample_start_time

  print(paste("Sampling is done.", "Use", round(sample_time[3], 2), "s"))

  return(my_samples)
}
