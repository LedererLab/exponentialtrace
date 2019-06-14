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

  diag(c(-log(diag(Tbar)[1:p1]), -log(1 / (1 - diag(Tbar)[(p1 + 1):p]) - 1)))
}

get_lambda <- function(p, graph_type, negative_prob, diagentry, offentry){
  # Generate random parameter matrix
  #
  # Args:
  #   p: Number of nodes.
  #   graph_type: Type of graphs, e.g., "random graph".
  #   negative_prob: Probability of negative interaction.
  #   diagentry: Value of diagonal entries.
  #   offentry: Value of off-diagonal entries.
  #
  # Returns:
  #   A p*p parameter matrix.

  # Generate adjacency matrix.
  lambda <- huge::huge.generator(n = 10, d = p, prob = 2 / p,
                                 graph = graph_type, vis = F)$theta

  random_sign <- sample(c(-1, 1), p ** 2, replace = TRUE,
                        prob = c(negative_prob, 1 - negative_prob))

  offentry_matrix <- matrix(random_sign, nrow = p) * offentry

  # Make the matrix symmetric
  offentry_matrix <- make_symmetry(offentry_matrix)

  lambda <- as.matrix(lambda) * offentry_matrix

  diag(lambda) <- diagentry

  return(lambda)
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

  for (i in 1:p1){
    MCsample[i, ] <- rpois(n = 10000, lambda = exp(-diag(lambda0)[i]))
  }
  for (j in (p1 + 1):p){
    MCsample[j, ] <- rbinom(n = 10000, size = 1,
                           prob = exp(-lambda0[j, j]) / (1 +
                                                           exp(-lambda0[j, j])))
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
    # Floor the first p1 is Poisson and the remaining p2 is Bernouli.
    x <- floor(x)
    Tx <- sqrt(x) %o% sqrt(x)
    exponent <- -sum(parameter_matrix * Tx) - sum(lfactorial(x[1:p1]))
    return(exponent)
  }

  myControl <- MfU.Control(p, slice.lower = c(rep(0, p)),
                           slice.upper = c(rep(Inf, p1), rep(2 - 1E-10, p2)))

  # Initial value
  initial_sample <- c(sample(1:100, p1, replace = TRUE),
                      sample(0:1, p2, replace = TRUE))

  my_samples <- MfU.Sample.Run(x = initial_sample, f = logTarget,
                               control = myControl, uni.sampler = "slice",
                               nsmp = sample_size)

  # Convert the generated sample to MCMC objective
  my_samples <- mcmc(my_samples)

  my_samples <- floor(my_samples)

  sample_time <- proc.time() - sample_start_time

  print(paste("Sampling is done.", "Use", round(sample_time[3], 2), "s"))

  return(my_samples)
}
