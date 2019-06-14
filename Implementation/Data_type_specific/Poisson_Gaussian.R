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
    stop("the simulated data is not a vector!")
  }

  T_matrix_value <- c(sqrt(input[1:p1]), input[(p1 + 1):p]) %*%
    t(c(sqrt(input[1:p1]), input[(p1 + 1):p]))

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

  diag(c(-log(diag(Tbar)[1:p1]), 1 / diag(Tbar)[(p1 + 1):p]))
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

  # Replace the Gaussian block of lambda to a pd matrix.
  edge_sign <- c(0, -1, 1)
  edge_vec <- sample(edge_sign, p2 ** 2, replace = TRUE,
                     prob = c(1 - 2 / p2, 1 / p2, 1 / p2))
  A <- matrix(edge_vec, nrow = p2, ncol = p2)
  diag(A) <- 0

  # Make the matrix symmetric
  A <- make_symmetry(A)

  max_degree <- max(apply(A, MARGIN = 1, FUN = function(x){
    sum(abs(x))
  }))
  A <- A * 1 / (max_degree + 0.1)
  diag(A) <- 1

  if (!corpcor::is.positive.definite(A)){
    stop("the parameter matrix is not positive definite!")
  }

  lambda[(p1 + 1):p, (p1 + 1):p] <- A

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

  # First p1 variables are Poisson and the remaining are Gaussian.
  MCsample <- matrix(nrow = p, ncol = 10000)
  for (i in 1:p1){
    MCsample[i, ] <- rpois(n = 10000, lambda = exp(-lambda0[i, i]))
  }
  for (j in (p1 + 1):p){
    MCsample[j, ] <- rnorm(n = 10000, mean = 0, sd = 1 / sqrt(lambda0[j, j]))
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
    # p1 is Poisson and the remaining p2 is Gaussian
    x[1:p1] <- floor(x[1:p1])
    Tx <- c(sqrt(x[1:p1]), x[(p1 + 1):p]) %o% c(sqrt(x[1:p1]), x[(p1 + 1):p])
    exponent <- -sum(parameter_matrix * Tx) - sum(lfactorial(x[1:p1]))
    return(exponent)
  }

  myControl <- MfU.Control(p, slice.lower = c(rep(0, p1), rep(-Inf, p2)))

  # Initial value
  initial_sample <- c(sample(0:10, p1, replace = TRUE),
                      sample(1:10, p2, replace = TRUE))

  my_samples <- MfU.Sample.Run(x = initial_sample, f = logTarget,
                               control = myControl,
                              uni.sampler = "slice", nsmp = sample_size)

  # Convert the generated sample to MCMC objective
  my_samples <- mcmc(my_samples)

  my_samples[, 1:p1] <- floor(my_samples[, 1:p1])

  sample_time <- proc.time() - sample_start_time

  print(paste("Sampling is done.", "Use", round(sample_time[3], 2), "s"))

  return(my_samples)
}
