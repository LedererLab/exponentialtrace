# Exponential trace model for discrete and continuous data

This repository provides an implementation of the sampling-based approximation for computing the maximum likelihood estimator (MLE) of exponential trace model. 

## Usage

The file `Implementation/exp_trace_model.R` contains functions for computing the approximate maximum likelihood estimator. Our upcoming paper describes the algorithm in details.

We include an example code `Implementation/example.R` for analyzing neuron spike data. The neuron spike data is from [Demas et al. 2003](http://www.jneurosci.org/content/23/7/2851.short).

> source("example.R")

## Simulation

Running the simulation, as described in the paper, takes a long time and is recommended to be implemented on a cluster. We include small-scale example code in the repository.

Data-type specific functions can be found under `Implementation/Data_type_specific.` 

### Poisson analog data
> source("Implementation/synthetic_Poisson_analog.R")

### Exponential analog data
> source("Implementation/synthetic_exponential_analog.R")

### Composite of Poisson analog and Bernoulli data
> source("Implementation/synthetic_Poisson_Bernoulli.R")

### Composite of Poisson analog and Gaussian data
> source("Implementation/synthetic_Poisson_Gaussian")

## Repository Authors

* **[Rui Zhuang](https://github.com/ruizhuanguw)** &mdash; Ph.D. candidate in Biostatistics, University of Washington &mdash; *methodology and `R` implementation*
* **[Noah Simon](http://faculty.washington.edu/nrsimon/)** &mdash; Assistant Professor in Biostatistics, University of Washington &mdash; *methodology*
* **[Johannes Lederer](https://johanneslederer.com/)** &mdash; Professor in Mathematical Statistics, Ruhr-University Bochum &mdash; *methodology*





