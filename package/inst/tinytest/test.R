library(tinytest)
library(dplyr)

## Load test data
data <- read.csv("testdata-sim-nurses-1-out.csv")
model <- stress ~ gender + age + experience + wardtype
center_name <- "hospital"
family <- "gaussian"
covariates <- attr(terms(model), "term.labels")
outcome <- all.vars(model)[1]
sites <- levels(data[[center_name]])
use_local_intercepts <- FALSE
use_local_variances <- FALSE
data_split <- data |>
  group_by(.data[[center_name]]) |>
  group_split()

n_centers <- length(data_split)
n_data <- sum(unlist(lapply(data_split, nrow)))
p <- length(covariates)

bstats <- bca_iterate_sites(
  outcome = outcome,
  covariates = covariates,
  model = model,
  family = family,
  data_split = data_split,
  use_local_intercepts = use_local_intercepts
)

params_oneshot <- bca_oneshot(
  sumstats = bstats,
  family = family,
  use_local_intercepts = use_local_intercepts,
  use_local_variances = use_local_variances,
  epsilon = 1e-10
)

### Test reduced params from Reverse-Bayes
center_identity <- 1
params_minus_l_reverse <- get_reduced_params(center_identity, params_oneshot, bstats, family)
params_minus_l_forward <- bca_oneshot(
  sumstats = bstats[-center_identity],
  use_local_intercepts = use_local_intercepts,
  use_local_variances = use_local_variances,
  family = family,
  epsilon = 1e-10
)

tinytest::expect_equal(params_minus_l_reverse$lambda_minus_l, params_minus_l_forward$lambda_l,
  info = "Reverse-Bayes should return identical results as Forward-Bayes minus l"
)
tinytest::expect_equal(params_minus_l_reverse$beta_minus_l, params_minus_l_forward$beta_l,
  info = "Reverse-Bayes should return identical results as Forward-Bayes minus l"
)
tinytest::expect_equal(params_minus_l_reverse$a_minus_l, params_minus_l_forward$a_l,
  info = "Reverse-Bayes should return identical results as Forward-Bayes minus l"
)
tinytest::expect_equal(params_minus_l_reverse$b_minus_l, params_minus_l_forward$b_l,
  info = "Reverse-Bayes should return identical results as Forward-Bayes minus l"
)
