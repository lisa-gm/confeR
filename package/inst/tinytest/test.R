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

bstats <- bca_iterate_sites(outcome, covariates, model, family, data_split, use_local_intercepts, center_name)
params_oneshot <- bca_oneshot(bstats, n_centers, use_local_intercepts, use_local_variances, family,
                              epsilon = 1e-10, center_name)

### Test reduced params from Reverse-Bayes

params_minus_l_reverse <- get_reduced_params(center_identity, params_oneshot, bstats, family)
params_minus_l_forward <- bca_oneshot(bstats[-center_identity], n_centers-1,
                                use_local_intercepts, use_local_variances, family,
                                epsilon = 1e-10, center_name)

tinytest::expect_equal(params_minus_l_reverse$lambda_minus_l, params_minus_l_forward$lambda_l,
            info = "Reverse-Bayes should return identical results as Forward-Bayes minus l")
tinytest::expect_equal(params_minus_l_reverse$beta_minus_l, params_minus_l_forward$beta_l,
            info = "Reverse-Bayes should return identical results as Forward-Bayes minus l")
tinytest::expect_equal(params_minus_l_reverse$a_minus_l, params_minus_l_forward$a_l,
            info = "Reverse-Bayes should return identical results as Forward-Bayes minus l")
tinytest::expect_equal(params_minus_l_reverse$b_minus_l, params_minus_l_forward$b_l,
            info = "Reverse-Bayes should return identical results as Forward-Bayes minus l")