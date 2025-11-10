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

data_split <- data |>
  group_by(.data[[center_name]]) |>
  group_split()

n_centers <- length(data_split)
n_data <- sum(unlist(lapply(data_split, nrow)))
p <- length(covariates)


### Test local BCA
center_identity <- 1
bstats_1 <- bca_iterate_sites(outcome, covariates, model, family, data_split[1],
                              use_local_intercepts = FALSE, center_name = NULL)
fit.bca.local <- bayes_local_glm(bstats_1, center_identity=center_identity, family=family, alpha=0.05)
fit.glm.local <- glm(data=data_split[[center_identity]], formula = model)
a <- as.vector(fit.bca.local$Estimate)
b <- as.vector(fit.glm.local$coefficients)

expect_equal(a, b,
             info = "bayes_local_glm() should return same as glm()")