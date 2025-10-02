library(dplyr)
library(tidyr)
library(invgamma)

# Pre-processing wrapper for BaySeq methods
bayseq_prepare <- function(target, covariates, Method, data_split, n_centers, use_local_intercepts, center_name=NULL) {

    bstats <- vector("list", n_centers)
    n_centers <- n_centers

    for (i in seq_along(data_split)) {

        # For GLM regression, we locally compute GLM parameters

        if (family != "gaussian") {

            if (use_local_intercepts) {
                stop("Not yet implemented")
            } else {
                res <- glm(Method, family, data_split[[i]])
                bstats[[i]] <- list("beta" = res$coefficients, "sigma" = vcov(res))
            }

        # For Bayesian linear regression, we locally compute sufficient statistics

        } else if (family == "gaussian") {
            mat <- data_split[[i]]

            if (use_local_intercepts) {
                stopifnot(!is.null(center_name))
                covariates_local <- covariates[covariates != center_name]
                bstats[[i]] <- bayes_lin_reg_stats(mat, target, covariates_local, k=i, n_centers=n_centers)
            } else {
                bstats[[i]] <- bayes_lin_reg_stats(mat, target, covariates)
            }

        } else {
            stop(paste0("Invalid family:", family))
        }
    }
    bstats
}

# Compute summary statistics for Bayesian linear regression
bayes_lin_reg_stats <- function(mat, target, covariates, weights=NULL, k=0, n_centers=0) {

    if (k == 0) {
        # Homogeneous setting: 1 global intercept
        mat["Intercept"] <- 1
        intercept_colnames <- "Intercept"
    } else if (k > 0 && n_centers > 0) {
        # Heterogeneous setting: 1 intercept per local center
        intercepts <- matrix(0, nrow(mat), n_centers)
        intercepts[, k] <- 1
        mat_new <- cbind(mat, intercepts)
        intercept_colnames <- paste0("Intercept_", seq_len(ncol(intercepts)))
        colnames(mat_new) <- c(colnames(mat), intercept_colnames)
        mat <- mat_new
    } else {
        stop("k must be >= 0 and n_centers must be > 0 if k > 0")
    }

    x <- mat[, c(intercept_colnames, covariates)]

    if (is.null(weights)) {
        weights <- diag(1, dim(x))
    }

    y <- mat[[target]]
    yy <- t(y) %*% weights %*% y

    x <- t(t(x))
    xx <- t(x) %*% weights %*% x

    xy <- t(x) %*% weights %*% y

    n <- nrow(mat)  # sample size

    list("xx" = xx,
        "xy" = xy,
        "yy" = yy,
        "n" = n)
}

# Compute posterior parameters for Bayesian linear regression
bayes_lin_reg_post_params <- function(bayes_stats_list, prior_params) {

    mu0 <- prior_params$mu0
    lambda0 <- prior_params$lambda0
    a0 <- prior_params$a0
    b0 <- prior_params$b0

    sum_xx <- Reduce("+", lapply(bayes_stats_list, function(x) x$xx))
    sum_xy <- Reduce("+", lapply(bayes_stats_list, function(x) x$xy))
    sum_yy <- Reduce("+", lapply(bayes_stats_list, function(x) x$yy))
    sum_n <- Reduce("+", lapply(bayes_stats_list, function(x) x$n))

    lambda_l <- lambda0 + sum_xx
    mu_l <- solve(lambda_l) %*% (lambda0 %*% mu0 + sum_xy)
    a_l <- a0 + sum_n / 2
    b_l <- b0 + 0.5 * sum_yy + 0.5 * (t(mu0) %*% lambda0 %*% mu0 - t(mu_l) %*% lambda_l %*% mu_l)

    list("lambda_l" = lambda_l,
        "mu_l" = mu_l,
        "a_l" = a_l,
        "b_l" = b_l
    )
}

# Compute maximum a posteriori estimates for Bayesian linear regression
bayes_lin_reg_post_map <- function(bayes_post_params, p) {
    beta_l <- bayes_post_params$mu_l
    a_l <- bayes_post_params$a_l
    b_l <- bayes_post_params$b_l
    lambda_l <- bayes_post_params$lambda_l
    tau_l <- as.numeric(a_l / b_l)
    sigma_l <- solve(tau_l * lambda_l)
    list("beta_l" = beta_l, "sigma_l" = sigma_l, "lambda_l"=lambda_l)
}

# Define prior for Bayesian linear regression
get_linreg_prior <- function(covariates, use_local_intercepts, n_centers, epsilon=1e-10) {

    if (use_local_intercepts) {
        p <- length(covariates) + n_centers
    } else {
        p <- length(covariates)
    }

    if (use_local_intercepts) {
        prior_params <- list("mu0" = as.matrix(rep(0, length(covariates) + n_centers)),
                "lambda0" = epsilon * diag(length(covariates) + n_centers),
                "a0" = epsilon -p/2,
                "b0" = epsilon
                )
    } else {
        prior_params <- list("mu0" = as.matrix(rep(0, length(covariates) + 1)),
                "lambda0" = epsilon * diag(length(covariates) + 1),
                "a0" = epsilon -p/2,
                "b0" = epsilon
                )
    }
    prior_params
}

# Compute BaySeq parameters using one-shot approach
bayseq_oneshot <- function(bstats, n_centers, use_local_intercepts, family,
                            covariates, epsilon=1e-10, center_name=NULL, CI="t") {
    params_seq <- list()

    if (family != "gaussian") {
        print("Normal Method known variance")

        update_normal_known_variance <- function(beta, sigma) {
            sigma_post <- solve(Reduce(`+`, lapply(sigma, solve)))
            beta_post <- sigma_post %*% Reduce(`+`, Map(function(s, b) solve(s) %*% b, sigma, beta))
            list("sigma" = sigma_post, "beta" = beta_post)
        }
        beta <- lapply(bstats, function(x) unlist(unname(x["beta"])))
        sigma <- lapply(bstats, function(x) x[["sigma"]])
        updated_params <- update_normal_known_variance(beta, sigma)
        params_seq$sigma <- updated_params$sigma
        params_seq$beta <- updated_params$beta

    } else {
        print("Bayesian linear regression")
        stopifnot(!is.null(center_name))
        covariates_local <- covariates[covariates != center_name]
        prior_params <- get_linreg_prior(covariates_local, use_local_intercepts, n_centers, epsilon=epsilon)

        if (use_local_intercepts) {
            p <- length(covariates) + n_centers
        } else {
            p <- length(covariates)
        }

        bayes_post_params <- bayes_lin_reg_post_params(bstats, prior_params)
        bayes_map <- bayes_lin_reg_post_map(bayes_post_params, p)
        params_seq$beta <- bayes_map$beta_l
        params_seq$sigma <- bayes_map$sigma_l
        params_seq$a_l <- bayes_post_params$a_l
        params_seq$b_l <- bayes_post_params$b_l
        params_seq$lambda_l <- bayes_map$lambda_l

        # use inverse of mean to get equivalent sigma2 to lm
        params_seq$dispersion <- as.numeric(bayes_post_params$b_l / (bayes_post_params$a_l))
        disp_ci <- qinvgamma(c(0.025, 0.975), shape = bayes_post_params$a_l, rate = bayes_post_params$b_l)
    }

    if (CI == "t") {
        params_seq$CI <- get_bayes_linreg_ci_t(params_seq)
    } else if (CI == "normal") {
        params_seq$CI <- get_bayes_linreg_ci_normal(params_seq)
    }

    if (family == "gaussian") {
        new_row <- data.frame(disp_ci[[1]], disp_ci[[2]], row.names = "sigma2")
        colnames(new_row) <- colnames(params_seq$CI)
        params_seq$CI <- rbind(params_seq$CI, new_row)
    }
    params_seq
}

# Compute credible intervals for Bayesian linear regression
get_bayes_linreg_ci_normal <- function(params_seq, alpha=0.05) {

    # normal
    z <- qnorm(1 - alpha / 2)  # ≈ 1.96 for 95% CI
    lower <- params_seq$beta - z * sqrt(diag(params_seq$sigma))
    upper <- params_seq$beta + z * sqrt(diag(params_seq$sigma))
    ci <- data.frame(
        lower = lower,
        upper = upper
    )
    ci
}

# Use marginal t-distributions for CI
get_bayes_linreg_ci_t <- function(params_seq, alpha=0.05) {
    mu <- params_seq$beta
    a <- as.vector(params_seq$a_l)
    b <- as.vector(params_seq$b_l)
    lambda_diag <- diag(params_seq$lambda_l)

    # Degrees of freedom
    nu <- 2 * a  # == n - m

    scale <- (lambda_diag * a) / b
    variance <- (1/scale) * nu / (nu-2) # B-S p. 435
    stdev <- sqrt(variance)

    t_crit <- qt(1-alpha/2, df = nu)
    lower <- mu - t_crit * stdev
    upper <- mu + t_crit * stdev

    data.frame(
        lower = lower,
        upper = upper
    )
}


tidy_results <- function(param_seq, use_local_intercepts) {
    params_seq_all <- rbind(params_seq$beta, sigma2 = params_seq$dispersion)
    if (use_local_intercepts) {
        params_seq_all <- params_seq_all[!grepl("Intercept_\\d+$", rownames(params_seq_all)), , drop=FALSE]
    }
    df <- data.frame(t(params_seq_all), check.names = FALSE)  # prevent renaming (Intercept) to X.Intercept.
    df$Method <- "BCA"
    df <- df %>% pivot_longer(-Method, names_to = "Covariate", values_to = "Estimate")

    params_seq$CI$Method <- "BCA"
    if (!("Covariate" %in% colnames(params_seq$CI)))
        params_seq$CI <- tibble::rownames_to_column(params_seq$CI, var = "Covariate")

    df_merged <- left_join(df, params_seq$CI, by = c("Method", "Covariate"))
    df_merged
}