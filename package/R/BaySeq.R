library(dplyr)
library(tidyr)
library(invgamma)

#' @title BCA iterate over local sites
#'
#' @description Pre-processing wrapper for BCA methods.
#'
#' @param outcome Outcome vector y
#'
#' @return List over all local sites with transmitted summary statistics
#'
#' @author Peter Degen
#'
#' @export
bca_iterate_sites <- function(outcome, covariates, Method, data_split,
                              use_local_intercepts, center_name = NULL) {

    n_sites <- length(data_split)
    bstats <- vector("list", n_sites)

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
                bstats[[i]] <- bayes_lin_reg_stats(mat, outcome, covariates_local, k = i, n_sites = n_sites)
            } else {
                bstats[[i]] <- bayes_lin_reg_stats(mat, outcome, covariates)
            }
        } else {
            stop(paste0("Invalid family:", family))
        }
    }
    bstats
}

bayes_lin_reg_stats <- function(mat, outcome, covariates, weights = NULL, k = 0, n_sites = 0) {
    # Compute summary statistics for Bayesian linear regression

    if (k == 0) {
        # Homogeneous setting: 1 global intercept
        mat["Intercept"] <- 1
        intercept_colnames <- "Intercept"
    } else if (k > 0 && n_sites > 0) {
        # Heterogeneous setting: 1 intercept per local center
        intercepts <- matrix(0, nrow(mat), n_sites)
        intercepts[, k] <- 1
        mat_new <- cbind(mat, intercepts)
        intercept_colnames <- paste0("Intercept_", seq_len(ncol(intercepts)))
        colnames(mat_new) <- c(colnames(mat), intercept_colnames)
        mat <- mat_new
    } else {
        stop("k must be >= 0 and n_sites must be > 0 if k > 0")
    }

    x <- mat[, c(intercept_colnames, covariates)]

    if (is.null(weights)) {
        weights <- diag(1, dim(x))
    }

    y <- mat[[outcome]]
    yy <- t(y) %*% weights %*% y

    x <- t(t(x))
    xx <- t(x) %*% weights %*% x

    xy <- t(x) %*% weights %*% y

    n <- nrow(mat) # sample size

    list(
        "xx" = xx,
        "xy" = xy,
        "yy" = yy,
        "n" = n
    )
}

bayes_lin_reg_post_params <- function(bayes_stats_list, prior_params, use_local_variances = FALSE) {
    # Compute posterior parameters for Bayesian linear regression

    mu0 <- prior_params$mu0
    lambda0 <- prior_params$lambda0
    a0 <- prior_params$a0
    b0 <- prior_params$b0

    sum_xx <- Reduce("+", lapply(bayes_stats_list, function(x) x$xx))
    sum_xy <- Reduce("+", lapply(bayes_stats_list, function(x) x$xy))
    sum_yy <- Reduce("+", lapply(bayes_stats_list, function(x) x$yy))
    sum_n <- Reduce("+", lapply(bayes_stats_list, function(x) x$n))

    lambda_l <- lambda0 + sum_xx
    beta_l <- solve(lambda_l) %*% (lambda0 %*% mu0 + sum_xy)

    if (use_local_variances) {
        # Assumes reference prior is used!
        a_l <- list()
        b_l <- list()
        for (l in seq_along(bayes_stats_list)) {
            a_l[[l]] <- bayes_stats_list[[l]]$n / 2
            xx <- bayes_stats_list[[l]]$xx
            xy <- bayes_stats_list[[l]]$xy
            beta_l_l <- solve(xx) %*% xy
            b_l[[l]] <- 0.5 * bayes_stats_list[[l]]$yy - 0.5 * t(beta_l_l) %*% xx %*% beta_l_l
        }
        a_l <- as.numeric(a_l)
        b_l <- as.numeric(b_l)
    } else {
        a_l <- a0 + sum_n / 2
        b_l <- b0 + 0.5 * sum_yy + 0.5 * (t(mu0) %*% lambda0 %*% mu0 - t(beta_l) %*% lambda_l %*% beta_l)
    }
    list(
        "lambda_l" = lambda_l,
        "beta_l" = beta_l,
        "a_l" = a_l,
        "b_l" = b_l
    )
}

bayes_lin_reg_post_map <- function(bayes_post_params, p) {
    # Compute maximum a posteriori estimates for Bayesian linear regression
    beta_l <- bayes_post_params$beta_l
    a_l <- bayes_post_params$a_l
    b_l <- bayes_post_params$b_l
    lambda_l <- bayes_post_params$lambda_l
    tau_l <- as.numeric(a_l / b_l)
    sigma_l <- solve(tau_l * lambda_l)
    list("beta_l" = beta_l, "sigma_l" = sigma_l, "lambda_l" = lambda_l)
}

get_linreg_prior <- function(covariates, use_local_intercepts, n_sites, epsilon = 1e-10) {
    # Define prior for Bayesian linear regression
    if (use_local_intercepts) {
        p <- length(covariates) + n_sites
    } else {
        p <- length(covariates)
    }

    if (use_local_intercepts) {
        prior_params <- list(
            "mu0" = as.matrix(rep(0, length(covariates) + n_sites)),
            "lambda0" = epsilon * diag(length(covariates) + n_sites),
            "a0" = epsilon - p / 2,
            "b0" = epsilon
        )
    } else {
        prior_params <- list(
            "mu0" = as.matrix(rep(0, length(covariates) + 1)),
            "lambda0" = epsilon * diag(length(covariates) + 1),
            "a0" = epsilon - p / 2,
            "b0" = epsilon
        )
    }
    prior_params
}

#' @title BCA one-shot
#'
#' @description Perform one-shot computation
#'
#' @param bstats List of summary statistics returned by bca_iterate_sites
#'
#' @return List of oneshot parameter estimates
#'
#' @author Peter Degen
#'
#' @export
bca_oneshot <- function(bstats, n_sites, use_local_intercepts, use_local_variances, family,
                           covariates, epsilon = 1e-10, center_name = NULL, CI = "normal",
                           return_post_params = FALSE, alpha=0.05) {
    # Compute BCA parameters using one-shot approach
    params_oneshot <- list()

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
        params_oneshot$sigma <- updated_params$sigma
        params_oneshot$beta <- updated_params$beta
    } else {
        print("Bayesian linear regression")
        stopifnot(!is.null(center_name))
        covariates_local <- covariates[covariates != center_name]
        prior_params <- get_linreg_prior(covariates_local, use_local_intercepts, n_sites, epsilon = epsilon)

        if (use_local_intercepts) {
            p <- length(covariates) + n_sites
        } else {
            p <- length(covariates)
        }

        bayes_post_params <- bayes_lin_reg_post_params(bstats, prior_params, use_local_variances)
        bayes_map <- bayes_lin_reg_post_map(bayes_post_params, p)
        params_oneshot$beta <- bayes_map$beta_l
        params_oneshot$sigma <- bayes_map$sigma_l
        params_oneshot$a_l <- bayes_post_params$a_l
        params_oneshot$b_l <- bayes_post_params$b_l
        params_oneshot$lambda_l <- bayes_map$lambda_l

        # use inverse of mean to get equivalent sigma2 to lm
        params_oneshot$dispersion <- as.numeric(bayes_post_params$b_l / (bayes_post_params$a_l))
        disp_ci <- invgamma::qinvgamma(c(alpha/2, 1-alpha/2),
                                       shape = bayes_post_params$a_l,
                                       rate = bayes_post_params$b_l)

        if (return_post_params) {
            params_oneshot$post_params <- bayes_post_params
        }
    }

    if (CI == "t") {
        params_oneshot$CI <- get_bayes_linreg_ci_t(params_oneshot)
    } else if (CI == "normal") {
        params_oneshot$CI <- get_bayes_linreg_ci_normal(params_oneshot)
    }

    if (family == "gaussian") {
        new_row <- data.frame(disp_ci[[1]], disp_ci[[2]], row.names = "sigma2")
        colnames(new_row) <- colnames(params_oneshot$CI)
        params_oneshot$CI <- rbind(params_oneshot$CI, new_row)
    }

    params_oneshot
}

get_bayes_linreg_ci_normal <- function(params_oneshot, alpha = 0.05) {
    # Compute credible intervals for Bayesian linear regression

    # normal
    z <- qnorm(1 - alpha / 2) # ≈ 1.96 for 95% CI
    lower <- params_oneshot$beta - z * sqrt(diag(params_oneshot$sigma))
    upper <- params_oneshot$beta + z * sqrt(diag(params_oneshot$sigma))
    ci <- data.frame(
        lower = lower,
        upper = upper
    )
    ci
}

get_bayes_linreg_ci_t <- function(params_oneshot, alpha = 0.05) {
    # Use marginal t-distributions for CI

    mu <- params_oneshot$beta
    a <- as.vector(params_oneshot$a_l)
    b <- as.vector(params_oneshot$b_l)
    lambda_diag <- diag(params_oneshot$lambda_l)

    # Degrees of freedom
    nu <- 2 * a # == n - m

    scale <- (lambda_diag * a) / b
    variance <- (1 / scale) * nu / (nu - 2) # B-S p. 435
    stdev <- sqrt(variance)

    t_crit <- qt(1 - alpha / 2, df = nu)
    lower <- mu - t_crit * stdev
    upper <- mu + t_crit * stdev

    data.frame(
        lower = lower,
        upper = upper
    )
}

#' @title BCA tidy results
#'
#' @description Tidy results from one-shot
#'
#' @param params_oneshot Parameters returned by bca_oneshot
#'
#' @return Dataframe
#'
#' @author Peter Degen
#'
#' @export
tidy_results <- function(params_oneshot, use_local_intercepts) {
    params_oneshot_all <- rbind(params_oneshot$beta, sigma2 = params_oneshot$dispersion)
    if (use_local_intercepts) {
        params_oneshot_all <- params_oneshot_all[!grepl("Intercept_\\d+$", rownames(params_oneshot_all)), ,
                                                drop = FALSE]
    }
    df <- data.frame(t(params_oneshot_all), check.names = FALSE) # prevent renaming (Intercept) to X.Intercept.
    df$Method <- "BCA"
    df <- df %>% pivot_longer(-Method, names_to = "Covariate", values_to = "Estimate")

    params_oneshot$CI$Method <- "BCA"
    if (!("Covariate" %in% colnames(params_oneshot$CI))) {
        params_oneshot$CI <- tibble::rownames_to_column(params_oneshot$CI, var = "Covariate")
    }

    df_merged <- left_join(df, params_oneshot$CI, by = c("Method", "Covariate"))
    df_merged
}

get_reduced_params <- function(center_identity, params_oneshot, bstats, family) {
    # Get params for reduced posterior by removing one center from the full posterior

    l <- center_identity

    if (family == "gaussian") {
        bayes_post_params <- params_oneshot$post_params
        lambda_minus_l <- bayes_post_params$lambda_l - bstats[[l]]$xx
        a_minus_l <- bayes_post_params$a_l - bstats[[l]]$n
        beta_minus_l <- solve(lambda_minus_l) %*% (bayes_post_params$lambda_l
            %*% bayes_post_params$beta_l - bstats[[l]]$xy)
        blb_full <- t(bayes_post_params$beta_l) %*% bayes_post_params$lambda_l %*% bayes_post_params$beta_l
        blb_local <- t(beta_minus_l) %*% lambda_minus_l %*% beta_minus_l
        b_minus_l <- bayes_post_params$b_l + 0.5 * (blb_full - blb_local - bstats[[l]]$yy)

        return(list(
            "lambda_minus_l" = lambda_minus_l,
            "a_minus_l" = a_minus_l,
            "beta_minus_l" = t(beta_minus_l),
            "b_minus_l" = b_minus_l
        ))
    } else if (family == "binomial") {
        delta_1s <- solve(params_oneshot$sigma)
        delta_l <- solve(bstats[[l]]$sigma)
        sigma_minus_l <- solve(delta_1s - delta_l)
        beta_minus_l <- sigma_minus_l %*% (delta_1s %*% params_oneshot$beta - delta_l %*% bstats[[l]]$beta)

        return(list(
            "sigma_minus_l" = sigma_minus_l,
            "beta_minus_l" = t(beta_minus_l)
        ))
    }
}

get_pred_probs <- function(center_identity, bstats, res_local, covariates_local,
                           n_sites, reduced_params, use_local_intercepts) {
    # Box prior predictive tail probability, check if center l estimate compatible with posterior from all other centers

    l <- center_identity
    rp <- reduced_params

    sigma2_pooled <- function(sigma2, nu, a, b) {
        (nu * sigma2 + 2 * b) / (nu + 2 * a)
    }

    fstat <- function(beta1, beta2, xx, lambda, m, sigma2, n, a, b) {
        nu <- n - m
        sp <- sigma2_pooled(sigma2, nu, a, b)
        diff <- beta1 - beta2
        as.numeric(t(diff) %*% solve(solve(xx) + solve(lambda)) %*% diff / (m * sp))
    }

    bl <- bstats[[l]]
    n_l <- bl$n
    sigma2 <- res_local$sigma2_list[[l]]
    m <- length(covariates_local)
    xx <- bl$xx

    # Assumes intercepts come first in covariate ordering
    if (use_local_intercepts) {
        xx <- xx[(n_sites + 1):(n_sites + m), (n_sites + 1):(n_sites + m)]
        xy <- bl$xy[(n_sites + 1):(n_sites + m)]
        beta_minus_l_cov <- rp$beta_minus_l[(n_sites + 1):(n_sites + m)]
        lambda_minus_l_cov <- rp$lambda_minus_l[(n_sites + 1):(n_sites + m), (n_sites + 1):(n_sites + m)]
    } else {
        xy <- bl$xy
        beta_minus_l_cov <- t(rp$beta_minus_l)
        lambda_minus_l_cov <- rp$lambda_minus_l
    }

    bhatl <- solve(xx) %*% xy
    f <- fstat(bhatl, beta_minus_l_cov, xx, lambda_minus_l_cov, m, sigma2, n_l, rp$a_minus_l, rp$b_minus_l)
    pf(f, df1 = m, df2 = (n_l - m) + 2 * rp$a_minus_l, lower.tail = FALSE)
}
