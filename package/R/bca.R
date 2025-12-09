library(dplyr)
library(tidyr)
library(tibble)
library(invgamma)

#' @title BCA iterate over local sites
#'
#' @description Pre-processing wrapper for testing BCA methods using individual participant data.
#'
#' @param outcome Name of outcome y.
#' @param covariates Vector of covariate names.
#' @param model Formula object relating outcome to covariates.
#' @param family Family to pass to glm() call, e.g. "gaussian".
#' @param data_split List, each element containing a dataframe with the IPD
#'   from a local site.
#' @param use_local_intercepts Logical. If true, use fixed site-specific
#'   intercepts for each local site.
#' @param center_name Character (optional). Name of covariate denoting site
#'   identity. Must be supplied if `use_local_intercepts` is TRUE.
#' @param alpha Numeric. Significance level for credible intervals. Default
#'   0.05.
#' @return List with transmitted summary statistics for each local site.
#'
#' @author Peter Degen
#'
#' @export
bca_iterate_sites <- function(model, family, data_split,
                              use_local_intercepts,
                              outcome=NULL,
                              covariates=NULL,
                              alpha=0.05
                              ) {

    if (is.null(outcome))
        outcome <- all.vars(model)[1]
    if (is.null(covariates))
        covariates <- attr(terms(model), "term.labels")

    # input checks
    if (family != "gaussian" && use_local_intercepts && family != "gaussian_forced_glm")
        stop("`family` must be `gaussian` when `use_local_intercepts` is TRUE")

    n_sites <- length(data_split)
    sumstats <- vector("list", n_sites)

    for (i in seq_along(data_split)) {

        # For GLM regression, we locally compute GLM parameters
        if (family != "gaussian" || family == "gaussian_forced_glm") {
            if (use_local_intercepts) {
                stop("Not yet implemented")
            } else {
                family_alt <- if (family == "gaussian_forced_glm") "gaussian" else family
                res <- stats::glm(model, family_alt, data_split[[i]])
                sumstats[[i]] <- list("beta" = res$coefficients, "sigma" = stats::vcov(res))

                if (family_alt %in% c("gaussian", "Gamma", "inverse.gaussian")) {
                    sumstats[[i]]$sigma2 <- summary(res)$dispersion
                    sigma2 <- summary(res)$dispersion
                    df <- res$df.residual
                    lower <- df * sigma2 / qchisq(1 - alpha/2, df)
                    upper <- df * sigma2 / qchisq(alpha/2, df)
                    sumstats[[i]]$sigma2 <- c(lower = lower, estimate = sigma2, upper = upper)
                }
            }

        # For Bayesian linear regression, we locally compute sufficient statistics
        } else if (family == "gaussian") {
            mat <- data_split[[i]]

            if (use_local_intercepts) {
                sumstats[[i]] <- bayes_lin_reg_stats(mat, outcome, covariates, k = i, n_sites = n_sites)
            } else {
                sumstats[[i]] <- bayes_lin_reg_stats(mat, outcome, covariates)
            }
        }
    }
    sumstats
}

#' @title BCA linear regression summary stats
#'
#' @description Takes in dataframe with (X, y) and returns sufficient statistics
#'   (X'X, X'y, y'y, n)
#'
#' @param df Dataframe containing covariates X and outcome y.
#' @param outcome Name of outcome y.
#' @param covariates Vector of covariate names.
#' @param weights Vector of weights for weighted linear regression. Default
#'   NULL.
#' @param k Site number in case local intercepts are used. Default 0 corresponds
#'   to not using local intercepts.
#' @param n_sites Number of sites. Default 0.
#'
#' @return List with transmitted summary statistics for each local site.
#'
#' @author Peter Degen
#'
#' @export
bayes_lin_reg_stats <- function(df, outcome, covariates, weights = NULL, k = 0, n_sites = 0) {
    # Compute summary statistics for Bayesian linear regression

    if (k == 0) {
        # Homogeneous setting: 1 global intercept
        df["Intercept"] <- 1
        intercept_colnames <- "Intercept"
    } else if (k > 0 && n_sites > 0) {
        # Heterogeneous setting: 1 intercept per local center
        intercepts <- matrix(0, nrow(df), n_sites)
        intercepts[, k] <- 1
        df_new <- cbind(df, intercepts)
        intercept_colnames <- paste0("Intercept_", seq_len(ncol(intercepts)))
        colnames(df_new) <- c(colnames(df), intercept_colnames)
        df <- df_new
    } else {
        stop("k must be >= 0 and n_sites must be > 0 if k > 0")
    }

    x <- df[, c(intercept_colnames, covariates)]
    x <- transform_X(x)

    if (is.null(weights)) {
        weights <- diag(1, dim(x))
    }

    y <- df[[outcome]]
    yy <- t(y) %*% weights %*% y

    x <- t(t(x))
    xx <- t(x) %*% weights %*% x

    xy <- t(x) %*% weights %*% y

    n <- nrow(df) # sample size

    list(
        "xx" = xx,
        "xy" = xy,
        "yy" = yy,
        "n" = n
    )
}

#' @export
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
        # TO DO: Assumes reference prior is used!
        a_l <- list()
        b_l <- list()
        for (l in seq_along(bayes_stats_list)) {
            a_l[[l]] <- a0 + bayes_stats_list[[l]]$n/2
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

#' @export
bayes_lin_reg_post_map <- function(params_oneshot) {
    # Compute maximum a posteriori estimates for Bayesian linear regression
    beta_l <- params_oneshot$beta_l
    a_l <- params_oneshot$a_l
    b_l <- params_oneshot$b_l
    lambda_l <- params_oneshot$lambda_l
    tau_l <- as.numeric(a_l / b_l)
    sigma_l <- solve(tau_l * lambda_l)
    list("beta_l" = beta_l, "sigma_l" = sigma_l, "lambda_l" = lambda_l)
}

#' @export
get_linreg_prior <- function(covariates, use_local_intercepts, n_sites, epsilon = 1e-10) {
    # Define prior for Bayesian linear regression
    if (use_local_intercepts) {
        p <- length(covariates) + n_sites
    } else {
        p <- length(covariates) + 1
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
#' @description Perform one-shot computation of parameter estimates given
#'   transmitted summary statistics.
#'
#' @param sumstats List of summary statistics from each local site, e.g. obbject
#'   returned by `bca_iterate_sites`.
#' @param family Family to pass to glm() call, e.g. "gaussian".
#' @param use_local_intercepts Logical. If true, use fixed site-specific
#'   intercepts for each local site.
#' @param use_local_variances Logical. If true, use fixed site-specific residual
#'   variances for each local site. Default FALSE.
#' @param epsilon Numeric. Regularization for prior hyperparameters. Default:
#'   1e-10.
#' @param center_name Character (optional). Name of covariate denoting site
#'   identity. Must be supplied if `family` is `gaussian`.
#' @param alpha Numeric. Significance level for credible intervals. Default:
#'   0.05.
#' @param glm_prior_lamda Numeric (optional). Prior for diagonal precision
#'   matrix in normal model, known variance case. See also Equation 7 in Jonker,
#'   Pazira and Coolen (2024). Default: 0.
#'
#' @return List of oneshot parameter estimates
#'
#' @author Peter Degen
#'
#' @export
bca_oneshot <- function(sumstats, family, use_local_intercepts, use_local_variances=FALSE,
                        epsilon = 1e-10, alpha=0.05, covariates=NULL, glm_prior_lamda=0) {

    params_oneshot <- list()
    n_sites <- length(sumstats)

    if (is.null(covariates)) {
        if (family == "gaussian") {
            ##print("Warning: Covariate names not supplied, trying to infer") # TO DO: print at higher log level
            covariates <- rownames(sumstats[[1]]$xy)
        } else {
            covariates <- rownames(sumstats[[1]]$sigma)
        }
        covariates <- covariates[!(covariates %in% c("Intercept", "(Intercept)"))]
        covariates <- covariates[!startsWith(covariates, "Intercept")]
    }

    if (family != "gaussian" || family == "gaussian_forced_glm") {
        m <- nrow(sumstats[[1]]$sigma)
        l <- diag(glm_prior_lamda, m, m)
        update_normal_known_variance <- function(beta, sigma) {
            sigma_post <- solve(Reduce(`+`, lapply(sigma, solve)) + l - n_sites*l)
            beta_post <- sigma_post %*% Reduce(`+`, Map(function(s, b) solve(s) %*% b, sigma, beta))
            list("sigma" = sigma_post, "beta" = beta_post)
        }
        beta <- lapply(sumstats, function(x) unlist(unname(x["beta"])))
        sigma <- lapply(sumstats, function(x) x[["sigma"]])
        updated_params <- update_normal_known_variance(beta, sigma)
        params_oneshot$sigma_l <- updated_params$sigma
        params_oneshot$beta_l <- updated_params$beta
    } else {
        prior_params <- get_linreg_prior(covariates, use_local_intercepts, n_sites, epsilon = epsilon)
        params_oneshot <- bayes_lin_reg_post_params(sumstats, prior_params, use_local_variances)
        bayes_map <- bayes_lin_reg_post_map(params_oneshot)
        params_oneshot$beta_l <- bayes_map$beta_l
        params_oneshot$sigma_l <- bayes_map$sigma_l
        params_oneshot$a_l <- params_oneshot$a_l
        params_oneshot$b_l <- params_oneshot$b_l
        params_oneshot$lambda_l <- bayes_map$lambda_l

        # use inverse of mean to get equivalent sigma2 to lm
        params_oneshot$dispersion <- as.numeric(params_oneshot$b_l / (params_oneshot$a_l))
        disp_ci <- invgamma::qinvgamma(c(alpha/2, 1-alpha/2),
                                       shape = params_oneshot$a_l,
                                       rate = params_oneshot$b_l)
    }

    params_oneshot$CI <- get_bayes_ci_normal(params_oneshot, alpha)

    if (family == "gaussian") {
        new_row <- data.frame(disp_ci[[1]], disp_ci[[2]], row.names = "sigma2")
        colnames(new_row) <- colnames(params_oneshot$CI)
        params_oneshot$CI <- rbind(params_oneshot$CI, new_row)
    } else if (family %in% c("gaussian_forced_glm", "Gamma", "inverse.gaussian")) {
        print("TO DO: sigma2 for glm")
    }
    params_oneshot
}

#' @export
get_bayes_ci_normal <- function(params_oneshot, alpha = 0.05) {
    # Compute credible intervals
    z <- stats::qnorm(1 - alpha / 2) # ≈ 1.96 for 95% CI
    lower <- params_oneshot$beta_l - z * sqrt(diag(params_oneshot$sigma_l))
    upper <- params_oneshot$beta_l + z * sqrt(diag(params_oneshot$sigma_l))
    ci <- data.frame(
        lower = lower,
        upper = upper
    )
    ci
}

#' @title BCA tidy results
#'
#' @description Tidy results from one-shot
#'
#' @param params_oneshot Parameters returned by `bca_oneshot`.
#' @param use_local_intercepts Logical. If true, use fixed site-specific
#'   intercepts for each local site.
#'
#' @return Dataframe
#'
#' @author Peter Degen
#'
#' @export
tidy_results <- function(params_oneshot, use_local_intercepts) {
    params_oneshot_all <- rbind(params_oneshot$beta_l, sigma2 = params_oneshot$dispersion)
    df <- data.frame(t(params_oneshot_all), check.names = FALSE) # prevent renaming (Intercept) to X.Intercept.
    df$Method <- "BCA"
    df <- df |> tidyr::pivot_longer(-.data$Method, names_to = "Covariate",
                                    values_to = "Estimate")

    params_oneshot$CI$Method <- "BCA"
    if (!("Covariate" %in% colnames(params_oneshot$CI))) {
        params_oneshot$CI <- tibble::rownames_to_column(params_oneshot$CI, var = "Covariate")
    }

    df_merged <- dplyr::left_join(df, params_oneshot$CI, by = c("Method", "Covariate"))
    df_merged["Covariate"][df_merged["Covariate"] == "Intercept"] <- "(Intercept)"
    df_merged
}

#' @title Get reduced params
#'
#' @description Use Reverse-Bayes to get reduced parameters obtained by leaving
#'   out one local site from the full posterior.
#'
#' @param center_identity The id of the site to be removed.
#' @param params_oneshot Parameters returned by bca_oneshot.
#' @param sumstats List of summary statistics from each local site, e.g. obbject
#'   returned by `bca_iterate_sites`.
#' @param family Family to pass to glm() call, e.g. "gaussian".
#'
#' @return List with reduced parameter estimates.
#'
#' @author Peter Degen
#'
#' @export
get_reduced_params <- function(center_identity, params_oneshot, sumstats, family) {

    l <- center_identity

    if (family == "gaussian") {
        lambda_minus_l <- params_oneshot$lambda_l - sumstats[[l]]$xx
        a_minus_l <- params_oneshot$a_l - sumstats[[l]]$n/2
        beta_minus_l <- solve(lambda_minus_l) %*% (params_oneshot$lambda_l
            %*% params_oneshot$beta_l - sumstats[[l]]$xy)
        blb_full <- t(params_oneshot$beta_l) %*% params_oneshot$lambda_l %*% params_oneshot$beta_l
        blb_local <- t(beta_minus_l) %*% lambda_minus_l %*% beta_minus_l
        b_minus_l <- params_oneshot$b_l + 0.5 * (blb_full - blb_local - sumstats[[l]]$yy)

        return(list( # nolint: return_linter.
            "lambda_minus_l" = lambda_minus_l,
            "a_minus_l" = a_minus_l,
            "beta_minus_l" = beta_minus_l,
            "b_minus_l" = b_minus_l
        ))
    } else {
        delta_1s <- solve(params_oneshot$sigma_l)
        delta_l <- solve(sumstats[[l]]$sigma)
        sigma_minus_l <- solve(delta_1s - delta_l)
        beta_minus_l <- sigma_minus_l %*% (delta_1s %*% params_oneshot$beta_l - delta_l %*% sumstats[[l]]$beta)

        return(list( # nolint: return_linter.
            "sigma_minus_l" = sigma_minus_l,
            "beta_minus_l" = beta_minus_l
        ))
    }
}

#' @title Box prior predictive tail probability
#'
#' @description Check if estimate from specific site is compatible with
#'   posterior from all other centers.
#'
#' @param center_identity The id of the site to be checked
#' @param sumstats List of summary statistics from each local site, e.g. obbject
#'   returned by `bca_iterate_sites`.
#' @param covariates Vector of covariate names.
#' @param n_sites Number of sites.
#' @param reduced_params Reduced posterior parameters returnd by
#'   `get_reduced_params`.
#' @param use_local_intercepts Logical. If true, use fixed site-specific
#'   intercepts for each local site.
#' @param remove_intercept Logical. If true, ignore intercept in pbox
#'   calculation. Default: FALSE.
#'
#' @return The Box prior predictive tail probability.
#'
#' @author Peter Degen
#'
#' @export
get_pred_prob <- function(center_identity, sumstats,
                           n_sites, reduced_params, use_local_intercepts, remove_intercept=FALSE) {

    if (use_local_intercepts && !remove_intercept) {
        print("Warning: removing intercept from pbox calculation when use_local_intercepts")
        remove_intercept <- TRUE
    }

    l <- center_identity
    rp <- reduced_params

    sigma2_pooled <- function(sigma2, nu, a, b) {
        (nu * sigma2 + 2 * b) / (nu + 2 * a)
    }

    fstat <- function(beta1, beta2, xx, lambda, nu, m_tilde, sigma2, n, a, b) {
        sp <- sigma2_pooled(sigma2, nu, a, b)
        diff <- beta1 - beta2
        as.numeric(t(diff) %*% solve(solve(xx) + solve(lambda)) %*% diff / (m_tilde * sp))
    }

    m <- nrow(sumstats[[1]]$xx)
    if (use_local_intercepts) {
        m <- m - n_sites  # subtract local intercepts
    } else {
        m <- m - 1  # subtract global intercept
    }

    bl <- sumstats[[l]]
    n_l <- bl$n

    m_tilde <- if (remove_intercept) m else m + 1  # = k in Box 1980
    xx <- bl$xx
    xy <- bl$xy
    yy <- bl$yy

    if (use_local_intercepts) {
        xxl <- sumstats[[l]]$xx
        last_n <- (nrow(xxl) - m + 1):nrow(xxl)
        keep <- sort(unique(c(l, last_n)))
        xxl <- xxl[keep, keep]
        xyl <- xy[keep]
        bhatl <- solve(xxl) %*% xyl
        sigma2 <- as.numeric((yy - t(bhatl)%*%xxl%*%bhatl) / (n_l - m - 1))
        # remove intercept
        bhatl <- bhatl[2:nrow(bhatl), , drop = FALSE]
    } else {
        bhatl <- solve(xx) %*% xy
        sigma2 <- as.numeric((yy - t(bhatl)%*%xx%*%bhatl) / (n_l - m - 1))
    }


    # Assumes intercepts come first in covariate ordering
    if (use_local_intercepts) {
        # double check this
        xx <- xx[(n_sites + 1):(n_sites + m), (n_sites + 1):(n_sites + m)]
        xy <- bl$xy[(n_sites + 1):(n_sites + m)]
        beta_minus_l_cov <- rp$beta_minus_l[(n_sites + 1):(n_sites + m)]
        lambda_minus_l_cov <- rp$lambda_minus_l[(n_sites + 1):(n_sites + m), (n_sites + 1):(n_sites + m)]
    } else {
        if (remove_intercept) {
            xx <- xx[1:m+1, 1:m+1]
            xy <- bl$xy[1:m+1]
            bhatl <- bhatl[1:m+1]
            beta_minus_l_cov <- rp$beta_minus_l[1:m+1]
            lambda_minus_l_cov <- rp$lambda_minus_l[1:m+1, 1:m+1]
        } else {
            xy <- bl$xy
            beta_minus_l_cov <- rp$beta_minus_l
            lambda_minus_l_cov <- rp$lambda_minus_l
        }
    }

    nu <- if (use_local_intercepts) n_l - m else n_l - m - 1
    f <- fstat(bhatl, beta_minus_l_cov, xx, lambda_minus_l_cov, nu, m_tilde, sigma2, n_l, rp$a_minus_l, rp$b_minus_l)
    stats::pf(f, df1 = m_tilde, df2 = (n_l - m - 1) + 2 * rp$a_minus_l, lower.tail = FALSE)
}
#' @title Box prior predictive tail probabilities for all sites
#'
#' @description Calculate p_Box for all sites.
#'
#' @param data_split The data split by local sites.
#' @param params_oneshot Parameters returned by `bca_oneshot`.
#' @param sumstats List of summary statistics from each local site, e.g. obbject
#'   returned by `bca_iterate_sites`.
#' @param family Family to pass to glm() call, e.g. "gaussian".
#' @param covariates Vector of covariate names.
#' @param n_sites Number of sites.
#' @param use_local_intercepts Logical. If true, use fixed site-specific
#'   intercepts for each local site.
#' @param center_name Character. Name of covariate denoting site identity.
#' @param remove_intercept Logical. If true, ignore intercept in pbox
#'   calculation. Default: FALSE.
#'
#' @return List with pboxes and reduced parameters for all sites.
#'
#' @author Peter Degen
#'
#' @export
box_check_all_sites <- function(params_oneshot, sumstats, family,
                                use_local_intercepts, center_name=NULL, remove_intercept=FALSE) {

    n_sites <- length(sumstats)
    results <- lapply(seq_along(sumstats), function(l) {
        reduced_params_l <- confeR::get_reduced_params(l, params_oneshot, sumstats, family)
        pbox <- get_pred_prob(
            l, sumstats, n_sites,
            reduced_params_l, use_local_intercepts, remove_intercept
        )
        list(
            pbox = pbox,
            beta_minus_l = t(reduced_params_l$beta_minus_l)
        )
    })

    reduced_params <- do.call(rbind, lapply(results, `[[`, "beta_minus_l"))
    pboxes <- lapply(results, `[[`, "pbox")
    pboxes <- unlist(pboxes)

    reduced_params <- as.data.frame(reduced_params)
    if (is.null(center_name))
        center_name <- "Site"
    reduced_params[[center_name]] <- seq_along(sumstats)
    list("pboxes" = pboxes, "reduced_params"=reduced_params)
}