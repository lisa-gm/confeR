library(metafor)
library(Matrix)
library(BFI)

### Meta ###

fit_local_glms <- function(data_split, Method, target, covariates, center_name) {

    coef_list <- list()
    se_list <- list()
    cov_list <- list()
    covariates_local <- covariates[covariates != center_name]

    for (i in seq_along(data_split)) {
        sub_x <- data[data[[center_name]] == i, c(target, covariates_local)]
        sub_fit <- glm(Method, family=family, data=sub_x)
        coef_summary <- summary(sub_fit)$coefficients
        coef_list[[i]] <- coef_summary[, "Estimate"]
        se_list[[i]] <- coef_summary[, "Std. Error"]
        cov_list[[i]] <- vcov(sub_fit)
    }

    list("coef_list" = coef_list,
        "se_list" = se_list,
        "cov_list" = cov_list
    )
}

fit_mv_meta_fixed <- function(coef_list, cov_list) {

    # Stack all coefficients
    y <- do.call(c, coef_list)
    names_y <- names(y)

    # Create block-diagonal covariance matrix
    V <- bdiag(cov_list)  # from Matrix package
    V <- as.matrix(V)

    # Design matrix X (maps each coefficient to its covariate)
    covariates_ <- unique(names_y)  # Preserves order of appearance
    X <- sapply(covariates_, function(cov) as.numeric(names_y == cov))
    colnames(X) <- covariates_

    fit <- rma.mv(y, V = V, mods = ~ X - 1, method = "FE")

    results_meta_mv <- data.frame(
                            Covariate = colnames(X),
                            estimate = coef(fit),
                            lower = fit$ci.lb,
                            upper = fit$ci.ub
                        )
    rownames(results_meta_mv) <- NULL
    results_meta_mv$Method <- "FE"
    names(results_meta_mv)[names(results_meta_mv) == "estimate"] <- "Estimate"
    results_meta_mv
}

fit_mv_meta_random <- function(coef_list, cov_list, method="REML") {

    stopifnot(method != "FE" && method != "EE")

    # Stack all coefficients
    y <- do.call(c, coef_list)
    names_y <- names(y)

    # Create block-diagonal covariance matrix
    V <- bdiag(cov_list)  # from Matrix package
    V <- as.matrix(V)

    # Design matrix X (maps each coefficient to its covariate)
    covariates_ <- unique(names_y)  # Preserves order of appearance
    X <- sapply(covariates_, function(cov) as.numeric(names_y == cov))
    colnames(X) <- covariates_

    center_id <- rep(seq_along(coef_list), times = sapply(coef_list, length))
    fit <- rma.mv(y, V = V, mods = ~ X - 1, random = ~ 1 | center_id, method = method)

    results_meta_mv <- data.frame(
                            Covariate = colnames(X),
                            estimate = coef(fit),
                            lower = fit$ci.lb,
                            upper = fit$ci.ub
                        )
    rownames(results_meta_mv) <- NULL
    results_meta_mv$Method <- toupper(method)
    names(results_meta_mv)[names(results_meta_mv) == "estimate"] <- "Estimate"
    results_meta_mv
}

### Combined ###

fit_combined_glm <- function(model, model_no_int, use_local_intercepts, center_name) {
    if (use_local_intercepts) {
        # Controlled for local center
        m_expanded <- update(model_no_int, as.formula(paste("~ . +", center_name)))
        fit_comb_glm <- glm(m_expanded, family, data)
    } else {
        fit_comb_glm <- glm(model, family, data)
    }

    sum_fit <- summary(fit_comb_glm)
    coefs <- sum_fit$coefficients
    rownames(coefs) <- sub(paste0("^", center_name), "Intercept_", rownames(coefs))
    sum_fit$coefficients <- coefs

    ci <- confint(fit_comb_glm, level = 0.95)
    df_combined <- data.frame(
        lower = ci[, 1],
        upper = ci[, 2],
        Method = "Combined"
    )

    df_combined <- tibble::rownames_to_column(df_combined, var = "Covariate")
    df_combined$Covariate <- rownames(coefs)
    df_combined$Estimate <- as.data.frame(coefs)$Estimate

    # Add dispersion
    row <- list("sigma2", NA_real_, NA_real_, "Combined", sigma(fit_comb_glm)^2)
    df_sigma2 <- as.data.frame(row, stringsAsFactors = FALSE)
    colnames(df_sigma2) <- colnames(df_combined)
    df_combined <- rbind(df_combined, df_sigma2)
    df_combined
}

### BFI ###

bfi_sub <- function(data_split, family, target, covariates, center_name) {
    
    covariates_local <- covariates[covariates != center_name]
    sub_X <- vector("list", length(data_split))
    sub_Lambda <- vector("list", length(data_split))
    sub_fit_bfi <- vector("list", length(data_split))

    if (family == "bernoulli") family <- "binomial"

    for (i in seq_along(data_split)) {

        sub_X[[i]] <- as.data.frame(subset(data, data[[center_name]] == i, select = covariates_local))

        sub_Lambda[[i]] <- inv.prior.cov(sub_X[[i]], lambda = 0.01,
                            L = length(data_split), family = family) # “gaussian”, “binomial”, “survival”
        sub_fit_bfi[[i]] <- MAP.estimation(y = data[[target]][data[[center_name]] == i], X = sub_X[[i]],
                                family = family, Lambda = sub_Lambda[[i]])
    }

    list("sub_X" = sub_X, "sub_Lambda" = sub_Lambda, "sub_fit_bfi" = sub_fit_bfi)
}

bfi_fit <- function(data_split, family, res_bfi_sub, use_local_intercepts) {

    if (family == "bernoulli") family <- "binomial"

    sub_fit_bfi <- res_bfi_sub$sub_fit_bfi
    sub_X <- res_bfi_sub$sub_X

    Ms <- fits <- thetahats <- Ahats <- priors <- Lambdas <- vector("list", length(data_split))
    for (l in seq_along(data_split)) {
        Ms[[l]] <- as.data.frame(data_split[[l]])
        Ms[[l]]$hospital <- as.double(as.character(Ms[[l]]$hospital))
        Lambdas[[l]] <- inv.prior.cov(Ms[[l]], lambda=0.01, family=family)
        fits[[l]] <- sub_fit_bfi[[l]]
        thetahats[[l]] <- fits[[l]]$theta_hat
        Ahats[[l]] <- fits[[l]]$A_hat
        priors[[l]] <- fits[[l]]$Lambda
    }

    if (use_local_intercepts) {
        Lambda_com <- inv.prior.cov(sub_X[[1]], lambda=0.01, L=length(data_split),
                            family=family, stratified=TRUE, strat_par=1)
        priors_all <- append(priors, list(Lambda_com))
        BFI_fit <- bfi(theta_hats=thetahats, A_hats=Ahats, Lambda=priors_all, family=family, stratified=TRUE, strat_par=1)
    } else {
        Lambda_com <- inv.prior.cov(sub_X[[1]], lambda=0.01, L=length(data_split), family=family)
        priors_all <- append(priors, list(Lambda_com))
        BFI_fit <- bfi(theta_hats=thetahats, A_hats=Ahats, Lambda=priors_all, family=family)
    }

    tidy_bfi_fit(BFI_fit)
}

tidy_bfi_fit <- function(BFI_fit) {
    df_bfi <- as.data.frame(summary(BFI_fit)$CI)
    rownames(df_bfi) <- sub("^\\(Intercept\\)_loc(\\d+)$", "Intercept_\\1", rownames(df_bfi))
    colnames(df_bfi)[colnames(df_bfi) == "2.5 %"] <- "lower"
    colnames(df_bfi)[colnames(df_bfi) == " 97.5 %"] <- "upper"

    df_bfi$Method <- "BFI"
    df_bfi$Covariate <- rownames(df_bfi)
    rownames(df_bfi) <- NULL
    df_bfi$Estimate <- BFI_fit$theta_hat[-length(BFI_fit$theta_hat)]

    # add dispersion
    row <- list(NA_real_, NA_real_, "BFI", "sigma2", BFI_fit$theta_hat[["sigma2"]])
    df_sigma2 <- as.data.frame(row, stringsAsFactors = FALSE)
    colnames(df_sigma2) <- colnames(df_bfi)
    df_bfi <- rbind(df_bfi, df_sigma2)

    df_bfi
}