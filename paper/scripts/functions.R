library(dplyr)
library(tibble)
library(metafor)
library(Matrix)
library(BFI)
library(lme4)

### Local fits ###

fit_local_glms <- function(data_split, Method, target, covariates, center_name) {

    coef_list <- list()
    se_list <- list()
    cov_list <- list()
    sigma2_list <- list()
    covariates_local <- covariates[covariates != center_name]

    for (i in seq_along(data_split)) {
        l <- levels(data[[center_name]])[i]
        sub_x <- data[data[[center_name]] == l, c(target, covariates_local)]
        sub_fit <- glm(Method, family=family, data=sub_x)
        coef_summary <- summary(sub_fit)$coefficients
        coef_list[[i]] <- coef_summary[, "Estimate"]
        se_list[[i]] <- coef_summary[, "Std. Error"]
        cov_list[[i]] <- vcov(sub_fit)
        sigma2_list[[i]] <- sigma(sub_fit)^2
    }

    list("coef_list" = coef_list,
        "se_list" = se_list,
        "cov_list" = cov_list,
        "sigma2_list" = as.numeric(sigma2_list)
    )
}

### Meta-analysis ###

fit_mv_meta_fixed <- function(coef_list, cov_list, use_local_intercepts) {

    if (use_local_intercepts) {
        cov_list <- lapply(cov_list, function(m) m[rownames(m) != "(Intercept)", colnames(m) != "(Intercept)"])
        local_intercepts <- lapply(coef_list, function(m) m[names(m) == "(Intercept)"])
        coef_list <- lapply(coef_list, function(m) m[names(m) != "(Intercept)"])
    }

    # Stack all coefficients
    y <- do.call(c, coef_list)
    names_y <- names(y)

    # Create block-diagonal covariance matrix
    V <- Matrix::bdiag(cov_list)
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

    if (use_local_intercepts) {
        local_intercepts_row <- lapply(seq_along(local_intercepts), function(i) {
            x <- local_intercepts[[i]]
            data.frame(
                Estimate = x,
                lower = NaN, # TO DO
                upper = NaN,
                Covariate = paste0("Intercept_", i),
                Method = "FE",
                row.names = NULL
            )
        })
        results_meta_mv <- bind_rows(results_meta_mv, local_intercepts_row)
    }
    results_meta_mv
}

fit_mv_meta_random <- function(coef_list, cov_list, method="REML", use_local_intercepts) {

    stopifnot(method != "FE" && method != "EE")

    # Fixed-effects site-specific estimates for local intercepts
    if (use_local_intercepts) {
        cov_list <- lapply(cov_list, function(m) m[rownames(m) != "(Intercept)", colnames(m) != "(Intercept)"])
        local_intercepts <- lapply(coef_list, function(m) m[names(m) == "(Intercept)"])
        coef_list <- lapply(coef_list, function(m) m[names(m) != "(Intercept)"])
    }

    # Stack all coefficients
    y <- do.call(c, coef_list)
    names_y <- names(y)

    # Create block-diagonal covariance matrix
    V <- Matrix::bdiag(cov_list)
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

    if (use_local_intercepts) {
        local_intercepts_row <- lapply(seq_along(local_intercepts), function(i) {
            x <- local_intercepts[[i]]
            data.frame(
                Estimate = x,
                lower = NaN, # TO DO
                upper = NaN,
                Covariate = paste0("Intercept_", i),
                Method = method,
                row.names = NULL
            )
        })
        results_meta_mv <- bind_rows(results_meta_mv, local_intercepts_row)
    }
    results_meta_mv
}

### Combined fit ###
# glm treated as special case of glmm
fit_combined_glmm <- function(data, family, model, model_no_int,
                            center_name, return_fit=FALSE, heterogeneity_effect="fixed") {
    require(tibble)
    require(dplyr)

    if (heterogeneity_effect == "none") {
        fit_comb_glmm <- glm(model, family, data)
    } else if (heterogeneity_effect=="fixed") {
        # Controlled for local center
        m_expanded <- update(model_no_int, as.formula(paste("~ . +", center_name)))
        fit_comb_glmm <- glm(m_expanded, family, data)
    } else if (heterogeneity_effect=="random") {
        m_mixed <- update(
            model,
            as.formula(paste("~ . + (1 |", center_name, ") - 1"))
        )
        if (family == "gaussian") {
            fit_comb_glmm <- lmer(m_mixed, data = data)
        } else {
            fit_comb_glmm <- glmer(m_mixed, family = family, data = data)
        }
    } else {
        stop(paste("Invalid heterogeneity effect:", heterogeneity_effect))
    }

    if (return_fit) return(fit_comb_glmm)

    sum_fit <- summary(fit_comb_glmm)
    coefs <- sum_fit$coefficients

    if (heterogeneity_effect=="random") {
        sum_fit <- summary(fit_comb_glmm)
        coefs <- sum_fit$coefficients[, "Estimate"]
        coefs <- data.frame("Estimate"=coefs)
        re <- ranef(fit_comb_glmm)$hospital
        colnames(re) <- "Estimate"
        sigma2 <- data.frame("Estimate"=sigma(fit_comb_glmm)^2, row.names="sigma2")
        coefs <- rbind(re, sigma2, coefs)
        rownames(coefs)[seq_len(nrow(re))] <- paste0("Intercept_", seq_len(nrow(re)))

        ci <- data.frame(confint(fit_comb_glmm,  oldNames=FALSE))
        colnames(ci) <- c("lower", "upper")
        ci <- rbind(ci, data.frame(ci["sigma", ]^2, row.names="sigma2"))

        coefs <- coefs |>
        rownames_to_column("Covariate") |>
        left_join(ci |> rownames_to_column("Covariate") |> select("Covariate", "lower", "upper"), by = "Covariate")
        coefs$Method <- "Combined_RE"
        df_combined <- coefs
    } else {

        rownames(coefs) <- sub(paste0("^", center_name), "Intercept_", rownames(coefs))
        sum_fit$coefficients <- coefs

        ci <- confint(fit_comb_glmm, level = 0.95)
        df_combined <- data.frame(
            lower = ci[, 1],
            upper = ci[, 2],
            Method = "Combined"
        )

        df_combined <- rownames_to_column(df_combined, var = "Covariate")
        df_combined$Covariate <- rownames(coefs)
        df_combined$Estimate <- as.data.frame(coefs)$Estimate

        # Add dispersion
        if (family == "gaussian") {
            df_resid <- fit_comb_glmm$df.residual
            sigma2 <- sigma(fit_comb_glmm)^2
            alpha <- 0.05

            # https://www.graphpad.com/support/faq/the-confidence-interval-of-a-standard-deviation/
            lower <- sigma2 * (df_resid) / qchisq(1 - alpha / 2, df = df_resid)
            upper <- sigma2 * (df_resid) / qchisq(alpha / 2, df = df_resid)

            row <- list("sigma2", lower, upper, "Combined", sigma(fit_comb_glmm)^2)
            df_sigma2 <- as.data.frame(row, stringsAsFactors = FALSE)
            colnames(df_sigma2) <- colnames(df_combined)
            df_combined <- rbind(df_combined, df_sigma2)
        }
    }
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
        l <- levels(data[[center_name]])[i]

        sub_X[[i]] <- as.data.frame(subset(data, data[[center_name]] == l, select = covariates_local))

        sub_Lambda[[i]] <- inv.prior.cov(sub_X[[i]], lambda = 0.01,
                            L = length(data_split), family = family) # “gaussian”, “binomial”, “survival”
        sub_fit_bfi[[i]] <- MAP.estimation(y = data[[target]][data[[center_name]] == l], X = sub_X[[i]],
                                family = family, Lambda = sub_Lambda[[i]])
    }

    list("sub_X" = sub_X, "sub_Lambda" = sub_Lambda, "sub_fit_bfi" = sub_fit_bfi)
}

fit_bfi <- function(data_split, family, res_bfi_sub, use_local_intercepts,
                    use_local_variances, return_fit = FALSE, alpha=0.05, lambda=0.01) {

    if (family == "bernoulli") family <- "binomial"

    sub_fit_bfi <- res_bfi_sub$sub_fit_bfi
    sub_X <- res_bfi_sub$sub_X

    Ms <- fits <- thetahats <- Ahats <- priors <- Lambdas <- vector("list", length(data_split))
    for (l in seq_along(data_split)) {
        Ms[[l]] <- as.data.frame(data_split[[l]])
        Ms[[l]]$hospital <- as.double(as.character(Ms[[l]]$hospital))
        Lambdas[[l]] <- inv.prior.cov(Ms[[l]], lambda=lambda, family=family)
        fits[[l]] <- sub_fit_bfi[[l]]
        thetahats[[l]] <- fits[[l]]$theta_hat
        Ahats[[l]] <- fits[[l]]$A_hat
        priors[[l]] <- fits[[l]]$Lambda
    }

    strat_par <- if (use_local_intercepts) 1 else NULL
    stratified <- use_local_variances || use_local_intercepts
    if (use_local_variances && use_local_intercepts) strat_par <- 1:2
    else if (use_local_variances) strat_par <- 2

    Lambda_com <- BFI::inv.prior.cov(sub_X[[1]], lambda=lambda, L=length(data_split),
                        family=family, stratified=stratified, strat_par=strat_par)
    priors_all <- append(priors, list(Lambda_com))
    fit.bfi <- BFI::bfi(theta_hats=thetahats, A_hats=Ahats, Lambda=priors_all,
                    family=family, stratified=stratified, strat_par=strat_par)

    if (return_fit) return(fit.bfi)

    tidy_bfi_fit(fit.bfi, use_local_intercepts, use_local_variances, alpha)
}

tidy_bfi_fit <- function(fit.bfi, use_local_intercepts, use_local_variances, alpha=0.05) {

    z <- qnorm(1 - alpha / 2)
    se <- sqrt(diag(solve(fit.bfi$A_hat)))

    if (!use_local_intercepts && !use_local_variances) {
        df_bfi <- setNames(as.data.frame(t(fit.bfi$theta_hat - z * se)), "lower")
        df_bfi$upper <- t(fit.bfi$theta_hat + z * se)
        df_bfi$Estimate <- t(fit.bfi$theta_hat)
    } else {
        df_bfi <- setNames(as.data.frame(fit.bfi$theta_hat - z * se), "lower")
        df_bfi$upper <- fit.bfi$theta_hat + z * se
        df_bfi$Estimate <- fit.bfi$theta_hat
    }

    rownames(df_bfi) <- sub("^\\(Intercept\\)_loc(\\d+)$", "Intercept_\\1", rownames(df_bfi))
    rownames(df_bfi) <- sub("^sigma2_loc(\\d+)$", "sigma2_\\1", rownames(df_bfi))

    df_bfi$Method <- "BFI"
    df_bfi$Covariate <- rownames(df_bfi)
    rownames(df_bfi) <- NULL

    df_bfi
}

### PDA ###

library(pda)

fit_pda <- function(target, covariates,
                    data_split, sites, dataname, family, dir, heterogeneity_effect="fixed") {

    model_pda <- if (family == "gaussian") "DLM" else "ODAL"

    # ############################  STEP 1: initialize  ###############################
    ## lead site1: please review and enter "1" to allow putting the control file to the server
    control <- list(project_name = dataname,
                    step = "initialize",
                    sites = sites,
                    heterogeneity = heterogeneity_effect != "none",
                    heterogeneity_effect = heterogeneity_effect, # no effect if heterogeneity false
                    model = model_pda,
                    family = family,
                    outcome = target,
                    variables = covariates,
                    optim_maxit = 100,
                    lead_site = sites[[1]],
                    upload_date = as.character(Sys.time()))

    ## run the example in local directory:
    ## assume lead site1: enter "1" to allow transferring the control file
    pda(site_id = sites[[1]], control = control, dir = dir)

    for (site in rev(sites)) {
        ##" assume remote site l: enter "1" to allow tranferring your local estimate
        pda(site_id = site, ipdata = data_split[[as.numeric(site)]], dir=dir)
    }

    #' ############################'  STEP 2: derivative  ###############################

    for (site in rev(sites)) {
        ##' assume remote site l: enter "1" to allow tranferring your derivatives
        pda(site_id = site, ipdata =  data_split[[as.numeric(site)]], dir=dir)
    }

    #' ############################'  STEP 3: estimate  ###############################
    ##' assume lead site1: enter "1" to allow tranferring the surrogate estimate
    pda(site_id = sites[[1]], ipdata = data_split[[1]], dir=dir)

    ##' the PDA is now completed!
    ##' All the sites can still run their own surrogate estimates and broadcast them.

    config <- getCloudConfig(site_id = sites[[1]], dir=dir)
    pdaGet(name = "1_estimate", config = config)
}


tidy_pda <- function(fit.pda, family, covariates, n_sites, n_data, heterogeneity_effect, alpha=0.05) {

    z <- qnorm(1 - alpha / 2)

    if (family == "gaussian") {

        if (heterogeneity_effect != "none") {
            # PDA FE always fits with global intercept and L-1 local intercepts
            # Therefore, use global intercept as intercept 1 and add to remaining intercepts
            if (heterogeneity_effect == "fixed") {
                fit.clean <- list()
                fit.clean$sigmahat <- fit.pda$sigmahat
                fit.clean$risk_factor <- c(
                    paste0("Intercept_", seq_len(n_sites)),
                    covariates
                )
                fit.clean$bhat <- c(fit.pda$bhat[[1]], # global intercept == intercept 1
                                    fit.pda$bhat[[1]] + fit.pda$uhat, # remaining intercepts
                                    fit.pda$bhat[2:length(fit.pda$bhat)] # covariate effects
                                    )
                fit.clean$sebhat <- c(fit.pda$sebhat[[1]], fit.pda$seuhat, fit.pda$sebhat[2:length(fit.pda$bhat)])
                fit.pda <- fit.clean

            # PDA RE always fits with global intercept and L local intercepts
            } else if (heterogeneity_effect == "random") {
                fit.clean <- list()
                fit.clean$sigmahat <- fit.pda$sigmahat
                fit.clean$risk_factor <- c(
                    paste0("Intercept_", seq_len(n_sites)),
                    covariates
                )
                fit.clean$bhat <- c(fit.pda$bhat[[1]] + fit.pda$uhat,
                                    fit.pda$bhat[2:length(fit.pda$bhat)]
                                    )
                print("Warning: ignoring global intercept sebhat (to fix)")
                fit.clean$sebhat <- c(fit.pda$seuhat, fit.pda$sebhat[2:length(fit.pda$bhat)])
                fit.pda <- fit.clean
            } else {
                stop(paste("Invalid heterogeneity effect:", heterogeneity_effect))
            }
        }

        lower <- fit.pda$bhat - z * fit.pda$sebhat
        upper <- fit.pda$bhat + z * fit.pda$sebhat
        df_pda <- data.frame(
            lower = lower,
            upper = upper
        )
        df_pda$Estimate <- fit.pda$bhat
        method <- if (heterogeneity_effect=="random") "PDA_RE" else "PDA_FE"
        df_pda$Method <- method
        df_pda$Covariate <- fit.pda$risk_factor
        row <- list(NA, NA, fit.pda$sigmahat^2, method, "sigma2")
        df_sigma2 <- as.data.frame(row, stringsAsFactors = FALSE)
        colnames(df_sigma2) <- colnames(df_pda)
        df_pda <- rbind(df_pda, df_sigma2)

    } else {
        # https://github.com/Penncil/pda/blob/master/R/ODAL.R
        # setilde=sqrt(diag(solve(sol$hessian))/N
        vcov <- solve(fit.pda$Htilde)
        se <- sqrt(diag(vcov) / n_data)

        # Wald confidence intervals
        lower <- fit.pda$btilde - z * se
        upper <- fit.pda$btilde + z * se

        df_pda <- data.frame(
            lower = lower,
            upper = upper
        )
        df_pda$Estimate <- fit.pda$btilde
        df_pda$Method <- "PDA_FE"

        if (!is.null(fit.pda$risk_factor))
            df_pda$Covariate <- fit.pda$risk_factor
        else
            df_pda$Covariate <- c("(Intercept)", covariates_local)
    }

    df_pda
}
