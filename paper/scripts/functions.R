library(metafor)
library(Matrix)
library(BFI)

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

fit_mv_meta_fixed <- function(coef_list, cov_list) {

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
    results_meta_mv
}

fit_mv_meta_random <- function(coef_list, cov_list, method="REML") {

    stopifnot(method != "FE" && method != "EE")

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
    results_meta_mv
}

### Combined fit ###

fit_combined_glm <- function(model, model_no_int, use_local_intercepts, center_name, return_fit=FALSE) {
    if (use_local_intercepts) {
        # Controlled for local center
        m_expanded <- update(model_no_int, as.formula(paste("~ . +", center_name)))
        fit_comb_glm <- glm(m_expanded, family, data)
    } else {
        fit_comb_glm <- glm(model, family, data)
    }

    if (return_fit) return(fit_comb_glm)

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
    if (family == "gaussian") {
        df_resid <- fit_comb_glm$df.residual
        sigma2 <- sigma(fit_comb_glm)^2
        alpha <- 0.05

        # https://www.graphpad.com/support/faq/the-confidence-interval-of-a-standard-deviation/
        lower <- sigma2 * (df_resid) / qchisq(1 - alpha / 2, df = df_resid)
        upper <- sigma2 * (df_resid) / qchisq(alpha / 2, df = df_resid)

        row <- list("sigma2", lower, upper, "Combined", sigma(fit_comb_glm)^2)
        df_sigma2 <- as.data.frame(row, stringsAsFactors = FALSE)
        colnames(df_sigma2) <- colnames(df_combined)
        df_combined <- rbind(df_combined, df_sigma2)
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

fit_bfi <- function(data_split, family, res_bfi_sub, use_local_intercepts, return_fit = FALSE) {

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
        BFI_fit <- bfi(theta_hats=thetahats, A_hats=Ahats, Lambda=priors_all,
                        family=family, stratified=TRUE, strat_par=1)
    } else {
        Lambda_com <- inv.prior.cov(sub_X[[1]], lambda=0.01, L=length(data_split), family=family)
        priors_all <- append(priors, list(Lambda_com))
        BFI_fit <- bfi(theta_hats=thetahats, A_hats=Ahats, Lambda=priors_all, family=family)
    }

    if (return_fit) return(return_fit)

    tidy_bfi_fit(BFI_fit, use_local_intercepts, family)
}

tidy_bfi_fit <- function(BFI_fit, use_local_intercepts, family) {
    df_bfi <- as.data.frame(summary(BFI_fit)$CI)
    rownames(df_bfi) <- sub("^\\(Intercept\\)_loc(\\d+)$", "Intercept_\\1", rownames(df_bfi))
    colnames(df_bfi)[colnames(df_bfi) == "2.5 %"] <- "lower"
    colnames(df_bfi)[colnames(df_bfi) == " 97.5 %"] <- "upper"

    df_bfi$Method <- "BFI"
    df_bfi$Covariate <- rownames(df_bfi)
    rownames(df_bfi) <- NULL

    # add dispersion
    if (family == "gaussian") {

        df_bfi$Estimate <- BFI_fit$theta_hat[-length(BFI_fit$theta_hat)]

        if (use_local_intercepts) {
            sigma2 <- BFI_fit$theta_hat[["sigma2"]]
        } else {
            sigma2 <- as.data.frame(BFI_fit$theta_hat)$sigma2
        }
        row <- list(NA_real_, NA_real_, "BFI", "sigma2", sigma2)
        df_sigma2 <- as.data.frame(row, stringsAsFactors = FALSE)
        colnames(df_sigma2) <- colnames(df_bfi)
        df_bfi <- rbind(df_bfi, df_sigma2)
    } else {
        df_bfi$Estimate <- t(BFI_fit$theta_hat)
    }

    df_bfi
}

### PDA ###

library(pda)

fit_pda <- function(target, covariates, use_local_intercepts, data_split, sites, dataname, family, dir) {

    model_pda <- if (family == "gaussian") "DLM" else "ODAL"

    # ############################  STEP 1: initialize  ###############################
    ## lead site1: please review and enter "1" to allow putting the control file to the server
    control <- list(project_name = dataname,
                    step = "initialize",
                    sites = sites,
                    heterogeneity = use_local_intercepts,
                    heterogeneity_effect = "fixed", # if (use_local_intercepts) "random" else "fixed",
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


tidy_pda <- function(fit.pda, family, use_local_intercepts, covariates_local, n_centers, alpha=0.05) {

    z <- qnorm(1 - alpha / 2)

    if (family == "gaussian") {

        # PDA FE always fits with global intercept and L-1 local intercepts
        # Therefore, use global intercept as intercept 1 and add to remaining intercepts
        if (use_local_intercepts) {
            fit.clean <- list()
            fit.clean$sigmahat <- fit.pda$sigmahat
            fit.clean$risk_factor <- c(
                paste0("Intercept_", seq_len(n_centers)),
                covariates_local
            )
            fit.clean$bhat <- c(fit.pda$bhat[[1]],
                                fit.pda$bhat[[1]] + fit.pda$uhat,
                                fit.pda$bhat[2:length(fit.pda$bhat)]
                                )
            fit.clean$sebhat <- c(fit.pda$sebhat[[1]], fit.pda$seuhat, fit.pda$sebhat[2:length(fit.pda$bhat)])
            fit.pda <- fit.clean
        }

        lower <- fit.pda$bhat - z * fit.pda$sebhat
        upper <- fit.pda$bhat + z * fit.pda$sebhat
        df_pda <- data.frame(
            lower = lower,
            upper = upper
        )
        df_pda$Estimate <- fit.pda$bhat
        df_pda$Method <- "PDA"
        df_pda$Covariate <- fit.pda$risk_factor
        row <- list(NA, NA, fit.pda$sigmahat, "PDA", "sigma2")
        df_sigma2 <- as.data.frame(row, stringsAsFactors = FALSE)
        colnames(df_sigma2) <- colnames(df_pda)
        df_pda <- rbind(df_pda, df_sigma2)

    } else {
        #vcov_mat <- solve(fit.pda$Htilde)
        se <- sqrt(diag(fit.pda$Htilde))

        # Wald confidence intervals
        lower <- fit.pda$btilde - z * se
        upper <- fit.pda$btilde + z * se

        df_pda <- data.frame(
            lower = lower,
            upper = upper
        )
        df_pda$Estimate <- fit.pda$btilde
        df_pda$Method <- "PDA"
        df_pda$Covariate <- fit.pda$risk_factor
    }

    df_pda
}
