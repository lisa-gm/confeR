library(metafor)
library(Matrix)

fit_local_glms <- function(data_split, model, target, covariates, center_name) {

    coef_list <- list()
    se_list <- list()
    cov_list <- list()
    covariates_local <- covariates[covariates != center_name]

    for (i in seq_along(data_split)) {
        sub_x <- data[data[[center_name]] == i, c(target, covariates_local)]
        sub_fit <- glm(model, family=family, data=sub_x)
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
    results_meta_mv$Model <- "FE"
    names(results_meta_mv)[names(results_meta_mv) == "estimate"] <- "Value"
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

    hospital_id <- rep(seq_along(coef_list), times = sapply(coef_list, length))
    fit <- rma.mv(y, V = V, mods = ~ X - 1, random = ~ 1 | hospital_id, method = method)

    results_meta_mv <- data.frame(
                            Covariate = colnames(X),
                            estimate = coef(fit),
                            lower = fit$ci.lb,
                            upper = fit$ci.ub
                        )
    rownames(results_meta_mv) <- NULL
    results_meta_mv$Model <- toupper(method)
    names(results_meta_mv)[names(results_meta_mv) == "estimate"] <- "Value"
    results_meta_mv
}