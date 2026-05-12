get_na_cols <- function(df) {
    names(df)[colSums(is.na(df)) > 0]
}

get_complete_sites <- function(df, center_name) {
    complete_rows <- complete.cases(df[, setdiff(names(df), center_name)])
    unique(df[[center_name]][complete_rows])
}

remove_covs_from_model <- function(model, covs_to_remove) {
    current_terms <- attr(terms(model), "term.labels")
    kept_terms <- setdiff(current_terms, covs_to_remove)
    response <- deparse(model[[2]])
    reformulate(termlabels = kept_terms, response = response)
}

reduce_data_split <- function(data_split, cols_to_remove) {
    lapply(data_split, function(df) {
        df[, setdiff(names(df), cols_to_remove), drop = FALSE]
    })
}

debiased_oneshot <- function(data_split,
                             model,
                             family,
                             center_name,
                             use_local_intercepts,
                             use_local_variances,
                             glm_prior_lamda = 0.1) {
    incomplete_covs <- get_na_cols(data)
    complete_sites <- get_complete_sites(data, center_name)
    model_reduced <- remove_covs_from_model(model, incomplete_covs)
    cat("Incomplete covariates:", incomplete_covs, "\n")
    cat("Complete sites:", length(complete_sites), "/", length(data_split), "\n")
    cat("Reduced model:", deparse(model_reduced))

    # eta: reduced model on all centers
    data_split_reduced <- reduce_data_split(data_split, incomplete_covs)
    sumstats_eta <- confeR::bca_iterate_sites(model_reduced, family, data_split_reduced, use_local_intercepts)
    params_oneshot_eta <- confeR::bca_oneshot(
        sumstats_eta, family, use_local_intercepts, use_local_variances,
        glm_prior_lamda
    )


    # eta_cc: reduced model on complete centers
    sumstats_eta_cc <- sumstats_eta[complete_sites]
    params_oneshot_eta_cc <- confeR::bca_oneshot(sumstats_eta_cc, family, use_local_intercepts, use_local_variances,
        glm_prior_lamda = glm_prior_lamda
    )

    # (beta_cc, alpha_cc): all covariates on complete centers
    data_split_cc <- data_split[complete_sites]
    sumstats_cc <- confeR::bca_iterate_sites(model, family, data_split_cc, use_local_intercepts)
    params_oneshot_cc <- confeR::bca_oneshot(sumstats_cc, family, use_local_intercepts, use_local_variances,
        glm_prior_lamda = glm_prior_lamda
    )

    # combine result: theta = (eta + beta_cc - eta_cc, alpha_cc)
    # note we call theta <=> beta_l (not same beta)
    params_oneshot <- list()
    x <- params_oneshot_cc$beta_l
    alpha_cc <- as.matrix(x[incomplete_covs, ])
    beta_cc <- x[!rownames(x) %in% incomplete_covs, ]
    beta <- params_oneshot_eta$beta_l + beta_cc - params_oneshot_eta_cc$beta_l
    params_oneshot$beta_l <- rbind(beta, alpha_cc)
    params_oneshot
}
