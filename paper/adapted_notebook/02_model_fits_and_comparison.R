##### slightly adapted version of paper/notebook.ipynb #######
# Script 2/3: Model fitting and method comparison
# include INLA comparison
# Purpose: run federated/combined model fits and write parameter summaries.
# Run after 01_setup_and_data_prep.R in the same R session.

if (!exists("heterogeneity") || is.null(heterogeneity)) {
    source("01_setup_and_data_prep.R", local = TRUE)
}

## NOTE: needs to be ADAPTED
folder_path_dalia = paste0(script_dir, "/../../../../../dalia-project/DALIA/examples")
print(folder_path_dalia)

# ---- Using Bayesian Conjguate Analysis ----
if (heterogeneity != Heterogeneity$RANDOM) {
    use_local_intercepts_bca <- heterogeneity == Heterogeneity$FIXED
    sumstats <- bca_iterate_sites(model, family_alt, data_split, use_local_intercepts_bca)
    params_oneshot <- bca_oneshot(sumstats, family_alt, use_local_intercepts_bca, use_local_variances, 
                                    glm_prior_lamda=0.1)  # setting 0 equivalent to FE in normal model
    df_bca <- tidy_results(params_oneshot, use_local_intercepts_bca)
}

# ---- Using BFI ----
if (heterogeneity != Heterogeneity$RANDOM) {
    use_local_intercepts_bfi <- heterogeneity == Heterogeneity$FIXED
    res_bfi_sub <- bfi_sub(data_split, family, outcome, covariates, center_name)
    df_bfi <- fit_bfi(data_split, family, res_bfi_sub, use_local_intercepts_bfi, use_local_variances, lambda=0.1)
}

# ---- Using Multivariate Meta-Analysis ----
# Fit local GLMs
local_glms <- fit_local_glms(data_split, model, outcome, covariates, center_name)
coef_list <- local_glms$coef_list
se_list <- local_glms$se_list
cov_list <- local_glms$cov_list
sigma2_list <- local_glms$sigma2_list

intercept_list = numeric(length(coef_list))
for (i in seq_along(coef_list)) {
    intercept_list[i] <- coef_list[[i]][5]
}

coef_list_with_sigma2 <- lapply(seq_along(coef_list), function(i) {c(coef_list[[i]], "sigma2"=sigma2_list[[i]])})
cov_list_with_sigma2 <- lapply(seq_along(cov_list), function(i) {
                    mat <- cov_list[[i]]
                    sigma2 <- log(sigma2_list[[i]])
                    df.res <- local_glms$df.res_list[[i]]
                    var_sigma2 <- 2 * sigma2^2 / (df.res)
                    p <- nrow(mat)
                    joint_cov <- matrix(0, nrow = p + 1, ncol = p + 1)
                    joint_cov[1:p, 1:p] <- mat
                    joint_cov[p + 1, p + 1] <- var_sigma2
                    rownames(joint_cov) <- colnames(joint_cov) <- c(rownames(mat), "sigma2")
                    joint_cov
})

use_joint_covariance <- FALSE # don't use joint, sigma2 is not normal distributed... or maybe try log transform

# Meta-analyse GLM parameters
if (heterogeneity != Heterogeneity$RANDOM) {

    use_local_intercepts_meta <- heterogeneity == Heterogeneity$FIXED

    if (!use_joint_covariance) {
        df_meta_mv_fe <- fit_mv_meta_fixed(coef_list, cov_list, use_local_intercepts_meta)
        df_meta_mv_reml <- fit_mv_meta_random(coef_list, cov_list, method="REML", use_local_intercepts_meta)
    } else {
        df_meta_mv_fe <- fit_mv_meta_fixed(coef_list_with_sigma2, cov_list_with_sigma2, use_local_intercepts_meta)
        df_meta_mv_reml <- fit_mv_meta_random(coef_list_with_sigma2, cov_list_with_sigma2, method="REML", use_local_intercepts_meta)
    }

    if (family == "gaussian" && !use_joint_covariance) {
        use_boot_CI <- FALSE
        if (use_boot_CI) {
            boot_means <- replicate(10000, mean(sample(local_glms$sigma2_list, replace = TRUE)))
            ci <- quantile(boot_means, c(0.025, 0.975))
            lower <- ci[[1]]
            upper <- ci[[2]]
        } else {
            # log transformed sigma2
            # result <- meta_sigma2(local_glms$sigma2_list, local_glms$df.res_list)
            # sigma2 <- result$pooled_sigma2 
            # lower <- result$ci_lower
            # upper <- result$ci_upper

            # weighted mean of error variances
            alpha = 0.05
            s2 <- local_glms$sigma2_list
            dfr <- sum(local_glms$df.res_list)
            #weights <- local_glms$n_list # sample size weighted FE, decent estimate
            weights <- local_glms$df.res_list # df.res weighted FE, decent estimate
            s2_var <- 2 * s2^2 / local_glms$df.res_list
            #weights <-  1 / s2_var # inverse variance weighted FE, terrible estimate, same as joint
            sigma2 <- sum(s2 * weights) / sum(weights)
            # https://www.graphpad.com/support/faq/the-confidence-interval-of-a-standard-deviation/
            lower <- sigma2 * (dfr) / qchisq(1 - alpha / 2, df = dfr)
            upper <- sigma2 * (dfr) / qchisq(alpha / 2, df = dfr)

            # simple mean of error variances (gives best estimate, bad CI)
            # sigma2 <- mean(local_glms$sigma2_list)
            # sigma2_sd <- sd(local_glms$sigma2_list)
            # n <- length(local_glms$sigma2_list)
            # error_margin <- qt(0.975, df = n - 1) * sigma2_sd / sqrt(n)
            # lower <- sigma2 - error_margin
            # upper <- sigma2 + error_margin
        }

        row <- list("sigma2", sigma2, lower, upper, "FE")
        df_sigma2 <- as.data.frame(row, stringsAsFactors = FALSE)
        colnames(df_sigma2) <- colnames(df_meta_mv_fe)
        df_meta_mv_fe <- rbind(df_meta_mv_fe, df_sigma2)

        row <- list("sigma2", sigma2, lower, upper, "REML")
        df_sigma2 <- as.data.frame(row, stringsAsFactors = FALSE)
        colnames(df_sigma2) <- colnames(df_meta_mv_fe)
        df_meta_mv_reml <- rbind(df_meta_mv_reml, df_sigma2)
    }
}

# ---- Using PDA ----
# PDA intentionally excluded. The installation didn't work for me and I didn't think it was so important now.

# -------------------- INLA call using combined data --------------------- # 
if (requireNamespace("INLA", quietly = TRUE)) {
    # Match BCA/BFI prior strength (lambda = 0.1) for fixed effects.
    inla_prior_lambda <- 0.1

    if (heterogeneity == Heterogeneity$NONE) {
        model_inla <- model
    } else if (heterogeneity == Heterogeneity$FIXED) {
        model_inla <- update(model_no_int, as.formula(paste("~ . +", center_name)))
    } else if (heterogeneity == Heterogeneity$RANDOM) {
        model_inla <- update(model, as.formula(paste("~ . + f(", center_name, ", model = \"iid\") - 1")))
    } else {
        stop("Invalid heterogeneity setting for INLA")
    }

    fit_combined_inla <- INLA::inla(
        formula = model_inla,
        family = family,
        data = data,
        control.fixed = list(
            mean = 0,
            prec = inla_prior_lambda,
            mean.intercept = 0,
            prec.intercept = inla_prior_lambda
        ),
        control.compute = list(dic = TRUE, waic = TRUE)
    )

    df_inla <- fit_combined_inla$summary.fixed |>
        tibble::rownames_to_column("Covariate") |>
        dplyr::transmute(
            lower = `0.025quant`,
            upper = `0.975quant`,
            Estimate = mean,
            Method = "INLA",
            Covariate = dplyr::if_else(
                startsWith(Covariate, center_name),
                sub(paste0("^", center_name), "Intercept_", Covariate),
                Covariate
            )
        )

    if (family == "gaussian") {
        hp <- fit_combined_inla$summary.hyperpar
        precision_idx <- grep("precision", rownames(hp), ignore.case = TRUE)

        if (length(precision_idx) > 0) {
            hp_row <- hp[precision_idx[1], , drop = FALSE]

            # INLA reports precision tau; convert directly to sigma^2 = 1 / tau.
            df_inla <- rbind(
                df_inla,
                data.frame(
                    lower = 1 / hp_row[["0.975quant"]],
                    upper =  1 / hp_row[["0.025quant"]], 
                    Estimate = 1 / hp_row[["mean"]],
                    Method = "INLA",
                    Covariate = "sigma2"
                )
            )
        } else {
            message("INLA hyperparameter precision not found; skipping INLA sigma2 extraction")
        }
    }

} else {
    warning("INLA package not available, skipping INLA fit.")
}


# ---- Using combined data (ground truth) ----

df_combined <- fit_combined_glmm(data, family, model, model_no_int, center_name, heterogeneity_effect=heterogeneity)

if (heterogeneity == Heterogeneity$FIXED) {
    df_combined_re <- fit_combined_glmm(data, family, model, model_no_int,
                                        center_name, heterogeneity_effect="random")
    int_re <- df_combined_re |> dplyr::filter(startsWith(Covariate, "Intercept_"))
    int_fe <- df_combined |> dplyr::filter(startsWith(Covariate, "Intercept_"))

    df_pda_re_fe <- data.frame(
    group     = rownames(int_fe),
    FE        = int_fe$Estimate,
    RE        = int_re$Estimate
    )

    ggplot(df_pda_re_fe, aes(x = FE, y = RE, label = group)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
    geom_point(size = 2) +
    geom_hline(yintercept=mean(int_re$Estimate)) +
    geom_vline(xintercept=mean(int_fe$Estimate)) +
    geom_text(nudge_x = 0.02, nudge_y = 0.02, size = 3) +
    labs(
        x = "Fixed-effects estimate",
        y = "Random-effects estimate",
        title = "Shrinkage: RE vs FE"
    ) +
    theme_minimal()
}


# ff <- fit_combined_glm(data, family, model, model_no_int, use_local_intercepts, center_name, return_fit=TRUE)

# X <- model.matrix(ff)
# b <- coef(ff)
# Xb <- X %*% b
# y <- data$stress

# t(y-Xb) %*% (y-Xb) / df.residual(ff)
# (t(y) %*% y - t(b) %*% t(X) %*% Xb) / df.residual(ff)
# sigma(ff)^2

# ---- Compare methods ----
## load dalia results
### make sure path is configured correctly
if(dataname == "trauma"){
    example_folder = "b_federated"
} else if(dataname == "nurses_hom"){
    example_folder = "g_federated"
} 

folder_path = file.path(folder_path_dalia, example_folder)

file_name_fed = paste0("dalia_summary_", dataname, "_global_intercept.csv")
df_dalia_fed = read.csv(file.path(folder_path, file_name_fed))

file_name_joint = paste0("dalia_summary_", dataname, "_joint.csv")
df_dalia_joint = read.csv(file.path(folder_path, file_name_joint))


df_merged <- rbind(df_combined, df_bca, df_inla,  df_dalia_joint, df_dalia_fed, df_bfi, df_meta_mv_fe, df_meta_mv_reml) #  df_dalia_joint, 
df_merged$Method <- factor(df_merged$Method, levels=unique(df_merged$Method))
levels(df_merged$Method)

library(ggplot2)
library(ggpubr)
library(stringr)

if (heterogeneity != Heterogeneity$RANDOM) {
    df_merged_all <- rbind(df_combined, df_bca, df_bfi, df_inla, df_meta_mv_fe, df_meta_mv_reml) # 
} else {
    df_merged_all <- rbind(df_combined, df_inla, df_dalia_joint, df_dalia_fed)
}

method_order <- unique(df_merged_all$Method)
df_merged_all$Method <- factor(df_merged_all$Method, levels = method_order)

if (heterogeneity != Heterogeneity$NONE) {

    # method_order <- c(method_order, "Combined_RE", "PDA_RE")
    # df_merged <- rbind(df_merged, df_combined_re, df_pda_re)

    df_merged <- df_merged_all |>
    dplyr::filter(
        Covariate != "(Intercept)",
        !str_starts(Covariate, "Intercept_")
    )
}

write.csv(df_merged_all, file.path(outpath, paste0("params.", dataname, ".csv")))

huemap <- setNames(get_palette(palette = "npg", length(levels(df_merged$Method))), levels(df_merged$Method))

plot_params <- function(dataname, df_merged) {
    require(ggplot2)
    p <- ggplot(df_merged, aes(x = Method, y = Estimate, color = Method)) +
    geom_point(position = position_dodge(width = 0.5)) +
    geom_errorbar(aes(ymin = lower, ymax = upper), 
                width = 0.2,
                position = position_dodge(width = 0.5)) +
    scale_color_manual(values = huemap) +
    facet_wrap(~ Covariate, scales = "free_y") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(y = "Estimate ± 95% CI", x = NULL)

    ggsave(paste0("figures_tmp/params.", dataname, ".png"), bg = "white")
    p
}

plot_params(dataname, df_merged)

# if (startsWith(dataname, "sim-")) {
#     simtruth <- read.csv(paste0("data/raw/", raw_dataname, ".truth.csv"))
#     simtruth
# }

# Intercepts

if (heterogeneity != Heterogeneity$NONE) {

    dfi <- df_merged_all |> dplyr::filter(startsWith(Covariate, "Intercept_") & !(Method %in% c("FE")))

    df_wide <- dfi |> dplyr::select(Method, Covariate, Estimate) |>
    tidyr::pivot_wider(
        names_from = Method,
        values_from = Estimate
    )

    # Long format for plotting y-values for non-Combined methods
    df_plot <- df_wide %>%
    tidyr::pivot_longer(
        cols = -c(Covariate, Combined),
        names_to = "Method",
        values_to = "Estimate_other"
    )

    # Scatterplot
    ggplot(df_plot, aes(x = Combined, y = Estimate_other, color = Method)) +
    geom_point(size = 2) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    labs(
        x = "Estimate (Combined)",
        y = "Estimate (Other Methods)",
        color = "Method",
        title = "Comparison of Combined vs. Other Methods by Covariate"
    ) +
    theme_bw()
}

if (family == "gaussian" && heterogeneity != Heterogeneity$NONE) {
    x <- df_combined[startsWith(df_combined$Covariate, "Intercept_"), "Estimate"]
    shapiro.test(x)
    par(bg = "white")
    qqnorm(x) + theme_minimal()
    qqline(x, col = "red")
}

