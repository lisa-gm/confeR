library(yardstick)
library(dplyr)

make_pred_logistic <- function(data, beta, formula) {
    x <- model.matrix(formula, data)
    logits <- x %*% beta
    probs <- 1 / (1 + exp(-logits))
    pred <- ifelse(probs >= 0.5, 1, 0)
    list("pred" = pred, "probs" = as.vector(probs))
}

eval_logistic <- function(truth, pred, probs = NULL)  {
    df_eval <- tibble(
        truth = as.factor(truth),
        pred = as.factor(pred)
    )
    df_eval <- df_eval |>
    mutate(across(c(truth, pred), ~ factor(.x, levels = c(0, 1))))
    metrics <- yardstick::metric_set(accuracy, precision, recall, f_meas)(df_eval, truth = truth, estimate = pred)

    if (!is.null(probs)) {
        df_eval$probs <- probs
        roc_auc_ <- yardstick::roc_auc(df_eval, truth, probs, event_level = "second")
        metrics <- rbind(metrics, roc_auc_)
    }
    metrics
}

eval_linear <- function(truth, data, beta, formula, use_local_intercepts) {

    if (use_local_intercepts) {
        formula <- update(formula, ~ 0 + . + hospital)
        x <- model.matrix(formula, data)
        colnames(x) <- sub("^hospital(\\d+)$", "Intercept_\\1", colnames(x))
        colnames(x) <- sub("\\(Intercept\\)", "Intercept_1", colnames(x))
    } else {
        x <- model.matrix(formula, data)
    }

    x <- x[, rownames(beta)]
    eta <- as.vector(x %*% beta)

    residuals <- truth - eta
    mse <- mean(residuals^2)
    rmse <- sqrt(mse)
    mae <- mean(abs(residuals))

    ss_total <- sum((truth - mean(truth))^2)
    ss_res <- sum((residuals)^2)
    r_squared <- 1 - ss_res / ss_total

    metrics <- list(
        RMSE = rmse,
        MAE = mae,
        R2 = r_squared,
        MSE = mse
    )
    list(metrics=metrics, x=x, eta=eta)
}

do_eval <- function(beta, data_split, model, outcome, use_local_intercepts) {
    require(ggplot2)
    if (family == "binomial") {
        l <- 1
        pred_probs <- make_pred_logistic(data_split[[l]], beta[[l]], model)
        pred <- pred_probs$pred
        probs <- pred_probs$probs
        truth <- data_split[[l]][[outcome]]
        table(Predicted = pred, Actual = truth)

        eval_logistic(truth, pred, probs)

    } else if (family == "gaussian") {
        res <- eval_linear(data[[outcome]], data, beta, model, use_local_intercepts)
        metrics <- res$metrics
        x <- res$x
        eta <- res$eta

        p <- ggplot(data, aes(x = .data[[outcome]], y = eta)) +
        geom_point() +
        # Identity line (y = x)
        geom_abline(slope = 1, intercept = 0, color = "grey", linetype = "dashed") +
        # Linear fit line (data-based)
        geom_smooth(method = "lm", se = FALSE, color = "red") +
        labs(
            x = outcome,
            y = "Predicted",
            title = paste("Predicted vs", outcome)
        ) +
        theme_minimal()
        list("p"=p, "metrics"=metrics)

    } else {
        stop(paste0("Invalid family:", family))
    }
}