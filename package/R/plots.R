# Reduce local intercepts extended covariate matrix
bca_remove_local_intercepts <- function(sumstats, center_identity, family, alpha = 0.05) {
    bs <- sumstats[[center_identity]]
    xx <- bs$xx
    xy <- bs$xy
    yy <- bs$yy
    n <- bs$n

    keep_intercept <- paste0("^Intercept_", center_identity, "$")
    drop_pattern   <- "^Intercept_.*$"
    drop_rows_cols <- grep(drop_pattern, colnames(xx), value = TRUE)
    drop_rows_cols <- setdiff(drop_rows_cols, grep(keep_intercept,
                                                drop_rows_cols,
                                                value = TRUE))

    if (length(drop_rows_cols) > 0) {
        xx <- xx[!rownames(xx) %in% drop_rows_cols, !colnames(xx) %in% drop_rows_cols, drop = FALSE]
        rnxy <- rownames(xy)[!rownames(xy) %in% drop_rows_cols]
        xy <- as.matrix(xy[rownames(xy) %in% rnxy])
        rownames(xy) <- rnxy
    }
    bs_clean <- list(list("xx"=xx, "xy"=xy, "yy"=yy, "n"=n))
    confeR::bca_oneshot(bs_clean, family, use_local_intercepts=FALSE, alpha=alpha)
}


# Build a diamond data frame for forest plot
make_diamond <- function(df, half_height = 0.15) {
  # df must contain: Estimate, lower, upper, site, Covariate
  # We create four vertices for each row:
  # (lower , y) → (Estimate , y‑offset) → (upper , y) → (Estimate , y+offset)

  df$Estimate_again <- df$Estimate
  x <- df |> dplyr::select(c("lower", "Estimate", "upper", "Estimate_again"))
  x <- c(t(x))

  data.frame(
    Covariate = rep(df$Covariate, each = 4),
    site      = rep(df$site,      each = 4),
    x = x,
    y         = c(0,
                  0 - half_height,
                  0,
                  0 + half_height)
  )
}

#' @title Prepare forest plot
#'
#' @description Get clean dataframe with local glm parameters from summary statistics for all
#'   sites plus the federated BCA estimates.
#'
#' @param df_bca Dataframe with BCA estimates returned by `tidy_results`.
#' @param sumstats List of summary statistics from each local site, e.g. obbject
#'   returned by `bca_iterate_sites`.
#' @param alpha Numeric. Significance level for credible intervals. Default
#'   0.05.
#'
#' @return Dataframe with parameter estimates and confidence intervals for all
#'   sites and covariates.
#'
#' @author Peter Degen
#'
#' @export
prepare_forest_plot <- function(df_bca, sumstats, family, alpha=0.05) {
    out_list <- vector("list", length(sumstats))

    for (i in seq_along(sumstats)) {
        bl <- bca_remove_local_intercepts(sumstats, i, family, alpha = alpha)
        bl <- confeR::tidy_results(bl, use_local_intercepts=FALSE)
        bl$Method <- NULL
        bl$site <- i
        out_list[[i]] <- bl
    }

    df_bca$site <- "Federated"
    out_list[[i+1]] <- dplyr::select(df_bca, !"Method")
    df_forest <- do.call(rbind, out_list)
    df_forest <- df_forest[!startsWith(df_forest$Covariate, "Intercept_"), ]
    df_forest
}

#' @title Multivariate forest plot
#'
#' @description Creates a forest plot with a facet for each covariate.
#'   Optionally prints Box's prior predictive tail probability for each site.
#'
#' @param df_forest Dataframe returned by `prepare_forest_plot`.
#' @param pboxes Optional numeric vector of Box p-values, contained in list
#'   returned by `box_check_all_sites`.
#' @param outfile Optional character. If not null, save figure in this file.
#' @param alpha Numeric. Significance level for credible intervals. Default
#'   0.05.
#' @param nrow Numeric. (default: 1)
#' @param inline_plot Logical. If TRUE, adjust width and height of inline plot
#'   in Jupyter notebooks using `options()`. (default: FALSE)
#' @param use_log_scale Logical. If TRUE, print pbox as an S-value, i.e. -log2(pbox) (default: FALSE)
#' @param order_box FALSE Logical. If TRUE, order sites by p-Box. (default: FALSE)
#'
#' @return ggplot2 object.
#'
#' @author Peter Degen
#'
#' @export
forest_plot <- function(df_forest,
                        pboxes=NULL,
                        outfile=NULL,
                        alpha=0.05,
                        nrow=1,
                        inline_plot=FALSE,
                        use_log_scale=TRUE,
                        order_box=TRUE
                        ) {

    require(ggplot2)
    require(ggsci)

    right_covariate <- sort(unique(df_forest$Covariate), decreasing = TRUE)[1]
    pbox_thresh <- alpha

    if (!is.null(pboxes)) {

        pbox_thresh <- alpha / length(pboxes)

        if (use_log_scale) {
            pbox_thresh <- -log2(pbox_thresh)
            pboxes <- -log2(pboxes)
        }

        margins <- margin(t = 5, r = 70, b = 5, l = 5, unit = "pt")
        p_box_df <- data.frame(
            site   = unique(dplyr::filter(df_forest, df_forest$site != "Federated")$site),
            p_Box  = pboxes,
            stringsAsFactors = FALSE
        )
        # Build a data frame that contains the p‑values *only* for that facet
        p_box_df_right <- p_box_df |>
            dplyr::mutate(Covariate = right_covariate) |>
            dplyr::select(.data$Covariate, .data$site, .data$p_Box)

        p_title_df <- data.frame(
        Covariate = right_covariate,
        site      = max(as.integer(factor(df_forest$site,
                                        levels = rev(unique(df_forest$site))))),
        label     = ifelse(use_log_scale, "~~~~s[Box]", "~~~~p[Box]") # will be parsed as p₍Box₎
        )
    } else {
        margins <- margin(t = 5, r = 5, b = 5, l = 5, unit = "pt")
    }

    color_site <- ggsci::pal_npg("nrc")(3)[3]
    color_diamond <- ggsci::pal_npg("nrc")(2)[1]

    if (order_box) {
        p_box_df_right <- p_box_df_right[order(p_box_df_right$p_Box), ]
        df_forest <- df_forest[order(match(df_forest$site, p_box_df_right$site, nomatch = Inf)), ]
    }

    ncovs <- length(unique(df_forest$Covariate))
    nsites <- length(unique(df_forest$site)) -1
    width <- ncovs * 3
    height <- length(unique(df_forest$site)) / 3

    if (inline_plot)
        options(repr.plot.width = width, repr.plot.height = height, repr.plot.res = 100)

    ## Plot
    p <- ggplot(df_forest,
        aes(x = .data$Estimate,
            ymax = nsites+0.5,
            ymin = -0.5,
            y = factor(.data$site, levels = rev(unique(.data$site))))) +

    geom_vline(data = subset(df_forest, df_forest$site == "Federated"),
        aes(xintercept = .data$Estimate), linetype = "dashed", color = color_diamond) +

    geom_vline(data = subset(df_forest, df_forest$site == "Federated"),
        aes(xintercept = 0), linetype = "solid", color = "grey") +

    # Site estimates
    geom_errorbarh(data = subset(df_forest, df_forest$site != "Federated"),
            aes(xmin = .data$lower, xmax = .data$upper), height = 0.3) +
    geom_point(
        data = subset(df_forest, df_forest$site != "Federated"),
        shape = 16, colour = color_site, size=3
    ) +

    ## Federated estimate
    geom_polygon(
        data = make_diamond(subset(df_forest, df_forest$site == "Federated")),
        aes(x = .data$x, y = .data$y, group = interaction(.data$Covariate, .data$site)),
        fill = color_diamond, colour = color_diamond
    )

    if (!is.null(pboxes)) {
        if (use_log_scale)
            fontface <- ifelse(p_box_df_right$p_Box > pbox_thresh, "bold", "plain")
        else
            fontface <- ifelse(p_box_df_right$p_Box < pbox_thresh, "bold", "plain")
        p <- p +
        ## p‑value labels only on the right‑most facet
        geom_text(data = p_box_df_right,
                #aes(label = format(p_Box, digits = 3),
                aes(
                    label = sprintf("    %.2f", .data$p_Box),
                    y     = factor(.data$site, levels = rev(unique(.data$site))),
                    x     = Inf,
                    fontface = fontface,
                ),
                hjust = 0,
                vjust = 0.5,
                size = 5,
                colour = "black",
                inherit.aes = FALSE) +
        # pBox label
        geom_text(
            data = p_title_df,
            aes(x = Inf,
                y = .data$site,
                label = .data$label),
            hjust   = 0,
            vjust   = -0.7,
            size    = 5,
            colour  = "black",
            parse   = TRUE,  # parse string
            inherit.aes = FALSE
        )
    }

    p <- p + facet_wrap(~ Covariate, scales = "free_x", nrow=nrow, labeller = label_parsed) +

    ## Make room on the right so the text isn't clipped
    coord_cartesian(clip = "off") +                     # allow drawing outside panels
    labs(x = "Effect estimate (95 % CI)", y = "Site") +
    theme_minimal(base_size = 16) +
    theme(panel.grid = element_blank(),
        strip.text   = element_text(face = "bold"),
        axis.text.y  = element_text(size = 15),
        panel.border = element_rect(colour = "black", fill=NA, linewidth=0.5),
        panel.spacing = unit(20, "pt"),
        plot.margin  = margins)

    if (!is.null(outfile)) {
        ggsave(outfile, bg = "white", width = width, height = height)
    }
    p
}

#' @title Multivariate forest plot wrapper
#'
#' @description Wrapper that computes pbox and makes forest plot in one call.
#'
#' @param params_oneshot Parameters returned by `bca_oneshot`.
#' @param sumstats List of summary statistics from each local site, e.g. obbject
#'   returned by `bca_iterate_sites`.
#' @param family Family to pass to glm() call, e.g. "gaussian".
#' @param use_local_intercepts Logical. If true, use fixed site-specific
#'   intercepts for each local site.
#' @param center_name Character. Name of covariate denoting site identity.
#' @param alpha Numeric. Significance level for credible intervals. Default
#'   0.05.
#' @param ... Further arguments passed to `prepare_forest_plot`.
#' @return ggplot2 object.
#'
#' @author Peter Degen
#'
#' @export
forest_box <- function(params_oneshot,
                            sumstats,
                            family,
                            use_local_intercepts,
                            center_name="site",
                            remove_intercept=FALSE,
                            alpha=0.05,
                            remove_sigma2=TRUE,
                            ...
                            ) {

    df_bca <- tidy_results(params_oneshot, use_local_intercepts)
    box_results <- confeR::pred_check(params_oneshot,
                                sumstats,
                                family = family,
                                use_local_intercepts,
                                center_name = center_name,
                                remove_intercept = remove_intercept)

    df_forest <- confeR::prepare_forest_plot(df_bca, sumstats, alpha=0.05, family="gaussian")
    if (remove_sigma2)
        df_forest <- df_forest |> dplyr::filter(df_forest$Covariate != "sigma2")
    forest_plot(df_forest, box_results$pboxes, ...)
}
