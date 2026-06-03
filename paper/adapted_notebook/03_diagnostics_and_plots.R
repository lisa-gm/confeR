# Script 3/3: Diagnostics and plots
# Purpose: run reverse-Bayes checks and generate diagnostic visualizations.
# Run after 01_setup_and_data_prep.R and 02_model_fits_and_comparison.R.

# ---- Reverse-Bayes ----

stopifnot(family=="gaussian")  # not yet implemented

use_local_intercepts <- heterogeneity != Heterogeneity$NONE
box_results <- pred_check(params_oneshot, sumstats, family,
                                   use_local_intercepts, center_name, remove_intercept = use_local_intercepts)
pboxes <- box_results$pboxes
write.csv(pboxes, file.path(outpath, paste0("pboxes", dataname, ".csv")))

cid <- 9
-log2(pboxes)
-log2(pboxes)[cid]
local_glms$sigma2_list
as.vector(t(lapply(data_split, nrow)))
mean(local_glms$sigma2_list)

ggplot(data=NULL, aes(y = -log2(pboxes), x = local_glms$sigma2_list)) +
  geom_point(color = "steelblue") +
  theme_minimal()

# ---- Box prior predictive check ----

gghistogram(as.numeric(pboxes), fill = "lightgray", bins=10)

df <- data.frame(Center = seq_along(pboxes), neg_log_p = -log10(pboxes))

ggplot(df, aes(x = Center, y = neg_log_p)) +
  geom_point(color = "steelblue") +
  theme_minimal() +
  labs(
    x = "Center",
    y = expression(-log[10](p)),
    title = "Manhattan plot for pbox"
  ) +
  geom_hline(yintercept = -log10(0.05 / n_sites), linetype = "dashed", color = "red")

# ---- Forest plot ----

if (heterogeneity == Heterogeneity$RANDOM) stop()
df_forest <- prepare_forest_plot(df_bca, sumstats, alpha=0.05, family=family)
df_forest <- df_forest |> filter(df_forest$Covariate != "sigma2")
# df_forest <- df_forest |> mutate(Covariate = ifelse(Covariate == "sigma2", "sigma^2", Covariate)) |> mutate(
#                             Covariate = factor(
#                             Covariate,
#                             levels = c(setdiff(unique(Covariate), "sigma^2"), "sigma^2")
#                         ))
forest_plot(df_forest, pboxes, outfile = paste0("figures_tmp/forest.", dataname, ".pdf"),
            inline_plot = TRUE, use_log_scale = TRUE, order_box = TRUE)

# sigma2: Forest bca_oneshot with 1 site equals local glm
df_forest <- prepare_forest_plot(df_bca, sumstats, alpha=0.05, family=family)
aa <- unname(unlist(as.vector(df_forest |> filter(Covariate == "sigma2" & site != "Federated") |> select(Estimate))))
bb <- local_glms$sigma2_list
all.equal(aa, bb)

# ---- SessionInfo ----

sessionInfo()

