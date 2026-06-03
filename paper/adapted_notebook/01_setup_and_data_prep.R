##### slightly adapted version of paper/notebook.ipynb #######
# Script 1/3: Setup and data preparation
# Purpose: load packages, configure dataset/options, and build split site data.
# Run first.

# ---- Load and pre-process data ----
library(this.path)
# current directory
script_dir <- this.path::this.dir()
# set working directory to parent of script dir
setwd(paste0(script_dir, "/.."))

# Skipped when running in Docker container
if (file.exists(".Rprofile"))
    source(".Rprofile", local = TRUE)

library(dplyr)
library(ggplot2)
library(tidyr)
library(confeR)

source("scripts/functions.R")

set.seed(123)

# set output file path for data, relative to script dir
outpath <- file.path(script_dir, "data/summarized")

# Create if it doesn't exist
if (!dir.exists(outpath)) {
  dir.create(outpath, recursive = TRUE)
}

# To reproduce Fig. 6 in the paper, select dataname `dataname <- "nurses_hom"` and `heterogeneity <- Heterogeneity$NONE`.
# To compute summary statistics for Fig. 2, select `dataname <- "trauma_shuffled"` and `heterogeneity <- Heterogeneity$NONE`. The notebook will throw an error in the Reverse-Bayes section, which is not an issue.
# To compute summary statistics for Fig. 3, select `dataname <- "Nurses"` and `heterogeneity <- Heterogeneity$NONE`.
# To compute summary statistics for Fig. 4, select `dataname <- "Nurses"` and `heterogeneity <- Heterogeneity$FIXED`.
# Fig. 5 uses summary statistics computed in the above steps.

###### Select data set ######

## Logistic regression

dataname <- "trauma"
#dataname <- "trauma_shuffled"

## Linear regression

#dataname <- "Nurses"
#dataname <- "nurses_hom"

#############################

# Name of the column designating the center id
center_name <- "hospital"

# Enum, don't change vals
Heterogeneity <- list(
    NONE = "none",     # global intercept
    FIXED = "fixed",   # local fixed intercepts
    RANDOM = "random"  # local random intercepts
)

heterogeneity <- Heterogeneity$NONE
#heterogeneity <- Heterogeneity$FIXED
#heterogeneity <- Heterogeneity$RANDOM

# Set TRUE to pool all rows into one synthetic site.
use_single_site <- FALSE

use_local_variances <- FALSE

if (startsWith(dataname, "trauma")) {
        data(trauma, package = "BFI")
        data_raw <- trauma
} else if (startsWith(dataname, "nurses_hom")) {
        data(nurses_hom, package = "confeR")
        data_raw <- nurses_hom$data_ipd
} else if (startsWith(dataname, "Nurses")) {
        data(Nurses, package = "BFI")
        data_raw <- Nurses
} else if (startsWith(dataname, "LOS")) {
        data(LOS, package = "pda")
        data_raw <- LOS
    data_raw <- dplyr::rename(data_raw, !!center_name := site)
        data_raw$hospital <- as.factor(sub("^site", "", data_raw$hospital))
} else {
        data_raw <- read.csv(paste0("data/raw/", dataname, ".csv"), row.names = 1)
}

# Shuffling center simulates homogeneous case
if (endsWith(dataname, "_shuffled")) {
        data_raw[[center_name]] <- sample(data_raw[[center_name]])
}

if (endsWith(dataname, "_sub")) {
    print("Subsampling data...")
    subsample_prop <- 0.2 # keep this fraction of patients
    data_raw <- data_raw |>
    dplyr::group_by(hospital) |>
    dplyr::slice_sample(prop = subsample_prop, replace = FALSE) |>
    dplyr::ungroup()
}

if (use_single_site) {
    data_raw[[center_name]] <- "site_all"
    if (heterogeneity != Heterogeneity$NONE) {
        message("use_single_site=TRUE: forcing heterogeneity <- Heterogeneity$NONE")
        heterogeneity <- Heterogeneity$NONE
    }
}

raw_dataname <- dataname
if (heterogeneity == Heterogeneity$FIXED) {
    dataname <- paste0(dataname, "_", "local_int")
} else if (heterogeneity == Heterogeneity$RANDOM) {
    dataname <- paste0(dataname, "_", "local_int_re")
}
if (use_local_variances) {
    dataname <- paste0(dataname, "_", "local_var")
}

data_raw[[center_name]] <- factor(data_raw[[center_name]],
                        levels = sort(unique(data_raw[[center_name]]),
                        decreasing = FALSE))

# Convert factors to int, only valid for two-level factors (0, 1)
data_raw[] <- lapply(data_raw, function(col) {
  if (is.factor(col)) {
    lev <- levels(col)
    if (setequal(lev, c("0", "1"))) {
      return(as.numeric(as.character(col)))
    } else {
      print("Warning: factors with more than two levels not yet supported; use dummy-encoding")
    }
  }
  col
})

dataname
head(data_raw)

table(data_raw[[center_name]])
median(table(data_raw[[center_name]]))
sort(unique(data_raw[[center_name]]), decreasing = FALSE)

scale_numeric_cols <- function(data) {
  data[] <- lapply(data, function(x) if (is.numeric(x)) scale(x)[, 1] else x)
  data
}

scale_all_cols <- function(data, center = TRUE, scale = TRUE) {
    require(dplyr)
    data |>
    dplyr::mutate(dplyr::across(where(is.numeric), ~ as.numeric(scale(.x, center = center, scale = scale))))
}

one_hot <- function(df, cols_to_encode) {
    dummies <- model.matrix(~ . - 1, data = df[cols_to_encode])
    cbind(df[setdiff(names(df), cols_to_encode)], dummies)
}

if (startsWith(dataname, "trauma")) {
    data_norm <- data_raw |>
        # don't need to normalize within group as combined mean and sd can be computed without combining data
        dplyr::mutate(
            age = scale(age),
            ISS = scale(ISS),
            GCS = scale(GCS)
        )

    model <- mortality ~ sex + age + ISS + GCS
    model_no_int <- as.formula(paste("mortality ~ 0 + sex + age + ISS + GCS +", center_name))
    family <- "binomial" #(link = "logit")

} else if (startsWith(dataname, "Nurses") || startsWith(dataname, "sim-nurses") || startsWith(dataname, "nurses_hom")) {
    if ("experience" %in% colnames(data_raw))
        data_raw <- data_raw |> dplyr::rename(experien=experience)
    data_norm <- data_raw |>
    dplyr::mutate(
            age = scale(age),
            experience = scale(experien),
            stress = scale(stress),

            # keep numeric for now; fine if binary 0, 1
            #gender = factor(gender),
            #wardtype = factor(wardtype)
    )
    #data_norm <- scale_numeric_cols(data_raw) |> mutate(experience=experien)
    model_no_int <- as.formula(paste("stress ~ 0 + gender + age + experience +", center_name))
    if (heterogeneity != Heterogeneity$NONE) {
        model <- stress ~ gender + age + experience  # BFI paper leaves out wardtype
    } else {
        model <- stress ~ gender + age + experience + wardtype
    }
    family <- "gaussian"
} else if (startsWith(dataname, "LOS")) {
    data_norm <- data_raw |>
    dplyr::mutate(
            lab = scale(lab),
            los = scale(los),
    )
    data_norm <- one_hot(data_norm, c("age", "sex")) |> dplyr::rename("sex"="sexM")
    model_no_int <- as.formula(paste("los ~ 0 + sex + ageyoung + ageold + lab +", center_name))
    model <- los ~ sex + ageyoung + ageold + lab
    family <- "gaussian"

} else if (startsWith(dataname, "sim-linear")) {
    data_norm <- scale_numeric_cols(data_raw)
    predictors <- setdiff(names(data_norm), c("y", "hospital"))
    model <- as.formula(paste("y ~", paste(predictors, collapse = " + ")))
    if (heterogeneity != Heterogeneity$NONE) {
        model_no_int <- as.formula(paste("y ~ 0 +", paste(predictors, collapse = " + "), "+ hospital"))
    }
    family <- "gaussian"
} else {
  stop()
}

covariates <- attr(terms(model), "term.labels")
#covariates_local <- covariates[covariates != center_name]
outcome <- all.vars(model)[1]

if (center_name %in% covariates) {
    data <- data_norm[c(outcome, covariates)]
} else {
    data <- data_norm[c(center_name, outcome, covariates)]
}

if (endsWith(dataname, "_glm")) {
    family_alt <- "gaussian_forced_glm"
} else {
    family_alt <- family
}

sites <- levels(data[[center_name]])

### 
data_split <- data |>
    dplyr::group_by(.data[[center_name]]) |>
    dplyr::group_split()

n_sites <- length(data_split)
head(data_split[[1]])
n_data <- sum(unlist(lapply(data_split, nrow)))

if (heterogeneity != Heterogeneity$NONE) {
    p <- length(covariates) + n_sites
} else {
    p <- length(covariates)
}
p


#### store to file to be able to rerun with DALIA in python 
file_name <- file.path(outpath, paste0("data_", dataname, "_", family, ".csv"))
write.csv(data, file_name, row.names = FALSE)

