# confeR

The **confeR** R package provides methods for performing conjugate federated inference in R, as described in Degen, Pawel, & Held (2026) [DOI].

## Installation

```r
## CRAN version
install.packages("XXX")

## development version from GitLab
install.packages("remotes") # requires remotes package
remotes::install_gitlab(repo = "crsuzh/confeR", subdir = "package", host = "gitlab.uzh.ch")
```

## Usage

See [example.ipynb](example.ipynb) for a more thorough explanation.

```r
library(confeR)

data(nurses_hom, package = "confeR")
summary_stats <- nurses_hom$summary_stats

# If true, use fixed local intercepts for each site, else use only a global intercept
use_local_intercepts <- TRUE 

params_oneshot <- bca_oneshot(summary_stats use_local_intercepts, family="gaussian")
df_bca <- tidy_results(params_oneshot, use_local_intercepts)
df_bca
```
