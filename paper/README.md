# Reproducing the results

## Set up the environment

Make sure to have Docker and Make installed, then run `make docker` from the root directory of this git repository. This will install all necessary dependencies. RStudio Server can then be opened from a browser (<http://localhost:8787>).

## Computing summary results and reproducing Figure 6

In RStudio Server, open `paper\notebook.qmd` and click "Render" (or Ctrl + Shift + K). The forest plot created at the end of the rendered document is Figure 6.

## Reproducing the remaining figures

To reprdouce Figures 2–5, you first need to create summary data using `paper\notebook.qmd` and then create the respective figure using `paper\figures.ipynb`.

### Figure 2

1. Set `dataname <- "trauma_shuffled"` and `heterogeneity <- Heterogeneity$NONE` in `notebook.qmd` and run it to create `data\summarized\params.trauma_shuffled.csv`.
2. Set `dataname = "trauma_shuffled"` in `figures.ipynb` and run the notebook.

### Figure 3 and Supplementary Figure 2

1. Set `dataname <- "Nurses"` and `heterogeneity <- Heterogeneity$NONE` in `notebook.qmd` and run it to create `data\summarized\params.Nurses.csv` and Supp. Fig. 2.
2. Set `dataname = "Nurses"` in `figures.ipynb` and run the notebook.

### Figure 4 and Supplementary Figure 3

1. Set `dataname <- "Nurses"` and `heterogeneity <- Heterogeneity$FIXED` in `notebook.qmd` and run it to create `data\summarized\params.Nurses_local_int.csv` and Supp. Fig. 3.
2. Set `dataname = "Nurses_local_int"` in `figures.ipynb` and run the notebook.

### Figure 5

Assuming you have created the summary data from the previous figures, run the cells in the section `Manhattan (Fig. 5).` in `figures.ipynb`.
