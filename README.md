# Bayesian conjugate analysis for federated statistical inference

This repository contains

1. `./package` The R package **confeR** to perform Bayesian conjugate analysis for federated inference

2. `./paper` Code and data to reproduce result from the paper: *Degen, P. M., Pawel, S., Held.
   L. (2025). Bayesian conjugate analysis for federated statistical inference. [DOI]*

To cite our work, use the following BibTeX reference

```BibTeX
@article{Degen2025,
  year = {2025},
  author = {Peter Methys Degen and Samuel Pawel and Leonhard Held},
  title = {Bayesian conjugate analysis for federated statistical inference},
  journal = {},
  doi = {}
}
```

## Reproducing the paper with Docker (TO DO)

Make sure to have Docker and Make installed, then run `make docker` from the
root directory of this git repository. This will install all necessary
dependencies. RStudio Server can then be opened from a browser
(<http://localhost:8787>), and the Quarto notebooks in `./paper`, which
contains all code for the results from the paper, can be rerun (make sure to set
the working directory to `./paper` when running R interactively).

## To do

- Dockerize
- Improve documentation
- More tests
- pBox check for normal model
- Random effects for normal model
