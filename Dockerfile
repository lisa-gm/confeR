FROM rocker/verse:4.3.3

## set up directories
WORKDIR /home/rstudio
RUN mkdir /home/rstudio/paper
COPY package /home/rstudio/package

## install R packages from CRAN the last day of the specified R version
RUN install2.r --error --skipinstalled --repos 'https://cloud.r-project.org' \
    remotes BFI ggpubr dplyr tidyr stringr tibble metafor Matrix lme4 invgamma && \
    R -e "remotes::install_url('https://cran.r-project.org/src/contrib/Archive/pda/pda_1.2.8.tar.gz', dependencies = TRUE)" && \
    R CMD INSTALL /home/rstudio/package/out/confeR_0.1.tar.gz

COPY requirements.txt . 

RUN apt-get update && apt-get install -y python3-pip && \
    pip3 install --no-cache-dir -r requirements.txt && \
    apt-get clean && rm -rf /var/lib/apt/lists/*

EXPOSE 8787

ENV USER=rstudio
ENV PASSWORD=rstudio

# Start RStudio Server
CMD ["/init"]