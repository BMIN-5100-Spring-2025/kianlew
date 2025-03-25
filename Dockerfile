FROM rocker/verse:4.4.2

# 1) Install system dependencies + tools needed for Google Chrome
RUN apt-get update && apt-get install -y --no-install-recommends \
    python3 \
    python3-pip \
    wget \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# 2) Download and install Google Chrome .deb package
RUN wget -O /tmp/google-chrome.deb https://dl.google.com/linux/direct/google-chrome-stable_current_amd64.deb \
    && apt-get update \
    && apt-get install -y --no-install-recommends /tmp/google-chrome.deb \
    && rm -f /tmp/google-chrome.deb \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /usr/src/app

# 3) Install required R packages from CRAN in one command
RUN Rscript -e "install.packages(c( \
    'webshot2', \
    'pagedown', \
    'readr', \
    'dplyr', \
    'tidyverse', \
    'tidyr', \
    'data.table', \
    'gtsummary', \
    'gt', \
    'ggplot2', \
    'scales', \
    'aws.s3' \
  ), repos='https://cran.rstudio.com/')"

# 4) Copy your R scripts (or entire app) to the container
COPY app/ app/

# 5) Default command: run your R script
CMD ["Rscript", "app/final_project_template_bmin5100.R"]
