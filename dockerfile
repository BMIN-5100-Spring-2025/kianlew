FROM rocker/verse:4.4.2

RUN apt-get update && apt-get install -y --no-install-recommends \
    python3 \
    python3-pip \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /usr/src/app

RUN Rscript -e "install.packages(c( \
    'readr', \
    'dplyr', \
    'tidyverse', \
    'tidyr', \
    'data.table', \
    'gtsummary', \
    'gt', \
    'ggplot2', \
    'scales' \
  ), repos='http://cran.rstudio.com/')"

COPY app/ app/
COPY data/ data/
COPY tests/ tests/

CMD ["Rscript", "app/final_project_template_bmin5100.R"]




