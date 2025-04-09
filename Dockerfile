FROM rocker/verse:4.4.2

RUN apt-get update && apt-get install -y --no-install-recommends \
    python3 \
    python3-pip \
    wget \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

RUN wget -O /tmp/google-chrome.deb https://dl.google.com/linux/direct/google-chrome-stable_current_amd64.deb \
    && apt-get update \
    && apt-get install -y --no-install-recommends /tmp/google-chrome.deb \
    && rm -f /tmp/google-chrome.deb \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /usr/src/app

RUN Rscript -e "install.packages('devtools', repos='https://cran.rstudio.com/')"

RUN Rscript -e "install.packages(c( \
    'aws.s3', \
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
    'scales' \
  ), repos='https://cran.rstudio.com/')" \
  && Rscript -e "devtools::install_github('cloudyr/aws.ec2metadata')"

COPY app/ app/
COPY entrypoint.sh /usr/src/app/entrypoint.sh

RUN chmod +x /usr/src/app/entrypoint.sh

ENTRYPOINT ["/usr/src/app/entrypoint.sh"]