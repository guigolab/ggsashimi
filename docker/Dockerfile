FROM debian:stable

LABEL maintainer "Emilio Palumbo <emilio.palumbo@crg.eu>" \
      version "1.0" \
      description "Docker image for ggsashimi"

# install needed tools
RUN apt-get update --fix-missing -qq && \
        apt-get install -y -q \
    curl \
    locales \
    libncurses5-dev  \
    libncursesw5-dev \
    build-essential \
    pkg-config \
    zlib1g-dev \
    bzip2 \
    r-base \
    python \
    libcairo2-dev \
    && apt-get clean \
    && apt-get purge \
    && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

# install SAMtools
RUN curl -fksSL https://github.com/samtools/samtools/releases/download/1.3.1/samtools-1.3.1.tar.bz2 | tar xj && \
    cd samtools-1.3.1 && \
    make all all-htslib && make install install-htslib

## Install R packages for ggplot2
RUN R -e 'install.packages( c("ggplot2", "gridExtra", "data.table", "svglite"), repos="http://cloud.r-project.org/");' 

# Copy ggsashimi in the docker image
ADD sashimi-plot.py /

# Run the container as an executable
ENTRYPOINT ["/sashimi-plot.py"]

