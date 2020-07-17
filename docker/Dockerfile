ARG R_VER=3.3.3
FROM rocker/r-ver:${R_VER} as builder

# install needed tools
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        libcurl4-openssl-dev \
        libssl-dev \
        zlib1g-dev \
        libncurses5-dev  \
        libncursesw5-dev \
        libcairo2-dev \
        libcurl4-openssl-dev \
        libssl-dev \
        curl \
        bzip2

# install SAMtools
RUN curl -fksSL https://github.com/samtools/samtools/releases/download/1.3.1/samtools-1.3.1.tar.bz2 | tar xj && \
    cd samtools-1.3.1 && \
    make all all-htslib && make install install-htslib

## Install R packages
ARG GGPLOT_VER=2.2.1
RUN R -e 'install.packages("devtools");' && \
    R -e "require(devtools); install_version('ggplot2', version = '${GGPLOT_VER}', repos = 'http://cloud.r-project.org');"

RUN R -e 'install.packages(c("gridExtra", "data.table", "svglite"), repos="http://cloud.r-project.org/");'

FROM rocker/r-ver:${R_VER}

LABEL maintainer "Emilio Palumbo <emilio.palumbo@crg.eu>" \
      version "1.0" \
      description "Docker image for ggsashimi"

# install needed tools
RUN apt-get update && apt-get install -y --no-install-recommends \
    locales \
    libncurses5 \
    bzip2 \
    python

COPY --from=builder /samtools-1.3.1/samtools /usr/local/bin
COPY --from=builder /usr/local/lib/R/site-library /usr/local/lib/R/site-library

# Copy ggsashimi in the docker image
ADD sashimi-plot.py /

# Run the container as an executable
ENTRYPOINT ["/sashimi-plot.py"]

