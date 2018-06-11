ARG R_VER=3.3.3
FROM rocker/r-ver:${R_VER} as builder

# install needed tools
RUN apt-get update \
    && apt-get install -y --no-install-recommends \
        libbz2-dev \
        libcairo2-dev \
        libcurl4-openssl-dev \
        libfontconfig1-dev \
        liblzma-dev \
        libssh2-1-dev \
        libssl-dev \
        python \
        python-dev \
        python-pip \
        zlib1g-dev

# Install pysam
RUN pip install pysam

## Install R packages
ARG GGPLOT_VER=2.2.1
RUN R -e 'install.packages("devtools");' && \
    R -e "require(devtools); install_version('ggplot2', version = '${GGPLOT_VER}', repos = 'https://cloud.r-project.org');"

RUN R -e 'install.packages(c("gridExtra", "data.table", "svglite"), repos="https://cloud.r-project.org/");'

FROM rocker/r-ver:${R_VER}

LABEL maintainer "Emilio Palumbo <emilio.palumbo@crg.eu>" \
      version "1.0" \
      description "Docker image for ggsashimi"

# install needed tools
RUN apt-get update && apt-get install -y --no-install-recommends \
    bzip2 \
    liblzma5 \
    locales \
    python \
    zlib1g

# Copy installed libraries from builder
COPY --from=builder /usr/local/lib/R/site-library /usr/local/lib/R/site-library
COPY --from=builder /usr/local/lib/python2.7/dist-packages/pysam /usr/local/lib/python2.7/dist-packages/pysam

# Copy ggsashimi in the docker image
ADD sashimi-plot.py /

# Run the container as an executable
ENTRYPOINT ["/sashimi-plot.py"]

