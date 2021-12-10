ARG UBUNTU_VERSION=20.04
FROM ubuntu:${UBUNTU_VERSION} as builder

# install needed tools
RUN apt-get update \
    &&  apt-get install -y --no-install-recommends \
         ca-certificates \
         curl \
         gdebi-core \
         libcairo2-dev \
         libfontconfig1-dev \
         lsb-release \
         python3-pip \
         python-is-python3

# install R
ARG R_VER=3.6.3
RUN export UBUNTU_VER=$(lsb_release -rs | tr -d '.') \
    && curl -O https://cdn.rstudio.com/r/ubuntu-${UBUNTU_VER}/pkgs/r-${R_VER}_1_amd64.deb \
    && DEBIAN_FRONTEND=noninteractve gdebi --non-interactive r-${R_VER}_1_amd64.deb

ENV PATH=/opt/R/${R_VER}/bin:$PATH

RUN echo 'options(repos = "https://cloud.r-project.org/")' > ~/.Rprofile

## Install R packages
ARG GGPLOT_VER=3.3.3
RUN R -e 'install.packages("remotes");' && \
    R -e 'remotes::install_version("ggplot2", version="'${GGPLOT_VER}'")' && \
    R -e 'remotes::install_cran(c("gridExtra", "data.table", "svglite"))'

# Install pysam
RUN pip3 install pysam

FROM ubuntu:${UBUNTU_VERSION}

LABEL maintainer "Emilio Palumbo <emilio.palumbo@crg.eu>" \
      version "1.0" \
      description "Docker image for ggsashimi"

# install needed tools
RUN apt-get update \
    && apt-get install -y --no-install-recommends \
        libblas3 \
        libgomp1 \
        libicu66 \
        liblapack3 \
        locales \
        python-is-python3

# set locale
RUN echo "en_US.UTF-8 UTF-8" >> /etc/locale.gen \
    && locale-gen en_US.utf8 \
    && /usr/sbin/update-locale LANG=en_US.UTF-8 LC_ALL=en_US.UTF-8

# set environment variables
ARG R_VER=3.6.3
ENV LANG=en_US.UTF-8 \
    LC_ALL=en_US.UTF-8 \
    PATH=/opt/R/${R_VER}/bin:$PATH

# Copy installed libraries from builder
COPY --from=builder /opt/R /opt/R
COPY --from=builder /usr/local/lib/python3.8/dist-packages /usr/local/lib/python3.8/dist-packages

# Copy ggsashimi in the docker image
ADD ggsashimi.py /

# Run the container as an executable
ENTRYPOINT ["/ggsashimi.py"]

