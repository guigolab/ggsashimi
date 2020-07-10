ARG R_VER=3.3.3
FROM rocker/r-ver:${R_VER}

RUN apt-get update && \
    apt-get install -y --no-install-recommends libcurl4-openssl-dev libssl-dev

ARG GGPLOT_VER=2.2.1
RUN R -e 'install.packages("devtools");' && \
    R -e "require(devtools); install_version('ggplot2', version = '${GGPLOT_VER}', repos = 'http://cloud.r-project.org');"
    
RUN R -e 'install.packages(c("gridExtra", "data.table", "svglite"), repos="http://cloud.r-project.org/");'

