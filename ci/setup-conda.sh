#!/bin/bash
set -e
set -u

export PS1=

# initialize variables
CONDA_PATH=${CONDA_PATH-./miniconda}
GGSASHIMI_ENV=${GGSASHIMI_ENV-"ggsashimi_${R_VER}_${GGPLOT_VER}"}

# get miniconda install script
if [ ! -e miniconda.sh ]; then
    if [[ "${TRAVIS_PYTHON_VERSION-3}" == "2.7" ]]; then
        COMDA_URL="https://repo.continuum.io/miniconda/Miniconda2-latest-Linux-x86_64.sh"
    else
        COMDA_URL="https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh";
    fi
    wget $COMDA_URL -O miniconda.sh
fi
# install miniconda
if [ ! -d $CONDA_PATH ]; then
    bash miniconda.sh -b -p $CONDA_PATH
fi

# source miniconda profile
source "$CONDA_PATH/etc/profile.d/conda.sh"

# configure miniconda
conda config --set always_yes yes --set changeps1 no
conda update -q conda
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge

# show miniconda info
conda info -a

# create environment
conda create -q -n $GGSASHIMI_ENV pysam r-base=$R_VER r-ggplot2=$GGPLOT_VER r-gridextra r-data.table r-svglite
