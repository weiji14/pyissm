FROM buildpack-deps:focal@sha256:c1557d7311f0f6ec8f47d179a5d93644a2b1b7234c1fbf47757e735909ad8b84
LABEL maintainer "https://github.com/weiji14"
ENV LANG C.UTF-8
ENV LC_ALL C.UTF-8

# Initiate docker container with user 'jovyan'
ARG NB_USER=jovyan
ARG NB_UID=1000
ENV NB_USER ${NB_USER}
ENV NB_UID ${NB_UID}
ENV HOME /home/${NB_USER}

RUN adduser --disabled-password \
    --gecos "Default user" \
    --uid ${NB_UID} \
    ${NB_USER}

# Setup conda
ENV CONDA_DIR ${HOME}/.conda
ENV NB_PYTHON_PREFIX ${CONDA_DIR}
ENV MINICONDA_VERSION 4.8.2

RUN cd /tmp && \
    wget --quiet https://repo.continuum.io/miniconda/Miniconda3-py38_${MINICONDA_VERSION}-Linux-x86_64.sh && \
    echo "cbda751e713b5a95f187ae70b509403f *Miniconda3-py38_${MINICONDA_VERSION}-Linux-x86_64.sh" | md5sum -c - && \
    /bin/bash Miniconda3-py38_${MINICONDA_VERSION}-Linux-x86_64.sh -f -b -p $CONDA_DIR && \
    rm Miniconda3-py38_${MINICONDA_VERSION}-Linux-x86_64.sh && \
    $CONDA_DIR/bin/conda config --system --prepend channels conda-forge && \
    $CONDA_DIR/bin/conda config --system --set auto_update_conda false && \
    $CONDA_DIR/bin/conda config --system --set show_channel_urls true && \
    $CONDA_DIR/bin/conda config --system --set pip_interop_enabled true && \
    $CONDA_DIR/bin/conda clean --all --quiet --yes && \
    $CONDA_DIR/bin/conda init --verbose

# Setup $HOME directory with correct permissions
USER root
RUN chown -R ${NB_UID} ${HOME}
USER ${NB_USER}
WORKDIR ${HOME}

USER root
ENV DEBIAN_FRONTEND noninteractive
RUN apt -qq update && apt install -y --no-install-recommends \
    git-svn \
    #python3-dev \
    #python3-matplotlib \
    #python3-minimal \
    #python3-numpy \
    #python3-pip \
    #python3-scipy \
    && rm -rf /var/lib/apt/lists/*
USER $NB_UID

# Use git svn to clone the ISSM svn trunk repository from revision 24686 (02 Apr 2020)
RUN echo 'anon' | git svn clone --username anon -r 24686 https://issm.ess.uci.edu/svn/issm/issm/trunk

# Change to interactive bash shell, so that `conda activate` works
SHELL ["/bin/bash", "-ic"]
RUN conda init --verbose

# Install dependencies in environment.yml file using conda
COPY environment.yml ${HOME}
RUN conda env update -n base -f environment.yml && \
    conda clean --all --yes && \
    conda list -n base

# Install dependencies shipped with ISSM
ENV ISSM_DIR ${HOME}/trunk
#RUN unset F90 && \
#    cd $ISSM_DIR/externalpackages/mpich && \
#    ./install-3.0-linux64.sh
#RUN . $ISSM_DIR/etc/environment.sh && \
#    cd $ISSM_DIR/externalpackages/petsc && \
#    ./install-3.11-linux64.sh
RUN . $ISSM_DIR/etc/environment.sh && \
    cd $ISSM_DIR/externalpackages/m1qn3 && \
    ./install.sh
RUN . $ISSM_DIR/etc/environment.sh && \
    cd $ISSM_DIR/externalpackages/triangle && \
    ./install-linux.sh

# Setup environment variables and configuration settings
RUN cd $ISSM_DIR && \
    . $ISSM_DIR/etc/environment.sh && \
    autoreconf -ivf

RUN cd $ISSM_DIR && \
    ./configure \
    --prefix="$ISSM_DIR" \
    --disable-static \
    --enable-development \
    --with-numthreads=8 \
    --with-python-version=3.7 \
    --with-python-dir="$CONDA_DIR" \
    --with-python-numpy-dir="$CONDA_DIR/lib/python3.7/site-packages/numpy/core/include/numpy" \
    --with-fortran-lib="-L$CONDA_DIR/lib/gcc/x86_64-conda_cos6-linux-gnu/7.3.0 -lgfortran" \
    --with-mpi-include="$CONDA_DIR/lib/include" \
    --with-mpi-libflags="-L$CONDA_DIR/lib -lmpi -lmpicxx -lmpifort" \
    --with-metis-dir="$CONDA_DIR/lib" \
    --with-scalapack-dir="$CONDA_DIR/lib" \
    --with-mumps-dir="$CONDA_DIR/lib" \
    --with-petsc-dir="$CONDA_DIR" \
    --with-triangle-dir="$ISSM_DIR/externalpackages/triangle/install" \
    --with-m1qn3-dir="$ISSM_DIR/externalpackages/m1qn3/install"

# Patch configure file to remove the 'm' for minimal (Python 3.8 and above only)
# RUN sed -i 's/-lpython${PYTHON_VERSION}m/-lpython${PYTHON_VERSION}/g' $ISSM_DIR/configure

# Compile ISSM
RUN cd $ISSM_DIR && \
    make --jobs=8 && \
    make install

# Make python packages available
ENV PYTHONPATH $ISSM_DIR/bin:$PYTHONPATH
ENV PYTHONPATH $ISSM_DIR/lib:$PYTHONPATH

# Ensure that the ISSM environment is loaded
RUN echo "source $ISSM_DIR/etc/environment.sh" >> $HOME/.bashrc

# Copy remaining files to $HOME
COPY --chown=1000:1000 . ${HOME}

EXPOSE 8888
RUN echo -e '#!/bin/bash -i\nset -e\nconda run "$@"' > .entrypoint.sh && \
    chmod +x .entrypoint.sh
ENTRYPOINT ["./.entrypoint.sh"]
CMD ["jupyter", "lab", "--ip", "0.0.0.0"]
