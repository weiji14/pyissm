# PyISSM

Experimental code on running the [Ice-sheet and Sea-level System Model (ISSM)](https://issm.jpl.nasa.gov) via Python instead of the default Matlab interface.
Focusing on Pine Island Glacier (based on this [tutorial](https://issm.jpl.nasa.gov/documentation/tutorials/pig/)).

## Quickstart

Launch in [Pangeo Binder](https://pangeo-binder.readthedocs.io) (Interactive jupyter lab environment in the cloud).

[![Binder](https://binder.pangeo.io/badge_logo.svg)](https://binder.pangeo.io/v2/gh/weiji14/pyissm/master)

## Installation

Start by cloning this [repo-url](/../../)

    git clone <repo-url>

Then I recommend [using conda](https://conda.io/projects/conda/en/latest/user-guide/install/index.html) to install the dependencies.

    cd pyissm
    conda env create -f environment.yml

Activate the conda environment first.

    conda activate pyissm

Then clone the ISSM svn trunk repository.
You will need to have [git-svn](https://git-scm.com/docs/git-svn) installed.

    echo 'anon' | git svn clone --username anon -r 25855 https://issm.ess.uci.edu/svn/issm/issm/trunk

After that, you will need to compile some of the dependencies shipped with ISSM.

    export ISSM_DIR=/path/to/pyissm/trunk
    source $ISSM_DIR/etc/environment.sh && cd $ISSM_DIR/externalpackages/m1qn3 && ./install.sh
    source $ISSM_DIR/etc/environment.sh && cd $ISSM_DIR/externalpackages/triangle && ./install-linux.sh

Setup some environment variables and configuration settings,
and patch the configure file to remove the 'm' for minimal (Python 3.8 and above only).

    export CONDA_DIR=/path/to/.conda/envs/pyissm
    cd $ISSM_DIR && source $ISSM_DIR/etc/environment.sh && autoreconf -ivf
    cd $ISSM_DIR && \
        ./configure \
        --prefix="$ISSM_DIR" \
        --disable-static \
        --enable-debugging \
        --enable-development \
        --with-numthreads=8 \
        --with-python-version=3.7 \
        --with-python-dir="$CONDA_DIR" \
        --with-python-numpy-dir="$CONDA_DIR/lib/python3.7/site-packages/numpy/core/include/numpy" \
        --with-fortran-lib="-L$CONDA_DIR/lib/gcc/x86_64-conda-linux-gnu/7.5.0/ -lgfortran" \
        --with-mpi-include="$CONDA_DIR/lib/include" \
        --with-mpi-libflags="-L$CONDA_DIR/lib -lmpi -lmpicxx -lmpifort" \
        --with-metis-dir="$CONDA_DIR/lib" \
        --with-scalapack-dir="$CONDA_DIR/lib" \
        --with-mumps-dir="$CONDA_DIR/lib" \
        --with-petsc-dir="$CONDA_DIR" \
        --with-triangle-dir="$ISSM_DIR/externalpackages/triangle/install" \
        --with-m1qn3-dir="$ISSM_DIR/externalpackages/m1qn3/install"
    # sed -i 's/-lpython${PYTHON_VERSION}m/-lpython${PYTHON_VERSION}/g' $ISSM_DIR/configure


For completeness, here's the equivalent commands for compiling under MATLAB.
The 'with-matlab-dir' filepath needs to be set to where MATLAB is installed.

    cd $ISSM_DIR && \
        ./configure \
        --prefix="$ISSM_DIR" \
        --disable-static \
        --enable-debugging \
        --enable-development \
        --with-numthreads=8 \
        --with-matlab-dir="/opt/Matlab/R2019a/" \
        --with-fortran-lib="-L$CONDA_DIR/lib/gcc/x86_64-conda-linux-gnu/7.5.0/ -lgfortran" \
        --with-mpi-include="$CONDA_DIR/lib/include" \
        --with-mpi-libflags="-L$CONDA_DIR/lib -lmpi -lmpicxx -lmpifort" \
        --with-metis-dir="$CONDA_DIR/lib" \
        --with-scalapack-dir="$CONDA_DIR/lib" \
        --with-mumps-dir="$CONDA_DIR/lib" \
        --with-petsc-dir="$CONDA_DIR" \
        --with-triangle-dir="$ISSM_DIR/externalpackages/triangle/install" \
        --with-m1qn3-dir="$ISSM_DIR/externalpackages/m1qn3/install"

Finally compile ISSM!

    cd $ISSM_DIR
    make --jobs=8
    make install

You'll need to export some more environment variables for things to work.
Plus add a line into your .bashrc file to apply the settings each time you restart.
See also official instructions at https://issm.jpl.nasa.gov/documentation/addpath/

    export PYTHONPATH=$ISSM_DIR/bin:$PYTHONPATH
    export PYTHONPATH=$ISSM_DIR/lib:$PYTHONPATH
    echo "source $ISSM_DIR/etc/environment.sh" >> $HOME/.bashrc

### Syncing/Updating to new dependencies

    conda env update -f environment.yml
    cd $ISSM_DIR && svn update
    make clean
    # now recompile ISSM following instructions above

## Running jupyter lab

    conda activate pyissm
    python -m ipykernel install --user --name pyissm  # to install conda env properly
    jupyter kernelspec list --json                    # see if kernel is installed
    jupyter lab &


## TODO

- [x] Use newer dataset inputs (e.g. BedMachine, ALBMAPv1, MEaSUREs Phase Map of Antarctic Ice Velocity)
- [x] Go from Python 2.7 to Python 3
- [x] Get reproducible [binder](https://mybinder.readthedocs.io) build to work (using Dockerfile)
- [ ] Jupyter notebooks, [PyGMT](https://pygmt.org) plots, and so on!

## Notes

- ISSM installed in $ISSM_DIR which is /home/user/pyissm/trunk, activated by running `source $ISSM_DIR/etc/environment.sh`
- Pine Island Tutorial located in $ISSM_DIR/examples/Pig, and can be copied to /home/user/pyissm/Pig
- Datasets found under /home/user/pyissm/Data

## References

- Larour, E., Seroussi, H., Morlighem, M., & Rignot, E. (2012). Continental scale, high order, high spatial resolution, ice sheet modeling using the Ice Sheet System Model (ISSM). Journal of Geophysical Research: Earth Surface, 117(F1). https://doi.org/10.1029/2011JF002140
