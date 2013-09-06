# Bayesic

Bayesic is an application for analysing pairwise genetic interactions in case/control studies. It provides several methods to choose from and is efficiently implemented in C++.

The following methods are currently implemented:

* Bayesian model posterior (supports covariates through heuristic, is very fast).
* Logistic regression with Wald test on cross-term (requires iteration, and is a bit slower).
* Loglinear model with analytic test (does not support covariates, but is very fast).

The following things will be added at some point:

* Support for imputed data (mach, IMPUTE2).
* Support for continuous phenotypes, like intermediate traits.
* SNP-environment interaction tests.
* Automatic power calculation based on MAF distribution.
* Helper-scripts for running large jobs on a cluster.
* More tests for genetic interactions.

## For the impatient

Make sure you have the Armadillo matrix library installed. Then to build and run (note --recursive flag to get submodules):

    > git clone --recursive https://github.com/fadern/bayesic
    > cd bayesic
    > git submodule init
    > git submodule update
    > mkdir build
    > cd build
    > cmake ../
    > make

    > ls /data/
    dataset.pair
    dataset.bim
    dataset.fam
    dataset.bed

    > cat dataset.pair
    rs412512 rs516161
    rs51512 rs151251
    rs51512 rs516163

    > ./src/bayesic -m bayes /data/dataset.pair /data/dataset
    snp1        snp2        Posterior   N
    rs412512    rs516161    0.9602      3418
    rs51512     rs151251    0.1012      3412
    rs51512     rs516163    0.2312      3416

## Dependencies

Bayesic requires you to have the following packages installed:

* [CMake](http://www.cmake.org/cmake/resources/software.html) (A build system for C and C++)
* [Armadillo](http://arma.sourceforge.net/download.html) (C++ matrix library which in turn depends on BLAS and LAPACK)

If you install them locally (as on a cluster) just specify -DCMAKE_PREFIX_PATH to tell cmake where the libraries are.

## Installing

Before starting it is important to note when building armadillo you need to edit the config.hpp, and add support for
blas and lapack, otherwise you may experience linking errors when building bayesic.

To build with dependencies without root access (note --recursive flag to get submodules):

    > mkdir ~/prefix

    > wget http://sourceforge.net/projects/arma/files/armadillo-3.910.0.tar.gz
    > tar -xf armadillo-3.910.0.tar.gz
    > cd armadillo-3.910.0
    > mkdir build
    > cd build
    > cmake ../ -DCMAKE_INSTALL_PREFIX=$HOME/prefix
    > make install

    > git clone --recursive https://github.com/fadern/bayesic
    > cd bayesic
    > mkdir build
    > cd build
    > cmake ../ -DCMAKE_PREFIX_PATH=$HOME/prefix

To build with dependencies with root access:
    
    > wget http://sourceforge.net/projects/arma/files/armadillo-3.910.0.tar.gz
    > tar -xf armadillo-3.910.0.tar.gz
    > cd armadillo-3.910.0
    > mkdir build
    > cd build
    > cmake ../
    > sudo make install

    > git clone --recursive https://github.com/fadern/bayesic
    > cd bayesic
    > mkdir build
    > cd build
    > cmake ../

## Running

Bayesic currently supports three modes *-m bayes*, *-m logistic* and *-m loglinear*, corresponding to the models above. It is invoked by:

    > ./src/bayesic -m bayes /data/dataset.pair /data/dataset
    snp1        snp2        Posterior   N
    rs412512    rs516161    0.9602      3418
    rs51512     rs151251    0.1012      3412
    rs51512     rs516163    0.2312      3416
    > ./src/bayesic -m logistic /data/dataset.pair /data/dataset
    snp1        snp2        Beta        SE          P       N
    rs412512    rs516161    0.9524      0.2041      0.0301  3418
    rs51512     rs151251    0.2141      0.5120      0.4121  3412
    rs51512     rs516163    -0.2331     0.6001      0.7122  3416
    > ./src/bayesic -m loglinear /data/dataset.pair /data/dataset
    snp1        snp2        P           N
    rs412512    rs516161    0.0112      3418
    rs51512     rs151251    0.3521      3412
    rs51512     rs516163    0.4001      3416
    
where

    > ls /data/
    dataset.pair
    dataset.bim
    dataset.fam
    dataset.bed

    > cat /data/dataset.pair
    rs412512 rs516161
    rs51512 rs151251
    rs51512 rs516163

