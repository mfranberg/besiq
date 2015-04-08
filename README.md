# Bayesic

Bayesic is an application for analysing pairwise genetic interactions in genome-wide association studies. It provides several methods to choose from and is efficiently implemented in C++.

The following methods are currently implemented for gene-gene tests:

* Stage-wise closed testing (efficient and powerful on a genome-wide scale)
* Bernoulli GLM regression with LR test (supports covariates).
* Logistic regression with Wald test (fast, closed form).
* Loglinear model with analytic test (does not support covariates or model main effects, but is very fast).
* Linear model for quantitative traits.
* Stage-wise closed testing for quantitative traits (experimental).
* Bayesian model posterior (experimental).

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

    > ./src/bayesic wald /data/dataset.pair /data/dataset
    snp1        snp2        LR  P   N
    rs412512    rs516161    4.1500    0.386 3418
    rs51512     rs151251    3.1200    0.5379    3412
    rs51512     rs516163    15.2661   0.0042    3416

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
    > make install

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
    > make install

## Running

### Genotype files

This software works with binary plink files .fam, .bim and .bam, and are specified using the path without the extension.

### Phenotype and covariate files

Missing values are specified with NA. All binary phenotypes and covariates should be coded with 0/1. The file format is

    FID IID pheno
    
or

    FID IID cov1 cov2 cov3 ...

# Stage-wise closed testing

Generate a list of all possible pairs with a combined maf of greater than 0.04, a marginal maf greater than 0.2 and with a distance of at least 1 Mbp between variants in the same pair.

    > bayesic pairs -c 0.04 -m 0.2 -d 1000000 /data/dataset > dataset.pair

Run the stepwise command, this should preferebly be parallelized on a cluster with more than 100k variants in the genotype file. It is also recommended to filter out pairs with a p-value higher than 0.05 / num_pairs. This command will use the phenotype in the plink file, if in another file specify with -p.

    > bayesic stagewise /data/dataset.pair /data/dataset > results.out
    
Use closed testing to perform multiple testing correction assuming 50 marginally associated variants and 100k variants.

    > bayesic correct --weight 0.25,0.25,0.25,0.25 --num-tests 4999950000,5000000,5000000,1225 results.out /data/dataset
    
This will output all pairs significant on at least one scale, and the adjusted p-values.

### GLM, Wald and loglinear

These are all run similiarly. First pairs should be created

    > bayesic pairs -c 0.04 -m 0.2 -d 1000000 /data/dataset > dataset.pair
    
Then these methods are simply run by the following. These commands will use the phenotype in the plink file, if in another file specify with -p.

    > bayesic wald /data/dataset.pair /data/dataset > result.wald.out
    > bayesic glm -f factor -l logistic /data/dataset.pair /data/dataset > results.logistic.out
    > bayesic loglinear /data/dataset.pair /data/dataset > results.loglinear.out

# Evaluation

If you want to evaluate your own method, or the methods implemented in bayesic under various simulation settings, then check out the Python packages [epibench](https://github.com/mfranberg/epibench) for benchmarking, and [epigen](https://github.com/mfranberg/epigen) for generating data.
