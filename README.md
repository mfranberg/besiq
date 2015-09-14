# Besiq

Besiq is an application for analysing pairwise genetic interactions in genome-wide association studies. It provides several methods to choose from and is efficiently implemented in C++.

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

    > git clone --recursive https://github.com/fadern/besiq
    > cd besiq
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

    > ./src/besiq wald /data/dataset.pair /data/dataset
    snp1        snp2        LR  P   N
    rs412512    rs516161    4.1500    0.386 3418
    rs51512     rs151251    3.1200    0.5379    3412
    rs51512     rs516163    15.2661   0.0042    3416

## Dependencies

Besiq requires you to have the following packages installed:

* [CMake](http://www.cmake.org/cmake/resources/software.html) (A build system for C and C++)
* [Armadillo](http://arma.sourceforge.net/download.html) (C++ matrix library which in turn depends on BLAS and LAPACK)

If you install them locally (as on a cluster) just specify -DCMAKE_PREFIX_PATH to tell cmake where the libraries are.

## Installing

Before starting it is important to note when building armadillo you need to edit the config.hpp, and add support for
blas and lapack, otherwise you may experience linking errors when building besiq.

To build with dependencies without root access (note --recursive flag to get submodules):

    > mkdir ~/prefix

    > wget http://sourceforge.net/projects/arma/files/armadillo-3.910.0.tar.gz
    > tar -xf armadillo-3.910.0.tar.gz
    > cd armadillo-3.910.0
    > mkdir build
    > cd build
    > cmake ../ -DCMAKE_INSTALL_PREFIX=$HOME/prefix
    > make install

    > git clone --recursive https://github.com/fadern/besiq
    > cd besiq
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

    > git clone --recursive https://github.com/fadern/besiq
    > cd besiq
    > mkdir build
    > cd build
    > cmake ../
    > make install

## Running

Besiq is run by simply typing *besiq*. This will display a list of available subcommands:

* **pairs** - Used to generate a list of pairs that will be analyzed. The reason why such a list generated beforehand is to allow simple distributed calculation.
* **view** - Show the binary result files.
* **correct** - Perform multiple testing correction.

The following analysis types are available:

* **glm** - Uses a generalized linear model with either a binary or continuous phenotype, and possible covariates. This method is relatively slow because of the underlying iterative algorithm. This implementation uses the likelihood ratio test between the null and alternative to compute a p-value.
* **wald** - Uses a generalized linear model with either a binary or continuous phenotype, without covariates. This method is fast because of the closed form solutions. Inference is performed by a Wald-test which is asymptotically equivalent to the likelihood ratio test.
* **stagewise** - Increases power by considering the pairs in stages, fast, without covariates. Inference is performed by either a closed testing scheme or an approximate adaptive method.
* **scaleinv** - Tests multiple link functions on the data using a generalized linear model and reports a p-value for each.
* **loglinear** - Fast, powerful, but assumes that there is at most a single main effect. Preferebly used on data where the significant variants have been filtered out beforehand.
* **caseonly** - A test based on LD, interaction generates LD in case/control cohorts. Here the LD is estimated using the covariance between variants. The specific test used depends on the -m flag (see command for more info).
* **separate** - Codes the variants into either dominant or rescessive encoding, creates the 4 possible models, and tests each one using a generalized linear model.
    
The following commands are experimental:
    
* **bayes** - Compute the model posterior using the same models as in the loglinear method, using a conjugate dirichlet prior.
* **imputed** - Perform interaction fine-mapping analysis using impute2 imputed data from two (small) regions.
* **env** - Perform stage-wise gene-environment analysis.
    
### Genotype files

This software works with binary plink files .fam, .bim and .bam, and are specified using the path without the extension.

### Phenotype and covariate files

These files differs slighly from the plink phenotype and covariates file. Missing values are specified with NA. All binary phenotypes and covariates should be coded with 0/1 (in contrast to 1/2). The file format is

    FID IID pheno [pheno2 ... phenon]
    
or

    FID IID cov1 cov2 cov3 ...

# Stage-wise closed testing

Generate a list of all possible pairs with a combined maf of greater than 0.04, a marginal maf greater than 0.2 and with a distance of at least 1 Mbp between variants in the same pair.

    > besiq pairs -c 0.04 -m 0.2 -d 1000000 /data/dataset > dataset.pair

Run the stepwise command, this should preferebly be parallelized on a cluster with more than 100k variants in the genotype file. It is also recommended to filter out pairs with a p-value higher than 0.05 / num_pairs. This command will use the phenotype in the plink file, if in another file specify with -p.

    > besiq stagewise /data/dataset.pair /data/dataset > results.out
    
Use closed testing to perform multiple testing correction assuming 50 marginally associated variants and 100k variants.

    > besiq correct --weight 0.25,0.25,0.25,0.25 --num-tests 4999950000,5000000,5000000,1225 results.out /data/dataset
    
This will output all pairs significant on at least one scale, and the adjusted p-values.

### GLM, Wald and loglinear

These are all run similiarly. First pairs should be created

    > besiq pairs -c 0.04 -m 0.2 -d 1000000 /data/dataset > dataset.pair
    
Then these methods are simply run by the following. These commands will use the phenotype in the plink file, if in another file specify with -p.

    > besiq wald /data/dataset.pair /data/dataset > result.wald.out
    > besiq glm -f factor -l logistic /data/dataset.pair /data/dataset > results.logistic.out
    > besiq loglinear /data/dataset.pair /data/dataset > results.loglinear.out

# Evaluation

If you want to evaluate your own method, or the methods implemented in besiq under various simulation settings, then check out the Python packages [epibench](https://github.com/mfranberg/epibench) for benchmarking, and [epigen](https://github.com/mfranberg/epigen) for generating data.
