import numpy as np
import argparse
import os
from math import sqrt, exp
from functools import partial

from plink.ped import PedFile
from plink.map import MapFile
from plink.pair import PairFile
from plink.cov import CovFile

##
# Represents a covariate in a logistic regression model.
#
class covariate:
    ##
    # Constructor.
    #
    # @param name Name of the covariate.
    # @param beta Beta of the covariate.
    # @param dist Function that generates samples for this covaraite.
    # @param mean Mean of this coviarate.
    #
    def __init__(self, name, beta, dist, mean):
        self.name = name
        self.beta = beta
        self.dist = dist
        self.mean = mean

##
# Converts the given value to a floating point value
# if it represents a probability, otherwise it raises
# a ValueError.
#
# @param prob_string String representing a probability.
#
# @return String as a floating point probability.
#
def probability(prob_string):
    prob_value = float( prob_string )
    if 0.0 <= prob_value <= 1.0:
        return prob_value
    else:
        raise ValueError( "Probability not between 0 and 1." )

    return prob_value

##
# Parses a coviarate string on the form
# "name,beta,distribution,param1,param2" and returns
# a covariate object for it.
#
# @param cov_str A string representing a covariate.
#
# @return A covariate object for the given covariate string.
#
def parse_cov(cov_str):
    params = cov_str.strip( ).split( "," )

    name = params[ 0 ]
    beta = float( params[ 1 ] )

    mean = None
    dist = None
    if params[ 2 ] == "normal":
        mean = float( params[ 3 ] )
        sd = float( params[ 4 ] )

        dist = partial( np.random.normal, mean, sd )
    elif params[ 2 ] == "binomial":
        n = int( params[ 3 ] )
        p = float( params[ 4 ] )
        mean = n * p

        dist = partial( np.random.binomial, n, p )

    return covariate( name, beta, dist, mean )

##
# Returns a beta0 such that individuals in the logistic
# regression model will have 0.5 chance of becoming a case.
#
# @param covariates List of covariates.
#
# @return A beta0 such that the probability is 0.5 of becoming
#         a case.
#
def get_beta0(covariates):
    cov_mean = sum( c.mean * c.beta for c in covariates )

    return -cov_mean

##
# Generates a SNP with the given maf according to
# Hardy-Weinburg equilibrium.
#
# @param maf Minor allele frequncy of the SNP.
# 
# @return The generated SNP.
#
def generate_genotype(maf):
    return np.random.binomial( 2, maf )

##
# Generates the risk from snps.
#
# @param snp1 Genotype at locus 1.
# @param snp2 Genotype at locus 2.
# @param model Model for genotype risk at two locus.
#
# @return The risk for the genotypes.
#
def generate_genotype_risk(snp1, snp2, model):
    return model[ snp1 * 3 + snp2 ]

##
# Generates a list of covariates according to the
# specified distribution.
#
# @param covariates List of covariates.
#
# @return A list of generated covariates.
#
def generate_covariates(covariates):
    return [ c.dist( ) for c in covariates ]

##
# Generates the risk from covariates.
#
# @param beta0 Beta0 of the model.
# @param cov Generated covariates.
# @param cov_beta Beta for covariates.
#
# @return The risk from the covariates.
#
def generate_covariate_risk(beta0, cov, cov_beta):
    if not cov:
        return 0.0

    score = beta0
    score += sum( b * c for b, c in zip( cov_beta, cov ) )
    p = 1 / ( 1 + exp( -score ) )

    return p

##
# Generates a phenotype under the assumption that genotypes
# and environmental factors act as sufficient independent
# causes.
#
# @param genotype_risk Risk from genotypes.
# @param covariate_risk Risk from covariates.
#
# @return 0 for controls and 1 for cases.
#
def generate_phenotype(genotype_risk, covariate_risk):
    p = 1 - ( 1 - genotype_risk ) * ( 1 - covariate_risk )

    return np.random.binomial( 1, p )

##
# Converts a snp to an allele string.
#
# 0 -> "A A"
# 1 -> "A C"
# 2 -> "C C"
#
# @return The corresponding allele string.
#
def to_allele(snp):
    if snp == 0:
        return "A A"
    elif snp == 1:
        return "A C"
    else:
        return "C C"

##
# Writes the plink data to the given .tped, .tfam and .pair files.
#
# @param args Parsed arguments.
# @param ped_file PedFile object.
# @param map_file MapFile object.
# @param pair_file PairFile object.
# @param cov_file CovFile object.
#
def write_genotypes(args, ped_file, map_file, pair_file, cov_file):
    beta0 = get_beta0( args.covariate )

    num_controls = 0
    num_cases = 0

    while num_cases < args.ncases or num_controls < args.ncontrols:
        snp1 = generate_genotype( args.maf[ 0 ] )
        snp2 = generate_genotype( args.maf[ 1 ] )
        snp_risk = generate_genotype_risk( snp1, snp2, args.model )
        
        covariates = generate_covariates( args.covariate )
        covariate_risk = generate_covariate_risk( beta0, covariates, [ c.beta for c in args.covariate ] )

        phenotype = generate_phenotype( snp_risk, covariate_risk )

        if phenotype == 0 and num_controls < args.ncontrols:
            num_controls += 1
        elif phenotype == 1 and num_cases < args.ncases:
            num_cases += 1
        else:
            continue

        ped_file.write( [ to_allele( snp1 ), to_allele( snp2 ) ], phenotype + 1 )

        if covariates:
            cov_file.write( covariates )
    
    pair_file.write( )
    map_file.write( )
    map_file.write( )

##
# Writes the plink data in the location specified by the
# given arguments.
#
# @param args Parsed arguments.
#
def write_data(args):
    path, ext = os.path.splitext( args.out )
    ped_file = PedFile( path + ".ped" )
    map_file = MapFile( path + ".map" )
    pair_file = PairFile( path + ".pair" )

    cov_file = None
    if args.covariate:
        cov_file = CovFile( path + ".cov", [ c.name for c in args.covariate ] )
    
    write_genotypes( args, ped_file, map_file, pair_file, cov_file )

    ped_file.close( )
    map_file.close( )
    pair_file.close( )

    if args.covariate:
        cov_file.close( )
    
    status = os.system( "plink --file {0} --make-bed --out {1} > /dev/null".format( path, path ) )
    if status == -1:
        print( "Could not run plink, is it installed?" )
        exit( 1 )


if __name__ == "__main__":
    arg_parser = argparse.ArgumentParser( description='Generates interaction pairs from a given model.' )
    arg_parser.add_argument( '--maf', metavar='maf', nargs=2, type=float, help='Minor allele frequency of the two snps.', required = True )
    arg_parser.add_argument( '--ncases', metavar='ncases', type=int, help='Number of cases.', required = True )
    arg_parser.add_argument( '--ncontrols', metavar='ncontrols', type=int, help='Number of controls.', required = True )
    arg_parser.add_argument( '--model', metavar='model', nargs=9, type=probability, help='Space-separated list of floating point numbers that represents the penetrance matrix, specified row-wise from left to right.', required = True )
    arg_parser.add_argument( '--covariate', metavar='covariate', nargs='+', type=parse_cov, help='List of name,beta,distribution,param1,param2. Valid distributions are "normal" or "binomial".', default = [] ) 
    arg_parser.add_argument( '--out', metavar='output_file', help='Output prefix for .bim, .bed, .fam, .pairs and .cov.', required = True )

    args = arg_parser.parse_args( )
    write_data( args )
