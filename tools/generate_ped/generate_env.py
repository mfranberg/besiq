import argparse
import random
from plinkio import plinkfile
from plink import generate
from math import sqrt

##
# Simple class for storing the parameters
# of the simulation.
#
class EnvParams:
    ##
    # @param env The environment factor.
    # @param model The mean level for each combination.
    # @param std The standard deviation.
    #
    def __init__(self, env, model, std):
        self.env = env
        self.model = model
        self.std = std

##
# Generates a phenotype for the given snp pair and penetrance
# matrix. 
#
# @param snp Genotype of variant.
# @param env Environmental factor.
# @param model Mean for each genotype and environmental level.
# @param std Standard deviation of phenotype.
#
# @return A continuous phenotype, or "NA" if genotype was missing.
#
def generate_phenotype(snp, env, model, std):
    if snp != 3:
        return random.normalvariate( model[ env * 3 + snp ], std )
    else:
        return "NA"

##
# Generates an environment variable
# from the level frequencies.
#
# @param level frequencies that sum to one.
#
# @return The generated level.
#
def generate_environment(env):
    u = random.random( )
    cumsum = 0.0
    for i, e in enumerate( env ):
        cumsum += e
        if u <= cumsum:
            return i

    raise Exception( "Environmental frequencies does not sum to 1." )

##
# Finds a variant and returns a name to it.
#
# @param loci A list of variants.
# @param snp The name of a variant.
#
# @return The index of the given variant.
#
def find_variant(loci, snp):
    snp_to_index = dict( zip( map( lambda x: x.name, loci ), range( len( loci ) ) ) )
    return snp_to_index[ snp1 ]

##
# Writes the phenotypes for two snps and a set of individuals to a file.
#
# @param sample_list List of plinkio.plinkfile.Sample.
# @param snp_row List of genotypes for the first snp.
# @param params Parameters for generation
# @param env_file Output file for environmental factors.
# @param pheno_file Output file for phenotype.
#
def write_phenotypes(sample_list, snp_row, params, env_file, pheno_file):
    env_file.write( "FID\tIID\tEnv\n" )
    pheno_file.write( "FID\tIID\tPheno\n" )

    for sample, snp in zip( sample_list, snp_row ):
        env = generate_environment( params.env )
        pheno = generate_phenotype( snp, env, params.model, params.std )

        env_file.write( "{0}\t{1}\t{2}\n".format( sample.fid, sample.iid, str( env ) ) )
        pheno_file.write( "{0}\t{1}\t{2}\n".format( sample.fid, sample.iid, str( pheno ) ) )

##
# Generates the index of random variant.
#
# @param loci List of locus.
#
# @return The index of the random variant.
#
def generate_random_variant(loci):
    loci_index = list( range( len( loci ) ) )
    random.shuffle( loci_index )

    return loci_index[ 0 ]

if __name__ == "__main__":
    arg_parser = argparse.ArgumentParser( description='Generates a phenotype and environmental factor.' )
    arg_parser.add_argument( 'plink_file', metavar='plink_file', type=str, help='Path to a plink file.' )
    arg_parser.add_argument( '--model', metavar='model', nargs='+', type=float, help='Space-separated list of mean for each cell, columns correspond to the genotype.', required = True )
    arg_parser.add_argument( '--std', metavar='std', type=float, help='Common standard deviation.', default = 1.0 )
    arg_parser.add_argument( '--env', metavar='env', nargs='+', type=float, help='Frequency of each environmental level, must sum to 1.', default = None )
    arg_parser.add_argument( '--variant', metavar='variant', type=str, help='Name of the variant the phenotype should be based on.', default = None )
    arg_parser.add_argument( '--pheno-file', metavar='pheno_file', help='Output phenotype file.', required = True )
    arg_parser.add_argument( '--env-file', metavar='env_file', help='Output environment file.', required = True )

    args = arg_parser.parse_args( )

    plink_file = plinkfile.open( args.plink_file ) 
    snp_index = generate_random_variant( plink_file.get_loci( ) )
    if args.variant:
        snp_index = find_variant( plink_file.get_loci( ), args.variant )

    print plink_file.get_loci( )[ snp_index ].name

    snp_row = None
    for row_num, row in enumerate( plink_file ):
        if row_num == snp_index:
            snp_row = list( row )
            break

    env = args.env
    if not env:
        num_env = len( args.model ) / 3
        env = [ 1.0 / num_env ] * num_env

    params = EnvParams( env, args.model, args.std )
    env_file = open( args.env_file, "w" )
    pheno_file = open( args.pheno_file, "w" )
    write_phenotypes( plink_file.get_samples( ), snp_row, params, env_file, pheno_file )
