import argparse
import random
from plinkio import plinkfile
from plink import generate
from math import sqrt, exp

##
# Generates a phenotype for the given snp pair and penetrance
# matrix. 
#
# @param snp1 Genotype of first snp
# @param snp2 Genotype of second snp
# @param model Penetrance matrix as a length 9 list
#
# @return 1 or 0 representing case control if no snp was missing, None
#         otherwise.
#
def generate_phenotype(variants, beta0, beta):
    if 3 in variants:
        return None

    score = beta0 + sum( [ b * v for b, v in zip( beta, variants ) ] )
    p = 1/(1 + exp(-score))
    
    return int( random.random( ) <= p )

def find_rows(plink_file, loci):
    loci_set = set( loci )
    rows = [ ]
    for i, row in enumerate( plink_file ):
        if i in loci_set:
            rows.append( list( row ) )

    return rows

def compute_maf(row):
    no_missing = filter( lambda x: x != 3, row )
    p = sum( no_missing ) / ( 2.0 * len( no_missing ) )

    return p

def find_beta0(rows, beta):
    maf = [ compute_maf( r ) for r in rows ]
    beta0 = -sum( b * 2 * m for b, m in zip( beta, maf ) )

    return beta0

##
# Writes the phenotypes for two snps and a set of individuals to a file.
#
# @param sample_list List of plinkio.plinkfile.Sample.
# @param snp1_row List of genotypes for the first snp.
# @param snp2_row List of genotypes for the second snp.
# @param model Penetrance matrix as a length 9 list
# @param output_file The phenotypes will be written to this file.
#
def write_phenotypes(sample_list, rows, beta, output_file, beta0 = None):
    output_file.write( "FID\tIID\tPheno\n" )

    if not beta0:
        beta0 = find_beta0( rows, beta )

    number_of_cases = number_of_controls = 0
    for i, sample in enumerate( sample_list ):
        variants = [ rows[ j ][ i ] for j in range( len( rows ) ) ]
        
        pheno = generate_phenotype( variants, beta0, beta )
        pheno_str = str( pheno )
        if pheno == 0:
            number_of_controls += 1
        elif pheno == 1:
            number_of_cases += 1
        else:
            pheno_str = "NA"

        output_file.write( "{0}\t{1}\t{2}\n".format( sample.fid, sample.iid, pheno_str ) )

    return number_of_cases, number_of_controls

##
# Generates a set of loci.
#
# @param loci List of locus.
#
# @return The indices of the selected loci.
#
def generate_loci_set(loci, n):
    loci_index = list( range( len( loci ) ) )
    random.shuffle( loci_index )

    return loci_index[:n]

##
# Generates effect sizes according to the effect size
# distribution.
#
# @param n The number of beta
# @param mean The mean of the beta distribution
# @param sd The standard deviation of the beta distibution
#
# @return list of effect sizes
#
def generate_beta(n, mean, sd):
    return [ random.normalvariate( mean, sd ) for i in range( n ) ]

if __name__ == "__main__":
    arg_parser = argparse.ArgumentParser( description='Generates interaction pairs from a given model.' )
    arg_parser.add_argument( 'plink_file', metavar='plink_file', type=str, help='Path to a plink file.' )
    arg_parser.add_argument( '--beta0', metavar='beta0', type=float, help='Sets the intercept, by default it is chosen to get 50/50 cases and controls.', default = None )
    arg_parser.add_argument( '--beta', metavar='beta', nargs=2, type=float, help='The mean and variance of the beta variables.', default = [ 0.2, 0.3 ] )
    arg_parser.add_argument( '--num-loci', metavar='num_loci', type=int, help='The number of loci that is involved in the phenotype.', default = 10 )
    arg_parser.add_argument( '--out', metavar='output_file', help='Output phenotype file.', required = True )

    args = arg_parser.parse_args( )

    plink_file = plinkfile.open( plink_file ) 
    loci = plink_file.get_loci( )
    snp_indices = generate_loci_set( loci, num_loci )
    beta = generate_beta( num_loci, beta[ 0 ], beta[ 1 ] )
    rows = find_rows( plink_file, snp_indices )

    print " ".join( [ loci[ i ].name for i in snp_indices ] )

    with open( out, "w" ) as output_file:
        number_of_cases, number_of_controls = write_phenotypes( plink_file.get_samples( ), rows, beta, output_file, beta0 = beta0 )
        print "Wrote {0} cases and {1} controls".format( number_of_cases, number_of_controls ) 
