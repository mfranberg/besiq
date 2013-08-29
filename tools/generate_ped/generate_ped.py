import numpy as np
import argparse
import os
from math import sqrt

from plink.tfam import TFamFile
from plink.tped import TPedFile
from plink.pair import PairFile

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
# Computes the joint Hardy-Weinberg model represented
# as a vector.
#
# @param maf Desired minor allele frequncy.
# @param ld The probability Pr[x2 = i| x1 = i].
#
# @return Joint probability model for the Hardy-Weinberg
#         model.
#
def joint_maf(maf, ld):
    maf1 = maf[ 0 ]
    maf2 = maf[ 1 ]

    p = [ ( 1 - maf[ 0 ] )**2, 2 * maf[ 0 ] * ( 1 - maf[ 0 ] ), ( maf[ 0 ] )**2 ]
    q = [ ( 1 - maf[ 1 ] )**2, 2 * maf[ 1 ] * ( 1 - maf[ 1 ] ), ( maf[ 1 ] )**2 ]

    if not ld:
        return   [ p[ 0 ] * q[ 0 ], p[ 0 ] * q[ 1 ], p[ 0 ] * q[ 2 ],
                   p[ 1 ] * q[ 0 ], p[ 1 ] * q[ 1 ], p[ 1 ] * q[ 2 ],
                   p[ 2 ] * q[ 0 ], p[ 2 ] * q[ 1 ], p[ 2 ] * q[ 2 ] ]
    else:
        return   [ p[ 0 ] * ld, p[ 0 ] * ( 1 - ld ) / 2.0, p[ 0 ] * ( 1 - ld ) / 2.0,
                   p[ 1 ] * ( 1 - ld ) / 2.0, p[ 1 ] * ld, p[ 1 ] * ( 1 - ld ) / 2.0,
                   p[ 2 ] * ( 1 - ld ) / 2.0, p[ 2 ] * ( 1 - ld ) / 2.0, p[ 2 ] * ld ]

    return joint_prob

##
# Computes the probability of a pair of snps x1, x2
# given the disease status.
#
# @param penetrance Penetrance of each genotype.
# @param maf Desired minor allele frequncy.
# @param ld The probability Pr[x2 = i| x1 = i].
#
# @return A tuple containing first the probabilities of 
#         x1, x2 | D = 1, and second the probabilityies
#         of x1, x2 | D = 2 represented as vectors.
#
def joint_snp(penetrance, maf, ld):
    joint_hw = np.array( joint_maf( maf, ld ) )

    penetrance = np.array( penetrance )

    joint_numerator = joint_hw * penetrance
    denom = sum( joint_numerator )
    sick_prob = ( joint_numerator / denom ).tolist( )

    joint_numerator = joint_hw * ( 1 - penetrance )
    denom = sum( joint_numerator )
    healthy_prob = ( joint_numerator / denom ).tolist( )

    return ( sick_prob, healthy_prob )

##
# Samples the given number of samples from
# the joint probability model.
#
# @param joint_prob Joint probability model represented as 
#                   a vector.
# @param num_samples Number of samples to take.
#
# @return Iterator over all samples.
#
def sample(joint_prob, num_samples):
    samples = np.random.multinomial( num_samples, joint_prob, 1 )[ 0 ]
    for i in range( 3 ):
        for j in range( 3 ):
            for k in range( samples[ 3 * i + j ] ):
                yield (i, j)

##
# Generates two lists of genotypes for two SNPs under
# the given joint probability model.
#
# @param joint_prob The probability model.
# @param num_samples Number of samples.
#
# @return A tuple containing first the genotypes of the
#         first snp, and second the genotypes of the second
#         snp.
#
def generate_genotypes(joint_prob, num_samples):
    TO_ALLELE = { 0 : 'A A', 1 : 'A T', 2 : 'T T' }
    snp1_list = [ 0 ] * num_samples
    snp2_list = [ 0 ] * num_samples
    i = 0
    for snp1, snp2 in sample( joint_prob, num_samples ):
        snp1_list[ i ] = TO_ALLELE[ snp1 ]
        snp2_list[ i ] = TO_ALLELE[ snp2 ]

        i += 1

    return ( snp1_list, snp2_list )

##
# Writes the plink data to the given .tped, .tfam and .pair files.
#
# @param args Parsed arguments.
# @param tped_file TPedFile object.
# @param tfam_file TFamFile object.
# @param pair_file PairFile object.
#
def write_genotypes(args, tped_file, tfam_file, pair_file):
    sick_prob, healthy_prob = joint_snp( args.model, args.maf, args.ld )

    for i in range( 0, args.npairs * 2, 2 ):
        snp1_sick, snp2_sick = generate_genotypes( sick_prob, args.ncases )
        snp1_healthy, snp2_healthy = generate_genotypes( healthy_prob, args.ncontrols )
    
        tped_file.write( snp1_sick + snp1_healthy )
        tped_file.write( snp2_sick + snp2_healthy )
        pair_file.write( i, i + 1 )
    
    for i in range( args.ncases + args.ncontrols ):
        tfam_file.write( int( i < args.ncases ) + 1 )

##
# Writes the plink data in the location specified by the
# given arguments.
#
# @param args Parsed arguments.
#
def write_data(args):
    path, ext = os.path.splitext( args.out )
    tped_file = TPedFile( path + ".tped" )
    tfam_file = TFamFile( path + ".tfam" )
    pair_file = PairFile( path + ".pair" )
    
    write_genotypes( args, tped_file, tfam_file, pair_file )

    tped_file.close( )
    tfam_file.close( )
    pair_file.close( )
    
    status = os.system( "plink --tfile {0} --make-bed --out {1} > /dev/null".format( path, path ) )
    if status == -1:
        print( "Could not run plink, is it installed?" )
        exit( 1 )


if __name__ == "__main__":
    arg_parser = argparse.ArgumentParser( description='Generates interaction pairs from a given model.' )
    arg_parser.add_argument( '--model', metavar='model', nargs=9, type=probability, help='Space-separated list of floating point numbers that represents the penetrance matrix, specified row-wise from left to right.', required = True )
    arg_parser.add_argument( '--maf', metavar='maf', nargs=2, type=probability, help='Minor allele frequency of the two snps.', required = True )
    arg_parser.add_argument( '--ncases', metavar='ncases', type=int, help='Number of cases.', required = True )
    arg_parser.add_argument( '--ncontrols', metavar='ncontrols', type=int, help='Number of controls.', required = True )
    arg_parser.add_argument( '--npairs', metavar='npairs', type=int, help='Number of interaction pairs', required = True )
    arg_parser.add_argument( '--ld', metavar='ld', type=probability, help='Strength of LD (ignores second maf).' )
    
    arg_parser.add_argument( '--out', metavar='output_file', help='Output .tped file.', required = True )

    args = arg_parser.parse_args( )
    write_data( args )
