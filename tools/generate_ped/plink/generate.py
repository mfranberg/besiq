import numpy as np
import argparse
import os
from math import sqrt

from .output import OutputFiles

##
# The parameters that does not changed between models.
# 
class FixedParams:
    def __init__(self, maf, ld, ncases, ncontrols):
        self.maf = maf
        self.ld = ld
        self.ncases = ncases
        self.ncontrols = ncontrols

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
# Writes the plink data to the given .tped, .case and .pair files.
#
# @param num_pairs The number of pairs to generate.
# @param penetrance The model penetrance.
# @param is_case Determines whether this is an interaction or not.
# @param fixed_params The simulation parameters.
# @param model_index Index of the model.
# @param output_file The output files.
#
def write_genotypes(num_pairs, penetrance, is_case, fixed_params, model_index, output_files):
    sick_prob, healthy_prob = joint_snp( penetrance, fixed_params.maf, fixed_params.ld )

    for i in range( num_pairs ):
        snp1_sick, snp2_sick = generate_genotypes( sick_prob, fixed_params.ncases )
        snp1_healthy, snp2_healthy = generate_genotypes( healthy_prob, fixed_params.ncontrols )
    
        output_files.tped_file.write( snp1_sick + snp1_healthy )
        output_files.tped_file.write( snp2_sick + snp2_healthy )
        output_files.pair_file.write( )
        output_files.case_file.write( is_case )
        output_files.model_file.write( model_index )

##
# Writes the indivduals to the given .tfam file, assuming that
# the cases are first followed by the controls.
#
# @param ncases The number of cases.
# @param ncontrols The number of controls.
# @param tfam_file An opened .tfam file.
#
def write_individuals(ncases, ncontrols, tfam_file):
    for i in range( ncases + ncontrols ):
        tfam_file.write( int( i < ncases ) + 1 )

##
# Writes the plink data in the location specified by the
# given arguments.
#
# @param fixed_params The simulation parameters.
# @param models The list of models to generate from.
# @param output_prefix The output prefix, different file endings will be generated.
#
def write_data(fixed_params, models, output_prefix):
    path, ext = os.path.splitext( output_prefix )
    output_files = OutputFiles( path )
  
    model_index = 1
    for num_pairs, params, is_case in models:
        write_genotypes( num_pairs, params, is_case, fixed_params, model_index, output_files )
        model_index += 1
    
    write_individuals( fixed_params.ncases, fixed_params.ncontrols, output_files.tfam_file )
    output_files.close( )

    status = os.system( "plink --tfile {0} --make-bed --out {1} > /dev/null".format( path, path ) )
    if status == -1:
        print( "Could not run plink, is it installed?" )
        exit( 1 )
