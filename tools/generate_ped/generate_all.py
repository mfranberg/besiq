import argparse
import itertools
from functools import partial

from plink import generate

##
# Class that generates all possible interactions under
# a given null model. It does so by generating all possible
# marginal vectors (p1, p2, p3) at two loci, and  combines
# them with the given null function to yield a complete
# penetrance (p1, p2, p3, p4, p5, p6, p7, p8, p9).
#
class InteractionGenerator:
    ##
    # @param mat_add The null model to use.
    #
    def __init__(self, mat_add):
        self.mat_add = mat_add

    ##
    # Converts a marginal null vector to a
    # complete penetrance matrix for two loci.
    #
    # e.g.
    #
    # (p1, p2, p3) -> ( p1, p2, p3, p1, p2, p3, p1, p2, p3 )
    #
    # or
    #
    # (p1, p2, p3) -> ( p1, p1, p1, p2, p2, p2, p3, p3, p3 )
    #
    # @param a The marginal vector. 
    # @param other If false each row has the same penetrance,
    #              otherwise each column.
    #
    # @return The penetrance matrix as a vector.
    #
    def to_mat(self, a, other = False):
        l = list( )
        for i in range(3):
            for j in range(3):
                index = i
                if other:
                    index = j

                l.append( a[ index ] )

        return tuple( l )

    ##
    # Generates all possible interaction and null models by
    # generating all possible marginal vectors and combining
    # them with the given null function, to see which null 
    # models can be reached. The interactions are then defined
    # as all models minus these null models.
    #
    def generate(self):
        all_mats = set( a for a in itertools.product( (0,1), repeat = 9 ) )
        null_mats = set( )
        for a in itertools.product( (0,1), repeat = 3 ):
            for b in itertools.product( (0,1), repeat = 3 ):
                a_mat = self.to_mat( a, True )
                b_mat = self.to_mat( b, True )
                null_mats.add( self.mat_add( a_mat, b_mat ) )

                a_mat = self.to_mat( a, True )
                b_mat = self.to_mat( b, False )
                null_mats.add( self.mat_add( a_mat, b_mat ) )
                
                a_mat = self.to_mat( a, False )
                b_mat = self.to_mat( b, True )
                null_mats.add( self.mat_add( a_mat, b_mat ) )

                a_mat = self.to_mat( a, False )
                b_mat = self.to_mat( b, False )
                null_mats.add( self.mat_add( a_mat, b_mat ) )

        int_mats = all_mats.difference( null_mats )

        return ( int_mats, null_mats )


##
# The OR null model, if either locus is one the disease
# develops.
#
def mat_or(a, b):
    return tuple( map( lambda x: x[ 0 ] | x[ 1 ], zip (a,b) ) )


##
# Computes the heritability of a given penetrance and minor
# allele frequency.
#
# @param penetrance 9-element vector of disease penetrance.
# @param maf 2-element vector of minor allele frequencies.
#
# @return The narrow sense heritability.
#
def heritability(penetrance, maf):
    maf1 = maf[ 0 ]
    maf2 = maf[ 1 ]

    p = [ ( 1 - maf[ 0 ] )**2, 2 * maf[ 0 ] * ( 1 - maf[ 0 ] ), ( maf[ 0 ] )**2 ]
    q = [ ( 1 - maf[ 1 ] )**2, 2 * maf[ 1 ] * ( 1 - maf[ 1 ] ), ( maf[ 1 ] )**2 ]

    joint_maf =  [ p[ 0 ] * q[ 0 ], p[ 0 ] * q[ 1 ], p[ 0 ] * q[ 2 ],
                   p[ 1 ] * q[ 0 ], p[ 1 ] * q[ 1 ], p[ 1 ] * q[ 2 ],
                   p[ 2 ] * q[ 0 ], p[ 2 ] * q[ 1 ], p[ 2 ] * q[ 2 ] ]

    pop_p = sum( p * m for p, m in zip( penetrance, joint_maf ) )

    h = 0.0
    for i in range( 3 ):
        for j in range( 3 ):
            cell = 3 * i + j
            h += ( penetrance[ cell ] - pop_p )**2 * joint_maf[ cell ]

    return h / ( pop_p * ( 1 - pop_p ) )

##
# Given a desired heritability, a fully penetrant disease model
# (penetrance is either 1 or 0), this function tries to find a
# non-fully penetrant disease model that has the desired
# heritability.
#
# @param desired_heritability The desired heritability
# @param base_risk The population risk.
# @param maf Minor allele frequency.
# @param model The fully penetrance disease model.
# 
# @return A non-fully penetrant disease model with as close
#         heritability to the desired as possible.
#
def find_penetrance(desired_heritability, base_risk, maf, model):
    if sum( model ) == 0 or sum( model ) == 9:
        return [ 0.5 ] * 9

    step = 0.01
    disease_p = base_risk + step
    disease_map = dict( { 0 : base_risk, 1 : 1.0 } )
    penetrance = list( map( lambda x: disease_map[ x ], model ) )
    if heritability( penetrance, maf ) < desired_heritability:
        print( "Warning: impossible to find penetrance with desired heritability under (will use maximum): " + str( model ) )

    while disease_p <= 1.0:
        disease_map = dict( { 0 : base_risk, 1 : disease_p } )
        penetrance = list( map( lambda x: disease_map[ x ], model ) )
        if heritability( penetrance, maf ) >= desired_heritability:
            return penetrance

        disease_p += step

    return penetrance

if __name__ == "__main__":
    arg_parser = argparse.ArgumentParser( description='Generates all possible interaction pairs.' )
    arg_parser.add_argument( '--maf', metavar='maf', nargs=2, type=generate.probability, help='Minor allele frequency of the two snps.', required = True )
    arg_parser.add_argument( '--ncases', metavar='ncases', type=int, help='Number of cases.', required = True )
    arg_parser.add_argument( '--ncontrols', metavar='ncontrols', type=int, help='Number of controls.', required = True )
    arg_parser.add_argument( '--ld', metavar='ld', type=generate.probability, help='Strength of LD (ignores second maf).' )
    arg_parser.add_argument( '--num-pairs', metavar='num_pairs', type=int, help='Number of pairs to generate from each model.' )
    arg_parser.add_argument( '--heritability', metavar='heritability', type=float, help='Approximate heritability of each model.', default = 0.02 )
    arg_parser.add_argument( '--base-risk', metavar='base_risk', type=float, help='The base risk of the neutral alleles.', default = 0.5 )
    
    arg_parser.add_argument( '--out', metavar='output_file', help='Output .tped file.', required = True )

    args = arg_parser.parse_args( )
    
    # Generate interaction models
    generator = InteractionGenerator( mat_or )
    interactions, nulls = generator.generate( )

    partial_find_penetrance = partial( find_penetrance, args.heritability, args.base_risk, args.maf )
    interaction_penetrances = list( map( partial_find_penetrance, interactions ) )

    models = [ ( args.num_pairs, p, 1 ) for p in interaction_penetrances ]

    fixed_params = generate.FixedParams( args.maf, args.ld, args.ncases, args.ncontrols )
    generate.write_data( fixed_params, models, args.out )
