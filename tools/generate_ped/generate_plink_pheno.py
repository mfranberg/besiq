import argparse
import random
from plinkio import plinkfile
from plink import generate

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
def generate_phenotype(snp1, snp2, model):
    if snp1 != 3 and snp2 != 3:
        return int( random.random( ) <= model[ 3 * snp1 + snp2 ] )
    else:
        return None

##
# Writes the phenotypes for two snps and a set of individuals to a file.
#
# @param sample_list List of plinkio.plinkfile.Sample.
# @param snp1_row List of genotypes for the first snp.
# @param snp2_row List of genotypes for the second snp.
# @param model Penetrance matrix as a length 9 list
# @param output_file The phenotypes will be written to this file.
#
def write_phenotypes(sample_list, snp1_row, snp2_row, model, output_file):
    output_file.write( "FID\tIID\tPheno\n" )

    number_of_cases = number_of_controls = 0
    for sample, snp1, snp2 in zip( sample_list, snp1_row, snp2_row ):
        pheno = generate_phenotype( snp1, snp2, model )
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
# Generates two distinct loci.
#
# @param loci List of locus.
#
# @return A list of two loci.
#
def generate_two_row_numbers(loci):
    loci_index = list( range( len( loci ) ) )
    random.shuffle( loci_index )

    return loci_index[:2]

if __name__ == "__main__":
    arg_parser = argparse.ArgumentParser( description='Generates interaction pairs from a given model.' )
    arg_parser.add_argument( 'plink_file', metavar='plink_file', type=str, help='Path to a plink file.' )
    arg_parser.add_argument( '--model', metavar='model', nargs=9, type=generate.probability, help='Space-separated list of floating point numbers that represents the penetrance matrix, specified row-wise from left to right.', required = True )
    arg_parser.add_argument( '--out', metavar='output_file', help='Output phenotype file.', required = True )

    args = arg_parser.parse_args( )

    plink_file = plinkfile.open( args.plink_file )
    snp1_index, snp2_index = generate_two_row_numbers( plink_file.get_loci( ) )
    snp1_row = snp2_row = None
    for row_num, row in enumerate( plink_file ):
        if row_num == snp1_index:
            snp1_row = list( row )
        elif row_num == snp2_index:
            snp2_row = list( row )

    with open( args.out, "w" ) as output_file:
        number_of_cases, number_of_controls = write_phenotypes( plink_file.get_samples( ), snp1_row, snp2_row, args.model, output_file )
        print "Wrote {0} cases and {1} controls".format( number_of_cases, number_of_controls ) 
