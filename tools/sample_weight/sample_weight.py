import sys
import os
import subprocess
import argparse
import random
from collections import defaultdict, Counter

def main():
    parser = argparse.ArgumentParser( description='Sample weights and determine frequency of significant pairs.' )
    parser.add_argument( 'result_file', type=str, help='The output file from bayesic -m stepwise.' )
    parser.add_argument( 'plink_file', type=str, help='The same plink file that was used to create the result file.' )
    parser.add_argument( 'output_file', type=str, help='The output prefix.' )
    parser.add_argument( '--args', type=str, help='Argument string to closed_correction.py.', default = "" )
    parser.add_argument( '--weight-file', type=str, help='Contains the weights to run', required = True )
    args = parser.parse_args( )

    closed_path = os.path.join( os.path.dirname( __file__ ), "../closed_correction/closed_correction.py" )

    weight_file = open( args.weight_file, "r" )
    weight_list = [ line.strip( ) for line in weight_file ]

    seen_count = Counter( )
    significant_count = Counter( )
    for w in weight_list:
        cmd = " ".join( [ "python", closed_path, "--weight", w, args.args, args.result_file, args.plink_file  ] )

        output_path = args.output_file
        with open( output_path, "w" ) as output_file:
            subprocess.check_call( cmd, shell = True, stdout = output_file )

        with open( output_path, "r" ) as output_file:
            try:
                next( output_file ) # Skip header
            except StopIteration:
                continue

            for line in output_file:
                column = line.strip( ).split( )
                seen_count[ column[ 0 ] + " " + column[ 1 ] ] += 1

                if column[ 7 ] != "NA" and float( column[ 7 ] ) <= 0.05:
                    significant_count[ column[ 0 ] + " " + column[ 1 ] ] += 1

        print seen_count

    all_pairs = set( seen_count.keys( ) ).union( set( significant_count.keys( ) ) )
    for pair in all_pairs:
        print "{0}\t{1}\t{2}".format( pair, float( seen_count[ pair ] ) / len( weight_list ), float( significant_count[ pair ] ) / len( weight_list ) )

if __name__ == "__main__":
    main( )
