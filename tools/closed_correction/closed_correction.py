import sys
import os
import subprocess
import argparse
from collections import defaultdict

##
# Perform correction for the log-linear models, storing the number of tests
# in each level, and a path containing the snp pairs that have survived correction
# so far.
#
# @param result_file The file that contains all results.
# @param alpha The significance threshold (should not be corrected for multiple tests).
# @param step1_tests The number of tests in the first step in case all were not included
#                    in the results file. Only used if not None.
# @param step2_tests The number of tests in the second step, used in simulated data to
#                    increase the multiple testing correction. Only used if not None.
#
# @return A tuple containing, a list of the number of tests for each level, path to a
#         a file that contains the remaining results.
#
def correct_non_regression_models(result_file, alpha, step1_tests, step2_tests):
    NUM_LEVELS = 3
    num_tests = [ step1_tests ]
    if not step1_tests:
        num_tests = [ sum( 1 for line in open( result_file, "r" ) if line[ :3 ] != "snp1" ) ]

    previous = result_file

    for k in range( NUM_LEVELS ):
        if num_tests[ -1 ] <= 0:
            break

        suffix = ".level{0}".format( k )
        cmd = [ "filter", str( 2 + k ), "le", str( alpha / num_tests[ -1 ] ), previous, ">", result_file + suffix ]
        subprocess.check_call( " ".join( cmd ), shell = True )

        if k == 0 and step2_tests:
            num_tests.append( step2_tests )
        else:
            num_tests.append( sum( 1 for line in open( result_file + suffix, "r" ) ) )

        previous = result_file + suffix

    return ( num_tests, previous )

##
# Reads the pairs and their p-value from a GLM output file.
#
# @param result_file A path to a file that contains the GLM p-values.
#
# @return A list of tuples containing, the first snp, the second snp, and
#         the p-value. The p-value is -1 if it could not be parsed.
#
def read_pairs(result_file):
    pairs = [ ]
    for line in open( result_file, "r" ):
        column = line.strip( ).split( )
        snp1 = column[ 0 ]
        snp2 = column[ 1 ]

        try:
            pvalue = float( column[ 3 ] )
            pairs.append( ( snp1, snp2, pvalue ) )
        except:
            pairs.append( ( snp1, snp2, -1.0 ) )
            continue
    
    return pairs

##
# Runs each regression model in the list and returns a corresonding
# list of pairs and p-values for each.
#
# @param method_list A list of glm to run.
# @param pairs_file A path to a file that contains the pairs to run with bayesic.
# @param plink_file A path to the plink-file to run with bayesic.
# @param result_file The output prefix of the result files.
#
# @return A list of lists, each containing the output form read_pairs.
#
def run_regression_models(method_list, pairs_file, plink_file, result_file):
    method_results = list( )
    for method in method_list:
        regression_result_file = result_file + "." + method
        cmd = [ "bayesic", "-m", "glm", "-l", method, pairs_file, plink_file, ">", regression_result_file ]
        subprocess.check_call( " ".join( cmd ), shell = True )
        method_results.append( read_pairs( regression_result_file ) )

    return method_results

##
# Parses the p-values from the saturated models and stores them as a mapping
# from pair to a list of p-values.
#
# @param last_result_file Path to the file that contains p-valeus for the saturated models.
#
# @return A mapping from pair to a list of p-values.
#
def read_simple_pvalues(last_result_file):
    pair_pvalue = dict( )
    with open( last_result_file, "r" ) as lrf:
        for line in lrf:
            column = line.strip( ).split( )

            pair = ( column[ 0 ], column[ 1 ] )
            
            try:
                pair_pvalue[ pair ] = list( map( float, column[ 2:5 ] ) )
            except:
                continue

    return pair_pvalue

##
# Calculates the adjusted p-values for each method.
#
# @param num_tests The number of tests for each stage.
# @param simple_pvalues The p-values for each saturated model.
# @param method_results The p-values for each glm.
#
# @return A mapping from each pair to a list of adjusted p-values, one
#         for each method.
#
def calculate_adjusted_results(num_tests, simple_pvalues, method_results):
    adjusted_simple_pvalues = dict( )
    for pair, pvalues in simple_pvalues.iteritems( ):
        adjusted_pvalues = [ min( n * p, 1 ) for n, p in zip( num_tests[:-1], pvalues ) ]
        adjusted_simple_pvalues[ pair ] = max( adjusted_pvalues )

    adjusted_method_results = defaultdict( list )
    for pairs in method_results:
        for snp1, snp2, pvalue in pairs:
            if pvalue >= 0.0:
                adjusted_pvalue = max( adjusted_simple_pvalues[ (snp1, snp2) ], num_tests[ -1 ] * pvalue )
                adjusted_method_results[ ( snp1, snp2 ) ].append( adjusted_pvalue ) 
            else:
                adjusted_method_results[ ( snp1, snp2 ) ].append( -1.0 )

    # Add final p-value
    for pair, results in adjusted_method_results.iteritems( ):
        adjusted_method_results[ pair ].append( max( results ) )

    return adjusted_method_results

##
# Print the snp pairs that are significant.
#
# @param method_list The list of methods that was run.
# @param adjusted_method_results The adjusted p-values for each method.
# @param alpha The significance threshold.
#
def print_significant(method_list, adjusted_method_results, alpha):
    header = "snp1 snp2\t" + "\t".join( method_list ) + "\tP-value"
    print header

    for pair, result in adjusted_method_results.iteritems( ):
        is_significant = map( lambda x: x >= 0.0 and x <= alpha, result )
        if any( is_significant ):
            sys.stdout.write( "{0} {1}".format( pair[ 0 ], pair[ 1 ] ) )
            for p in result:
                if p >= 0.0:
                    sys.stdout.write( "\t{0:.4g}".format( p ) )
                else:
                    sys.stdout.write( "\tNA" )

            sys.stdout.write( "\n" )

def main():
    parser = argparse.ArgumentParser( description='Computed adjusted p-values for each scale.' )
    parser.add_argument( 'result_file', type=str, help='The output file from bayesic -m stepwise.' )
    parser.add_argument( 'plink_file', type=str, help='The same plink file that was used to create the result file.' )
    parser.add_argument( '--alpha', type=float, help='The significance level', default = 0.05 )
    parser.add_argument( '--step1-tests', type=int, help='The number of tests to adjust the first step for.', default = None )
    parser.add_argument( '--step2-tests', type=int, help='The number of tests to adjust the first step for.', default = None )
    args = parser.parse_args( )

    num_tests, last_result_file = correct_non_regression_models( args.result_file, args.alpha, args.step1_tests, args.step2_tests )
    if min( num_tests ) <= 0:
        exit( )

    last_pairs = last_result_file + ".pairs"
    cmd = [ "cut", "-f", "1", last_result_file, ">", last_pairs ]
    subprocess.check_call( " ".join( cmd ), shell = True )
    
    method_list = [ "odds-additive", "logistic", "penetrance-additive", "penetrance-multiplicative", "log-complement" ]
    method_results = run_regression_models( method_list, last_pairs, args.plink_file, last_result_file )
    simple_pvalues = read_simple_pvalues( last_result_file )
    adjusted_method_results = calculate_adjusted_results( num_tests, simple_pvalues, method_results )
    print_significant( method_list, adjusted_method_results, args.alpha )

if __name__ == "__main__":
    main( )
