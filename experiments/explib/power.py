from math import sqrt, log
from . import fdr_power

def compute_from_file(csv_file, column, threshold, num_tests, is_pvalue = True, include = None):
    probs = fdr_power.read_from_file( csv_file, column, include )
    
    total = len( probs )
    num_significant = 0.0
    if is_pvalue:
        num_significant = compute_from_pvalues( probs, threshold, num_tests )
    else:
        num_significant = compute_from_posterior( probs, threshold )

    return fdr_power.get_confidence_interval( total, num_significant )

def compute_from_file_stepwise(csv_file, column, threshold, num_tests, include = None):
    probs = fdr_power.read_from_file_stepwise( csv_file, column, include = include )

    total = len( probs[ 0 ] )
    num_significant = compute_from_pvalues_stepwise( probs, threshold, num_tests )

    return fdr_power.get_confidence_interval( total, num_significant )

def compute_from_pvalues_stepwise(pvalues, threshold, num_tests):
    corrected_threshold = threshold / num_tests
    indices = range( len( pvalues[ 0 ] ) )
    levels = range( len( pvalues ) )

    for l in levels:
        next_indices = []
        for i in indices:
            if pvalues[ l ][ i ] <= corrected_threshold:
                next_indices.append( i )

        if len( next_indices ) <= 0:
            return 0

        corrected_threshold = threshold / min( len( next_indices ), num_tests )
        indices = next_indices

    return len( next_indices )

def compute_from_posterior(posterior, threshold):
    return len( [ p for p in posterior if p >= threshold ] )

def compute_from_pvalues(pvalues, threshold, num_tests):
    return len( [ p for p in pvalues if p <= threshold / num_tests ] )

def get_ranks(csv_file, column, is_pvalue = True):
    return fdr_power.get_ranks( csv_file, column, is_pvalue )

def get_ranks_stepwise(csv_file, column):
    valid_pairs = [ ]
    probs = fdr_power.read_from_file_stepwise( csv_file, column, pairs = valid_pairs )

    total = len( probs[ 0 ] )
    ranks = [ 0 ] * total
    levels = range( len( probs ) )

    num_tests = [ ]
    indices = list( range( total ) )
    for l in levels:
        if len( indices ) <= 0:
            break

        alpha = 0.05 / len( indices )
        num_tests.append( len( indices ) )

        next_indices = [ ]
        for i in indices:
            if probs[ l ][ i ] <= alpha:
                next_indices.append( i )
            else:
                adjusted_p = 0.0
                for k in range( l + 1 ):
                    adjusted_p = max( adjusted_p, num_tests[ k ] * probs[ k ][ i ] )

                ranks[ i ] = -adjusted_p

        indices = next_indices

    result = [ ]
    for pair, rank in zip( valid_pairs, ranks ):
        result.append( ( pair[ 0 ], pair[ 1 ], rank ) )

    return result
