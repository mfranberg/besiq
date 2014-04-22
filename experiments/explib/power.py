from math import sqrt, log


def read_significance_value_from_file(csv_file, column, include = None):
    csv_file.seek( 0 )
    value_list = list( )
    for line in csv_file:
        column_list = line.strip( ).split( )

        if include and not ( column_list[ 0 ], column_list[ 1 ] ) in include:
            continue

        value = 0.0
        try:
            value =  float( column_list[ column ] )
        except:
            continue

        if value >= 0.0 and value <= 1.0:
            value_list.append( value )

    return value_list

def compute_from_file(csv_file, column, threshold, num_tests, is_pvalue = True, include = None, correction = "HB"):
    probs = read_significance_value_from_file( csv_file, column, include )
    
    total = len( probs )
    num_significant = 0.0
    if is_pvalue:
        if correction == "HB":
            num_significant = compute_from_pvalues_hb( probs, threshold, num_tests )
        elif correction == "B":
            num_significant = compute_from_pvalues_b( probs, threshold, num_tests )
    else:
        num_significant = compute_from_posterior( probs, threshold )

    return ( total, num_significant )

def compute_from_posterior(posterior, threshold):
    return len( [ p for p in posterior if p >= threshold ] )

##
# Determines the number of significant p-value using the
# Holm-Bonferroni method.
#
# @param pvalues List of p-values.
# @param threshold The significant threshold.
# $param num_tests The total number of tests.
#
# @return The number of significant p-values.
#
def compute_from_pvalues_hb(pvalues, threshold, num_tests):
    num_significant = 0
    for i, p in enumerate( sorted( pvalues ) ):
        if p > ( threshold / ( num_tests - i ) ):
            break

        num_significant += 1
    
    return num_significant

##
# Determines the number of significant p-value using the
# Bonferroni method.
#
# @param pvalues List of p-values.
# @param threshold The significant threshold.
# $param num_tests The total number of tests.
#
# @return The number of significant p-values.
#
def compute_from_pvalues_b(pvalues, threshold, num_tests):
    return len( [ p for p in pvalues if p <= ( threshold / num_tests ) ] )

def get_ranks(csv_file, column, is_pvalue = True):
    ranks = list( )
    for line in csv_file:
        column_list = line.strip( ).split( )

        value = 0.0
        try:
            value =  float( column_list[ column ] )
        except:
            continue

        if value >= 0.0 and value <= 1.0:
            if is_pvalue:
                value = -value
            ranks.append( ( column_list[ 0 ], column_list[ 1 ], value ) )

    return ranks

def get_ranks_stepwise(csv_file, column):
    valid_pairs = [ ]
    probs = read_significance_value_from_file_stepwise( csv_file, column, pairs = valid_pairs )

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
