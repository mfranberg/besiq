from math import sqrt, log

def compute_from_file(csv_file, column, threshold, num_tests, is_pvalue = True, include = None):
    probs = read_from_file( csv_file, column, threshold, include )
    
    total = len( probs )
    num_significant = 0.0
    if is_pvalue:
        num_significant = compute_from_pvalues( probs, threshold, num_tests )
    else:
        num_significant = compute_from_posterior( probs, threshold )

    return get_confidence_interval( total, num_significant )

def get_confidence_interval(total, num_significant):
    if total != 0:
        # Approximate 95% Wilson score interval
        power = ( num_significant + 2.0 ) / ( total + 4.0 )
        err = 1.96 * sqrt( ( 1.0 / total ) * power * ( 1 - power ) )

        lower = max( power - err, 0.0 )
        upper = min( power + err, 1.0 )

        return ( power, lower, upper )
    else:
        return ( 0.0, 0.0, 0.0 )

def read_from_file(csv_file, column, include = None):
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

def read_from_file_stepwise(csv_file, column_indices, pairs = [], include = None):
    csv_file.seek( 0 )
    value_list = list( )
    for i in range( len( column_indices ) ):
        value_list.append( [ ] )

    for line in csv_file:
        column_list = line.strip( ).split( )
       
        if include and not ( column_list[ 0 ], column_list[ 1 ] ) in include:
            continue

        no_exception = True
        line_values = [ ]
        for column in column_indices:
            value = 0.0
            try:
                value =  float( column_list[ column ] )
            except:
                no_exception = False
                break

            if value >= 0.0 and value <= 1.0:
                line_values.append( value )
        
        if no_exception:
            for i, value in enumerate( line_values ):
                value_list[ i ].append( value )

            pairs.append( ( column_list[ 0 ], column_list[ 1 ] ) )

    return value_list

def compute_from_posterior(posterior, threshold):
    cur_fdr = 0.0
    significant = 0.0
    for i, p in enumerate( sorted( posterior, reverse = True ) ):
        prev_fdr = cur_fdr
        cur_fdr = ( prev_fdr * i + ( 1.0 -  p ) ) / ( i + 1 )
        if cur_fdr > threshold:
            break
        else:
            significant += 1
   
    return significant

def compute_from_pvalues(pvalues, threshold, num_tests):
    significant = 0.0
    for i, p in enumerate( sorted( pvalues ) ):
        if p > min( ( (i + 1.0 ) / num_tests) * threshold, threshold ):
            break
        else:
            significant += 1

    return significant

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
