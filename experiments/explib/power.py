from math import sqrt

##
# Estimates power and confidence interval from a file
# that contains ranked entities.
#
# @param csv_file File to read from.
# @param column Column that contains rank.
# @param operator Floating point operator to use when comparing
#                 the rank to the threshold.
# @param threshold Threshold to compare the rank with.
#
# @return A tuple containing estimated power, lower confidence
#         limit, and upper confidence limit.
#
def compute_from_file(csv_file, column, operator, threshold):
    num_significant = 0
    total = 0
    for line in csv_file:
        column_list = line.strip( ).split( )

        try:
            value = float( column_list[ column ] )
        except:
            continue

        if value < 0.0 or value > 1.0:
            continue

        if operator( value, threshold ):
            num_significant += 1

        total += 1

    if total != 0:
        # Approximate 95% Wilson score interval
        power = ( num_significant + 2.0 ) / ( total + 4.0 )
        err = 1.96 * sqrt( ( 1.0 / total ) * power * ( 1 - power ) )

        lower = max( power - err, 0.0 )
        upper = min( power + err, 1.0 )

        return ( power, lower, upper )
    else:
        return ( 0.0, 0.0, 0.0 )
