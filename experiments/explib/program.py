import os
import subprocess
import sys

from math import sqrt
from . import power
from scipy.stats import chi2

class BayesicMethodFirstStep:
    def __init__(self, name, path, alpha = 10.0, beta = 10.0, num_tests = None, num_single = None, p = None):
        self.name = name
        self.path = path
        self.alpha = alpha
        self.beta = beta
        self.num_tests = num_tests
        self.num_single = num_single
        self.p = p

    def run(self, data_prefix, params, output_file, include_covariates = False):
        num_tests = params.num_tests
        if self.num_tests:
            num_tests = self.num_tests

        cmd = [ self.path,
                "-m", "bayes",
                "-n", str( num_tests ),
                "-a", str( self.alpha ),
                "-b", str( self.beta ),
                data_prefix + ".pair",
                data_prefix ]

        if self.num_single:
            cmd.extend( [ "-s", str( self.num_single ) ] )

        if self.p and self.p > 0.0:
            cmd.extend( [ "-t", str( self.p ) ] )

        if include_covariates:
            cmd.extend( [ "-c", data_prefix + ".cov" ] )

        print " ".join( cmd )
        subprocess.call( cmd, stdout = output_file )

    def num_significant(self, output_path, params, include = None):
        return power.compute_from_file( output_path, 2, params.threshold, params.num_tests, False, include )

    def get_ranks(self, output_path):
        return power.get_ranks( output_path, 2, False )

##
# Wrapper for running Bayesic.
#
class BayesicMethod:
    def __init__(self, name, path, alpha = 10.0, beta = 10.0):
        self.name = name
        self.path = path
        self.alpha = alpha
        self.beta = beta

    def run(self, daa_prefix, params, output_file, include_covariates = False):
        tmp_file = open( output_file.name + ".tmp", "w+" )

        cmd = [ self.path,
                "-m", "bayes",
                "-n", str( params.num_tests ),
                "-a", str( self.alpha ),
                "-b", str( self.beta ),
                data_prefix + ".pair",
                data_prefix ]

        if include_covariates:
            cmd.extend( [ "-c", data_prefix + ".cov" ] )

        print " ".join( cmd )
        subprocess.call( cmd, stdout = tmp_file )

       # # Find significant pairs
        tmp_pair_path = output_file.name + ".tmp.pair"
        tmp_pair = open( tmp_pair_path, "w" )
        tmp_file.seek( 0 )

        for line in tmp_file:
            columns = line.strip( ).split( )
            try:
                posterior = float( columns[ 2 ] )
                if posterior >= 0.0:
                    tmp_pair.write( "{0} {1}\n".format( columns[ 0 ], columns[ 1 ] ) )
            except:
                continue

        tmp_pair.close( )
        tmp_file.close( )

        # Run bayesic on this pair
        cmd = [ self.path,
                "-m", "bayes-fine",
                "-n", str( params.num_tests ),
                "-a", str( self.alpha ),
                "-b", str( self.beta ),
                tmp_pair_path,
                data_prefix ]

        if include_covariates:
            cmd.extend( [ "-c", data_prefix + ".cov" ] )

        print " ".join( cmd )
        subprocess.call( cmd, stdout = output_file )

    def num_significant(self, output_path, params, include = None):
        return power.compute_from_file( output_path, 2, params.threshold, params.num_tests, False, include )

    def get_ranks(self, output_path):
        return power.get_ranks( output_path, 2, False )

    
##
# Wrapper for running log-linear method.
#
class LogLinearMethod:
    def __init__(self, name, path):
        self.name = name
        self.path = path

    def run(self, data_prefix, params, output_file, include_covariates = False):
        cmd = [ self.path,
                "-m", "loglinear",
                "-n", str( params.num_tests ),
                data_prefix + ".pair",
                data_prefix ]

        if include_covariates:
            cmd.extend( [ "-c", data_prefix + ".cov" ] )

        print " ".join( cmd )
        subprocess.call( cmd, stdout = output_file )

    def num_significant(self, output_path, params, include = None):
        return power.compute_from_file( output_path, 2, params.threshold, params.num_tests, True, include )
    
    def get_ranks(self, output_path):
        return power.get_ranks( output_path, 2, True )

##
# Wrapper for running stepwise log-linear method.
#
class ScaleInvarianceMethod:
    def __init__(self, name, path, pvalue_column = 7):
        self.name = name
        self.path = path
        self.pvalue_column = pvalue_column

    def run(self, data_prefix, params, output_file, include_covariates = False):
       
        cmd =[ "python",
               self.path,
               "--only-glm",
               "--alpha", str( params.threshold ),
               data_prefix + ".pair",
               data_prefix
               ]

        num_tests = [ params.num_tests, 0, 0, 0 ]
        cmd.append( "--num-tests" )
        cmd.extend( map( str, num_tests ) )

        if include_covariates:
            cmd.extend( [ "-c", data_prefix + ".cov" ] )

        print " ".join( cmd )
        subprocess.call( cmd, stdout = output_file )

    def num_significant(self, output_path, params, include = None):
        return power.compute_from_file( output_path, self.pvalue_column, params.threshold, 1, True, include, correction = "B" )
    
    def get_ranks(self, output_path):
        return power.get_ranks_stepwise( output_path, [2,3,4,9] )

##
# Wrapper for running stepwise log-linear method.
#
class ClosedMethod:
    def __init__(self, name, path, correction_tool_path, pvalue_column = 7, weight = [ 0.25, 0.25, 0.25, 0.25 ] ):
        self.name = name
        self.path = path
        self.correction_tool_path = correction_tool_path
        self.pvalue_column = pvalue_column
        self.weight = weight

    def run(self, data_prefix, params, output_file, include_covariates = False):
        step1_file = open( output_file.name + "_step1", "w+" )
        
        cmd = [ self.path,
                "-m", "stepwise",
                "-n", str( params.num_tests ),
                data_prefix + ".pair",
                data_prefix ]

        if include_covariates:
            cmd.extend( [ "-c", data_prefix + ".cov" ] )

        print " ".join( cmd )
        subprocess.call( cmd, stdout = step1_file )
        step1_file.close( )

        num_tests = [ n for n in params.closed_num_tests ]
        if params.closed_num_tests[ 0 ] == 0 and params.num_tests > 1:
            num_tests[ 0 ] = params.num_tests
       
        cmd =[ "python",
               self.correction_tool_path,
               "--alpha", str( params.threshold ),
               step1_file.name,
               data_prefix
               ]

        cmd.append( "--weight" )
        cmd.extend( map( str, self.weight ) )

        cmd.append( "--num-tests" )
        cmd.extend( map( str, num_tests ) )

        if include_covariates:
            cmd.extend( [ "-c", data_prefix + ".cov" ] )

        print " ".join( cmd )
        subprocess.call( cmd, stdout = output_file )

    def num_significant(self, output_path, params, include = None):
        return power.compute_from_file( output_path, self.pvalue_column, params.threshold, 1, True, include, correction = "B" )
    
    def get_ranks(self, output_path):
        return power.get_ranks_stepwise( output_path, [2,3,4,9] )

##
# Wrapper for running stepwise log-linear method.
#
class ClosedPartMethod:
    def __init__(self, name, path, pvalue_column = 7):
        self.name = name
        self.path = path
        self.pvalue_column = pvalue_column

    def run(self, data_prefix, params, output_file, include_covariates = False): 
        cmd = [ self.path,
                "-m", "stepwise",
                "-n", str( params.num_tests ),
                data_prefix + ".pair",
                data_prefix ]

        print " ".join( cmd )
        subprocess.call( cmd, stdout = output_file )

    def num_significant(self, output_path, params, include = None):
        return power.compute_from_file( output_path, self.pvalue_column, params.threshold, 1, True, include, correction = "B" )
    
    def get_ranks(self, output_path):
        return power.get_ranks_stepwise( output_path, [2,3,4,9] )

##
# Wrapper for running logistic method.
#
class LogisticMethod:
    def __init__(self, name, path):
        self.name = name
        self.path = path

    def run(self, data_prefix, params, output_file, include_covariates = False):
        cmd = [ self.path,
                "-m", "glm",
                "-l", "logistic",
                "-f", "additive",
                "-n", str( params.num_tests ),
                data_prefix + ".pair",
                data_prefix ]

        if include_covariates:
            cmd.extend( [ "-c", data_prefix + ".cov" ] )
        
        print " ".join( cmd )
        subprocess.call( cmd, stdout = output_file )

    def num_significant(self, output_path, params, include = None):
        return power.compute_from_file( output_path, 4, params.threshold, params.num_tests, True, include )
    
    def get_ranks(self, output_path):
        return power.get_ranks( output_path, 4, True )

##
# Wrapper for running logistic method.
#
class GLMFactorMethod:
    def __init__(self, name, path, method = "log-complement"):
        self.name = name
        self.path = path
        self.method = method

    def run(self, data_prefix, params, output_file, include_covariates = False):
        cmd = [ self.path,
                "-m", "glm",
                "-l", self.method,
                "-n", str( params.num_tests ),
                data_prefix + ".pair",
                data_prefix ]

        if include_covariates:
            cmd.extend( [ "-c", data_prefix + ".cov" ] )
        
        print " ".join( cmd )
        subprocess.call( cmd, stdout = output_file )

    def num_significant(self, output_path, params, include = None):
        return power.compute_from_file( output_path, 3, params.threshold, params.num_tests, True, include, correction = "B" )
    
    def get_ranks(self, output_path):
        return power.get_ranks( output_path, 3, True )

##
# Wrapper for running case only method.
#
class CaseOnlyMethod:
    def __init__(self, name, path):
        self.name = name
        self.path = path

    def run(self, data_prefix, params, output_file, include_covariates = False):
        cmd = [ self.path,
                "-m", "caseonly",
                "-n", str( params.num_tests ),
                data_prefix + ".pair",
                data_prefix ]

        if include_covariates:
            cmd.extend( [ "-c", data_prefix + ".cov" ] )
        
        print " ".join( cmd )
        subprocess.call( cmd, stdout = output_file )

    def num_significant(self, output_path, params, include = None):
        return power.compute_from_file( output_path, 3, params.threshold, params.num_tests, True, include )
    
    def get_ranks(self, output_path):
        return power.get_ranks( output_path, 3, True )

class StepwiseRegression:
    def __init__(self, name, path):
        self.name = name
        self.path = path
        self.num_single_significant = 0

    def run(self, data_prefix, params, output_file, include_covariates = False):
        cmd = [ "plink2", "--bfile", data_prefix, "--model", "--out", output_file.name ]
        subprocess.check_call( cmd, stdout = open( os.devnull, "w" ) )

        significant = set( )
        with open( output_file.name + ".model", "r" ) as assoc_file:
            next( assoc_file )

            for line in assoc_file:
                column = line.strip( ).split( )
                if column[ 4 ] != "GENO":
                    continue

                try:
                    pvalue = float( column[ 9 ] )
                    if pvalue < 0.10:
                        significant.add( column[ 1 ] )
                except:
                    continue

        cmd = [ self.path,
                "-m", "glm",
                "-f", "factor",
                "-l", "logistic",
                data_prefix + ".pair",
                data_prefix ]

        if include_covariates:
            cmd.extend( [ "-c", data_prefix + ".cov" ] )
        
        print " ".join( cmd )
        subprocess.call( cmd, stdout = output_file )

        output_file.seek( 0 )
        lines = output_file.readlines( )
        output_file.seek( 0 )
        output_file.truncate( 0 )

        num_single = max( int( 0.5 + sqrt( 2 * params.num_tests + 0.25 ) ), 1 )
        expected_num_single = int( 0.10 * num_single + 0.5 )
        self.num_single_significant = expected_num_single * ( expected_num_single - 1 ) / 2
        for line in lines:
            column = line.strip( ).split( )

            if not column[ 0 ] in significant or not column[ 1 ] in significant:
                column[ 3 ] = "1.0"
            
            output_file.write( "\t".join( column ) + "\n" )

    def num_significant(self, output_path, params, include = None):
        return power.compute_from_file( output_path, 3, params.threshold, max( self.num_single_significant, 1 ), True, include )
     
    def get_ranks(self, output_path):
        return power.get_ranks( output_path, 3, True )

#g_methods = [ ClosedMethod( "ScaleInv equal", "bayesic", "../tools/closed_correction/closed_correction.py" ),
#              ClosedMethod( "ScaleInv stoh", "bayesic", "../tools/closed_correction/closed_correction.py", weight = [ 0.8, 0.09, 0.09, 0.02 ] ),
#              ClosedMethod( "ScaleInv htos", "bayesic", "../tools/closed_correction/closed_correction.py" , weight = [ 0.02, 0.09, 0.09, 0.8 ] ),
#              ClosedMethod( "ScaleInv estoh", "bayesic", "../tools/closed_correction/closed_correction.py", weight = [ 0.799, 0.1, 0.1, 0.001 ] ),
#              ClosedMethod( "ScaleInv htoes", "bayesic", "../tools/closed_correction/closed_correction.py" , weight = [ 0.001, 0.1, 0.1, 0.799 ] )
#              ]

g_methods = [ LogLinearMethod( "Log-linear", "bayesic" ),
             ClosedMethod( "Stage-wise (static)", "bayesic", "../tools/closed_correction/closed_correction.py" ),
             ScaleInvarianceMethod( "ScaleInv (bonferroni)", "../tools/closed_correction/closed_correction.py" ),
             StepwiseRegression( "Marginal first", "bayesic" ) ]

#g_methods = [ ClosedPartMethod( "Logistic H1", "bayesic", 2 ),
#              GLMFactorMethod( "Logistic H4", "bayesic", "logistic" ) ]

#g_methods = [ StepwiseRegression( "Stepwise", "bayesic" ), ClosedMethod( "ScaleInv (closed)", "bayesic", "../tools/closed_correction/closed_correction.py" ) ]

##
# Returns a list of methods that can run plink files.
#
# @return A list of methods that can run plink files.
#
def get_methods( ):
    global g_methods
    return g_methods

##
# Runs all defined methods on the given plink file.
#
# @param params General experiment parameters.
# @param plink_path Prefix path to the plink, .pair and .cov files.
# @param method_handler Object that handles method output.
# @param include_covariates Determines whether to include covariates or not.
#
def run_methods(params, plink_path, method_handler, include_covariates = False):
    method_list = get_methods( )
    for method in method_list:
        method_output_file = method_handler.get_output_file( method.name )
        method.run( plink_path, params, method_output_file, include_covariates )

def get_confidence_interval(total, num_significant):
    if total > 0.0:
        # Approximate 95% Wilson score interval
        power = ( num_significant + 2.0 ) / ( total + 4.0 )
        err = 1.96 * sqrt( ( 1.0 / total ) * power * ( 1 - power ) )

        lower = max( power - err, 0.0 )
        upper = min( power + err, 1.0 )

        return ( power, lower, upper )
    else:
        return ( 0.0, 0.0, 0.0 )

##
# Calculates power for all the output files found
# in the method handler, and returns a dict mapping method
# names to a tuple that contains power, lower confidence limit,
# and upper confidence limit.
#
# @param param General experiment parameters.
# @param method_handler Object that handles method output.
# @param include Only include the given pairs in the power calculation.
#
# @return A dict containing power estimates for each method.
#
def calculate_power(params, method_handler, include = None):
    method_power = dict( )
    method_list = get_methods( )
    for method in method_list:
        method_output_file = method_handler.get_output_file( method.name )
        num_significant = method.num_significant( method_output_file, params, include )

        method_power[ method.name ] = get_confidence_interval( params.num_pairs, num_significant )

    return method_power

##
# Calculates the number of significant pairs.
#
def num_significant(params, method_handler, include = None):
    method_num_significant = dict( )
    for method in get_methods( ):
        method_output_file = method_handler.get_output_file( method.name )
        method_num_significant[ method.name ] = method.num_significant( method_output_file, params, include )

    return method_num_significant

##
# Calculates a ranking of all SNP pairs.
#
def get_ranks(params, method_handler):
    method_rank = dict( )
    method_list = get_methods( )
    for method in method_list:
        method_output_file = method_handler.get_output_file( method.name )
        method_rank[ method.name ] = method.get_ranks( method_output_file )

    return method_rank

