import subprocess

from . import fdr_power

class BayesicMethodFirstStep:
    def __init__(self, name, path, alpha = 10.0, beta = 10.0, num_tests = None, num_single = None):
        self.name = name
        self.path = path
        self.alpha = alpha
        self.beta = beta
        self.num_tests = num_tests
        self.num_single = num_single

    def run(self, data_prefix, num_tests, output_file, include_covariates = False):
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

        if include_covariates:
            cmd.extend( [ "-c", data_prefix + ".cov" ] )

        print " ".join( cmd )
        subprocess.call( cmd, stdout = output_file )

    def compute_power(self, output_path, num_tests, threshold = 0.05):
        return fdr_power.compute_from_file( output_path, 2, threshold, num_tests, False )

    def get_ranks(self, output_path):
        return fdr_power.get_ranks( output_path, 2, False )

##
# Wrapper for running Bayesic.
#
class BayesicMethod:
    def __init__(self, name, path, alpha = 10.0, beta = 10.0):
        self.name = name
        self.path = path
        self.alpha = alpha
        self.beta = beta

    def run(self, daa_prefix, num_tests, output_file, include_covariates = False):
        tmp_file = open( output_file.name + ".tmp", "w+" )

        cmd = [ self.path,
                "-m", "bayes",
                "-n", str( num_tests ),
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
                "-n", str( num_tests ),
                "-a", str( self.alpha ),
                "-b", str( self.beta ),
                tmp_pair_path,
                data_prefix ]

        if include_covariates:
            cmd.extend( [ "-c", data_prefix + ".cov" ] )

        print " ".join( cmd )
        subprocess.call( cmd, stdout = output_file )

    def compute_power(self, output_path, num_tests, threshold = 0.05):
        return fdr_power.compute_from_file( output_path, 2, threshold, num_tests, False )

    def get_ranks(self, output_path):
        return fdr_power.get_ranks( output_path, 2, False )

    
##
# Wrapper for running log-linear method.
#
class LogLinearMethod:
    def __init__(self, name, path):
        self.name = name
        self.path = path

    def run(self, data_prefix, num_tests, output_file, include_covariates = False):
        cmd = [ self.path,
                "-m", "loglinear",
                "-n", str( num_tests ),
                data_prefix + ".pair",
                data_prefix ]

        if include_covariates:
            cmd.extend( [ "-c", data_prefix + ".cov" ] )

        print " ".join( cmd )
        subprocess.call( cmd, stdout = output_file )

    def compute_power(self, output_path, num_tests, threshold = 0.05):
        return fdr_power.compute_from_file( output_path, 2, threshold, num_tests, True )
    
    def get_ranks(self, output_path):
        return fdr_power.get_ranks( output_path, 2, True )

##
# Wrapper for running logistic method.
#
class LogisticMethod:
    def __init__(self, name, path):
        self.name = name
        self.path = path

    def run(self, data_prefix, num_tests, output_file, include_covariates = False):
        cmd = [ self.path,
                "-m", "logistic",
                "-n", str( num_tests ),
                data_prefix + ".pair",
                data_prefix ]

        if include_covariates:
            cmd.extend( [ "-c", data_prefix + ".cov" ] )
        
        print " ".join( cmd )
        subprocess.call( cmd, stdout = output_file )

    def compute_power(self, output_path, num_tests, threshold = 0.05):
        return fdr_power.compute_from_file( output_path, 4, threshold, num_tests, True )
    
    def get_ranks(self, output_path):
        return fdr_power.get_ranks( output_path, 4, True )

##
# Wrapper for running logistic method.
#
class LogisticFactorMethod:
    def __init__(self, name, path):
        self.name = name
        self.path = path

    def run(self, data_prefix, num_tests, output_file, include_covariates = False):
        cmd = [ self.path,
                "-m", "logistic-factor",
                "-n", str( num_tests ),
                data_prefix + ".pair",
                data_prefix ]

        if include_covariates:
            cmd.extend( [ "-c", data_prefix + ".cov" ] )
        
        print " ".join( cmd )
        subprocess.call( cmd, stdout = output_file )

    def compute_power(self, output_path, num_tests, threshold = 0.05):
        return fdr_power.compute_from_file( output_path, 3, threshold, num_tests, True )
    
    def get_ranks(self, output_path):
        return fdr_power.get_ranks( output_path, 3, True )

##
# Returns a list of methods that can run plink files.
#
# @return A list of methods that can run plink files.
#
def get_methods( ):
    return [ BayesicMethod( "Bayesic", "../build/src/bayesic" ),
             LogLinearMethod( "Log-linear", "../build/src/bayesic" ),
             LogisticFactorMethod( "Logistic", "../build/src/bayesic" ) ]

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
        method.run( plink_path, params.num_tests, method_output_file, include_covariates )

##
# Calculates power for all the output files found
# in the method handler, and returns a dict mapping method
# names to a tuple that contains power, lower confidence limit,
# and upper confidence limit.
#
# @param param General experiment parameters.
# @param method_handler Object that handles method output.
#
# @return A dict containing power estimates for each method.
#
def calculate_power(params, method_handler):
    method_power = dict( )
    method_list = get_methods( )
    for method in method_list:
        method_output_file = method_handler.get_output_file( method.name )
        power, lower, upper = method.compute_power( method_output_file, params.num_tests, params.threshold )

        method_power[ method.name ] = ( power, lower, upper )

    return method_power

def get_ranks(params, method_handler):
    method_rank = dict( )
    method_list = get_methods( )
    for method in method_list:
        method_output_file = method_handler.get_output_file( method.name )
        method_rank[ method.name ] = method.get_ranks( method_output_file )

    return method_rank

