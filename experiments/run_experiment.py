import argparse
import json
import os
import tempfile
import subprocess

from explib import program, plot, plinkdata, util, expfile

##
# Contains general parameters for the experiment.
#
class ExperimentParams:
    ##
    # Constructor.
    #
    def __init__(self, maf = [ 0.0, 0.0 ], sample_size = [ 0, 0 ], num_pairs = 0, num_tests = 0):
        ##
        # Minor allele frequency.
        #
        self.maf = maf

        ##
        # List of sample sizes for controls and cases.
        #
        self.sample_size = sample_size

        ##
        # Number of pairs to generate when estimating power.
        #
        self.num_pairs = num_pairs

        ##
        # Number of statistical tests to adjust for.
        #
        self.num_tests = num_tests

        ##
        # The false positive rate.
        #
        self.threshold = 0.05

def updated_params(params, json_object):
    new_params = ExperimentParams( params.maf, params.sample_size, params.num_pairs, params.num_tests )
    
    new_params.maf = json_object.get( "maf", params.maf )
    new_params.sample_size = json_object.get( "sample_size", params.sample_size )
    new_params.num_pairs = json_object.get( "num_pairs", params.num_pairs )
    new_params.num_tests = json_object.get( "num_tests", params.num_tests )
    new_params.threshold = json_object.get( "threshold", params.threshold )

    return new_params

##
# Computes the heritability V(P|G) / V(P).
#
# @param penetrance The probability of disease for each genotype as a vector.
# @param maf The minor allele frequency for both locus.
#
# @return The heritability.
#
def heritability(penetrance, maf):
    maf1 = maf[ 0 ]
    maf2 = maf[ 1 ]

    p = [ ( 1 - maf[ 0 ] )**2, 2 * maf[ 0 ] * ( 1 - maf[ 0 ] ), ( maf[ 0 ] )**2 ]
    q = [ ( 1 - maf[ 1 ] )**2, 2 * maf[ 1 ] * ( 1 - maf[ 1 ] ), ( maf[ 1 ] )**2 ]

    joint_maf =  [ p[ 0 ] * q[ 0 ], p[ 0 ] * q[ 1 ], p[ 0 ] * q[ 2 ],
                   p[ 1 ] * q[ 0 ], p[ 1 ] * q[ 1 ], p[ 1 ] * q[ 2 ],
                   p[ 2 ] * q[ 0 ], p[ 2 ] * q[ 1 ], p[ 2 ] * q[ 2 ] ]

    pop_p = sum( p * m for p, m in zip( penetrance, maf ) )

    h = 0.0
    for i in range( 3 ):
        for j in range( 3 ):
            cell = 3 * i + j
            h += ( penetrance[ cell ] - pop_p )**2 * joint_maf[ cell ]

    return h / ( pop_p * ( 1 - pop_p ) )

##
# Main class for simulating data for a set of parameter
# combinations.
#
# An experiment is defined by a json file of the following
# structure.
#
# {
#   "experiments" =
#   [
#       {
#           "type" : "covariate or genetic",
#           "name" : "anything",
#           "title" : "Plot title",
#           "xlabel" : "Label on x-axis",
#           "ylabel" : "Label on y-axis",
#           "models" :
#            [
#               {
#                   "params" : [ 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1 ],
#                   "snpbeta" : [ 0.1, 0.1, 0.1 ],
#                   "covariates" : [ "age,0.8,normal,65,5" ]
#               }
#            ]
#       }
#   ]
# }
#
# Since the parameter combinations, data generation and power 
# estimation can be very different in nature for different
# experiments, the type of experiment is mapped to a Python class
# that handles that particular type of experiment.
#
class GeneticExperiment:
    ##
    # Constructor.
    #
    # @param params General parameters.
    # @param output_dir The output directory.
    #
    def __init__(self, params, output_dir):
        ##
        # General experiment parameters.
        #
        self.params = params

        ##
        # Object that handles paths.
        #
        self.path_handler = expfile.PathHandler( output_dir, [ "plots", "tmp", "method", "power" ] )

    ##
    # Initializes the output directory, and sets up sub directories.
    #
    def init_output_dir(self):
        for subdir in self.path_handler.get_subdirs( ):
            if not os.path.exists( subdir ):
                os.makedirs( subdir )

    ##
    # Create the different type of strategies for different kind
    # of experiments.
    #
    # @param method_handler MethodHandler object.
    #
    def create_strategies(self, method_handler):
        return { "genetic" : NoCovariateStrategy( method_handler ),
                 "covariate" : CovariateStrategy( method_handler ),
                 "sufficient" : SufficientCovariateStrategy( method_handler ),
                 "mixed" : MixedStrategy( method_handler ) }

    ##
    # Runs the experiments defined by the given json array of
    # experiment objects.
    #
    # @param experiment_list Json array of experiment objects.
    #
    def run(self, experiment_list):
        self.init_output_dir( )

        plink_path = self.path_handler.get_plink_path( )
        method_handler = expfile.MethodHandler( self.path_handler.get_method_base_path( ) )
        strategies = self.create_strategies( method_handler )

        for experiment in experiment_list:
            strategy = strategies.get( experiment.get( "type", "" ), None )
            if not strategy:
                print( "run_experiment.py: warning: Experiment {0} has no type, continuing.".format( experiment[ 'name' ] ) )

            power_path = self.path_handler.get_power_path( experiment[ 'name' ] )
            with open( power_path, "w" ) as power_file:
                strategy.calculate_power( self.params, experiment, plink_path, power_file )
                
            plot_path = self.path_handler.get_plot_path( experiment[ 'name' ] )
            strategy.plot_power( self.params, experiment, power_path, plot_path )

##
# Abstract class for handling data generation, power calculation
# and plotting for a particular kind of experiment.
#
class ExperimentStrategy:
    pass

##
# Handles experiments that lacks covariates.
#
class NoCovariateStrategy(ExperimentStrategy):
    def __init__(self, method_handler):
        self.method_handler = method_handler

    def calculate_power(self, params, experiment, plink_path, power_file):
        ld = experiment.get( "ld", 0.0 )
        for model_id, model in enumerate( experiment[ 'models' ] ):
            params = updated_params( params, model )

            plinkdata.generate_data( params,
                                     model[ 'params' ],
                                     ld, plink_path )

            self.method_handler.start_experiment( experiment[ 'name' ], model_id )
            program.run_methods( params, plink_path, self.method_handler )
            self.method_handler.reset_files( )

            xvalue = model.get( 'xvalue', None )
            if not xvalue:
                xvalue = str( heritability( model[ 'params' ], params.maf ) )

            method_power = program.calculate_power( params, self.method_handler )
            for method_name, power_data in method_power.iteritems( ):
                power, lower, upper = power_data
                line = "{0}\t{1}\t{2}\t{3}\t{4}\n".format( method_name, xvalue, power, lower, upper )
                power_file.write( line )

    def plot_power(self, params, experiment, power_path, plot_path):
        plot.plot_power( power_path,
                         experiment[ 'title' ],
                         experiment[ 'xlabel' ],
                         experiment[ 'ylabel' ],
                         params,
                         plot_path )


class MixedStrategy(ExperimentStrategy):
    def __init__(self, method_handler):
        self.method_handler = method_handler
    
    def calculate_power(self, params, experiment, plink_path, power_file):
        all_model_params = list( )
        for model_id, model in enumerate( experiment[ 'models' ] ):
            params = updated_params( params, model )

            all_model_params.append( ( params.num_pairs, model[ 'is_case' ], model[ 'params' ] ) )

        plinkdata.generate_mixed_data( params, all_model_params, plink_path )
        
        self.method_handler.start_experiment( experiment[ 'name' ], 0 )
        program.run_methods( params, plink_path, self.method_handler )
        self.method_handler.reset_files( )

        is_case_map = dict( )
        with open( plink_path + ".case", "r" ) as case_file:
            for line in case_file:
                snp1, snp2, is_case = line.strip( ).split( )
                is_case_map[ ( snp1, snp2 ) ] = is_case

        method_rank = program.get_ranks( params, self.method_handler )
        for method_name, snp_rank_data in method_rank.iteritems( ):
            for snp1, snp2, rank in snp_rank_data:
                is_case = is_case_map.get( ( snp1, snp2 ), -1 )
                if is_case == -1:
                    continue

                line = "{0}\t{1}\t\"{2}\"\t{3}\t{4}\n".format( snp1, snp2, method_name, rank, is_case )
                power_file.write( line )

    def plot_power(self, params, experiment, power_path, plot_path):
        plot.plot_roc( power_path, experiment[ 'title' ], plot_path )

##
# Handles experiments that contains covariates
# generated from the logistic model.
#
class CovariateStrategy(ExperimentStrategy):
    def __init__(self, method_handler):
        self.method_handler = method_handler

    ##
    # Calculate power with or without covariates.
    #
    def run_experiment(self, params, experiment, plink_path, power_file, with_cov):
        postfix = "with_cov"
        if not with_cov:
            postfix = "without_cov"

        for model_id, model in enumerate( experiment[ 'models' ] ): 
            covariates = model[ 'covariates' ]
            if not with_cov:
                covariates = None


            params = updated_params( params, model )
            self.method_handler.start_experiment( experiment[ 'name' ], model_id, postfix )
            for i in range( params.num_pairs ):
                plinkdata.generate_cov_data( params,
                                             model[ 'snpbeta' ],
                                             covariates,
                                             plink_path )

                program.run_methods( params, plink_path, self.method_handler, with_cov )
            
            xvalue = model.get( 'xvalue', None )
            if not xvalue:
                xvalue = model_id

            self.method_handler.reset_files( )
            method_power = program.calculate_power( params, self.method_handler )
            for method_name, power_data in method_power.iteritems( ):
                power, lower, upper = power_data
                line = "{0}\t{1}\t{2}\t{3}\t{4}\t{5}\n".format( method_name, xvalue, power, lower, upper, int( with_cov ) )
                power_file.write( line )


    def calculate_power(self, params, experiment, plink_path, power_file):
        self.run_experiment( params, experiment, plink_path, power_file, True )
        self.run_experiment( params, experiment, plink_path, power_file, False )
    
    def plot_power(self, params, experiment, power_path, plot_path):
        plot.plot_cov_power( power_path,
                         experiment[ 'title' ],
                         experiment[ 'xlabel' ],
                         experiment[ 'ylabel' ],
                         params,
                         plot_path )

class SufficientCovariateStrategy(ExperimentStrategy):
    def __init__(self, method_handler):
        self.method_handler = method_handler

    ##
    # Calculate power with or without covariates.
    #
    def run_experiment(self, params, experiment, plink_path, power_file, with_cov):
        postfix = "with_cov"
        if not with_cov:
            postfix = "without_cov"

        for model_id, model in enumerate( experiment[ 'models' ] ): 
            covariates = model[ 'covariates' ]
            if not with_cov:
                covariates = None


            params = updated_params( params, model )
            self.method_handler.start_experiment( experiment[ 'name' ], model_id, postfix )
            for i in range( params.num_pairs ):
                plinkdata.generate_sufficient_cov_data( params,
                                             model[ 'params' ],
                                             covariates,
                                             plink_path )

                program.run_methods( params, plink_path, self.method_handler, with_cov )
            
            xvalue = model.get( 'xvalue', None )
            if not xvalue:
                xvalue = model_id

            self.method_handler.reset_files( )
            method_power = program.calculate_power( params, self.method_handler )
            for method_name, power_data in method_power.iteritems( ):
                power, lower, upper = power_data
                line = "{0}\t{1}\t{2}\t{3}\t{4}\t{5}\n".format( method_name, xvalue, power, lower, upper, int( with_cov ) )
                power_file.write( line )


    def calculate_power(self, params, experiment, plink_path, power_file):
        self.run_experiment( params, experiment, plink_path, power_file, True )
        self.run_experiment( params, experiment, plink_path, power_file, False )
    
    def plot_power(self, params, experiment, power_path, plot_path):
        plot.plot_cov_power( power_path,
                         experiment[ 'title' ],
                         experiment[ 'xlabel' ],
                         experiment[ 'ylabel' ],
                         params,
                         plot_path )



##
# Argparse json file loader.
#
def json_file(json_path):
    return json.load( open( json_path ) )

if __name__ == "__main__":
    
    parser = argparse.ArgumentParser( description = 'Runs all methods on the given models.' )
    parser.add_argument( 'experiment_file', type=json_file, help='JSON file describing the experiments to run.' )
    parser.add_argument( 'output_dir', type=str, help='Path where all output and temporary files will be placed.' )
    parser.add_argument( '-m', type=float, nargs=2, help="Minor allele frequency of SNPs.", default = [ 0.2, 0.2 ] )
    parser.add_argument( '-n', type=int, nargs=2, help="Number of controls and cases.", default = [ 1893, 1525 ] )
    parser.add_argument( '-s', type=int, help="Number of model instances to use when estimating power", default = 200 )
    parser.add_argument( '-i', type=int, help="Number of interactions to adjust for", default = 4000000 )

    args = parser.parse_args( )

    params = ExperimentParams( args.m, args.n, args.s, args.i )
    experiment = GeneticExperiment( params, args.output_dir )
    experiment.run( args.experiment_file[ 'experiments' ] )
