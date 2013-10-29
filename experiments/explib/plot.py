import subprocess
import os

##
# Returns the path of this script.
#
def get_script_path():
    return os.path.join( os.path.realpath( os.path.dirname( __file__ ) ), "external" )

##
# Creates a title that contains all general experiment
# parameters.
#
# @param title The title.
# @param params General experiment parameters.
#
# @return The modified title.
#
def create_title(title, params):
    return "{0}\n(cases = {1}, controls = {2},\nmaf1 = {3}, maf2 = {4}, interactions = {5:g}).".format( title, params.sample_size[ 1 ], params.sample_size[ 0 ], params.maf[ 0 ], params.maf[ 1 ], params.num_tests )

##
# Plots the power of an experiment without covariates.
#
# @param power_path Path to the file containing power estimates.
# @param title Title of the plot.
# @param xlabel Label on the x-axis.
# @param ylabel Label on the y-axis.
# @param params General experiment parameters.
# @param output_path Output path of plot.
#
def plot_power(power_path, title, xlabel, ylabel, params, output_path):
    full_title = create_title( title, params )

    cmd = [ "Rscript", os.path.join( get_script_path( ), "plot_model_power.r" ),
            xlabel, ylabel, full_title, power_path, output_path ]

    subprocess.call( cmd )

##
# Plots the power of an experiment with covariates.
#
# @param power_path Path to the file containing power estimates.
# @param title Title of the plot.
# @param xlabel Label on the x-axis.
# @param ylabel Label on the y-axis.
# @param params General experiment parameters.
# @param output_path Output path of plot.
#
def plot_cov_power(power_path, title, xlabel, ylabel, params, output_path):
    full_title = create_title( title, params )

    cmd = [ "Rscript", os.path.join( get_script_path( ), "plot_covariate_power.r" ),
            xlabel, ylabel, full_title, power_path, output_path ]

    subprocess.call( cmd )

##
# Plots a roc curve showing how well interaction models
# are separated from non-interaction models.
#
# @param roc_path Path to a file containg pairs of snps, their rank,
#                 and whether they are an interaction or not.
# @param title Title of the plot.
# @param output_path The path to the output pdf.
#
def plot_roc(roc_path, title, output_path):
    cmd = [ "Rscript", os.path.join( get_script_path( ), "plot_roc.r" ),
            title, roc_path, output_path ]

    subprocess.call( cmd )
