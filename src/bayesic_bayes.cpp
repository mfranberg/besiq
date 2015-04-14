#include <iostream>

#include <armadillo>

#include <cpp-argparse/OptionParser.h>

#include <bayesic/prior.hpp>
#include <bayesic/method/bayesic_method.hpp>
#include <bayesic/method/bayesic_fine_method.hpp>

#include "common_options.hpp"

using namespace arma;
using namespace optparse;

const std::string USAGE = "bayesic-bayes [OPTIONS] pairs genotype_plink_prefix";
const std::string DESCRIPTION = "Bayesian inference of genetic interactions.";

int
main(int argc, char *argv[])
{
    OptionParser parser = create_common_options( USAGE, DESCRIPTION, false );
    
    OptionGroup group = OptionGroup( parser, "Options for bayes", "These options will change the behaviour of bayes and fine." );
    group.add_option( "-n", "--num-interactions" ).type( "int" ).help( "The number of interactions to correct for, this is used in the model prior (default: all)." );
    group.add_option( "-s", "--num-single" ).type( "int" ).help( "The number of snps to consider when correcting (default: proportional to square of the number of interactions)." );
    group.add_option( "-t", "--single-prior" ).type( "float" ).help( "The probability that a single snp is associated (default: %default)." ).set_default( 0.0 );
    group.add_option( "-i", "--mc-iterations" ).type( "int" ).help( "The number of monte carlo iterations to use in the fine method (default: %default)." ).set_default( 4000000 );
    group.add_option( "-a", "--beta-prior-param1" ).type( "float" ).help( "First shape parameter of beta prior (default: %default)." ).set_default( 2.0 );
    group.add_option( "-b", "--beta-prior-param2" ).type( "float" ).help( "Second shape parameter of beta prior (default: %default)." ).set_default( 2.0 );
    group.add_option( "-e", "--estimate-prior-params" ).action( "store_true" ).help( "Estimate prior parameters from data by permuting phenotype (default: off)." );
    parser.add_option( "--additive" ).action( "store_true" ).help( "Use an additive model (slow)." );
    parser.add_option_group( group );

    Values options = parser.parse_args( argc, argv );
    if( parser.args( ).size( ) != 2 )
    {
        parser.print_help( );
        exit( 1 );
    }

    shared_ptr<common_options> parsed_data = parse_common_options( options, parser.args( ) );

    /* Read prior parameters */
    arma::vec alpha = arma::ones<arma::vec>( 2 );
    alpha[ 0 ] = (float) options.get( "beta_prior_param1" );
    alpha[ 1 ] = (float) options.get( "beta_prior_param2" );
    if( options.is_set( "estimate_prior_params" ) )
    {
        alpha = estimate_prior_parameters( parsed_data->genotypes, parsed_data->data->phenotype, parsed_data->data->missing, 5000 );
    }

    /* Count the number of interactions to adjust for */
    parsed_data->data->single_prior = (float) options.get( "single_prior" );
    parsed_data->data->num_single = (unsigned int) options.get( "num_single" );
    parsed_data->data->num_interactions = (unsigned int) options.get( "num_interactions" );

    method_type *m = new bayesic_method( parsed_data->data, alpha );
    if( options.is_set( "additive" ) )
    {
        m = new bayesic_method( parsed_data->data, alpha );
    }
    else
    {
        m = new bayesic_fine_method( parsed_data->data, (int) options.get( "mc_iterations" ), alpha );
    }
    
    run_method( *m, parsed_data->genotypes, *parsed_data->pairs, *parsed_data->result_file );

    delete m;

    return 0;
}
