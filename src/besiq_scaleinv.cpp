#include <iostream>

#include <armadillo>

#include <glm/irls.hpp>
#include <glm/models/binomial.hpp>
#include <cpp-argparse/OptionParser.h>

#include <besiq/method/scaleinv_method.hpp>
#include <besiq/method/boxcox_method.hpp>
#include <besiq/method/method.hpp>

#include "common_options.hpp"

using namespace arma;
using namespace optparse;

const std::string USAGE = "besiq-scaleinv [OPTIONS] pairs genotype_plink_prefix";
const std::string DESCRIPTION = "Scale invariant tests (multiple link function tests).";

int
main(int argc, char *argv[])
{
    OptionParser parser = create_common_options( USAGE, DESCRIPTION, true );
    
    char const* const model_choices[] = { "binomial", "normal" };
    char const* const factor_choices[] = { "factor", "additive", "tukey" };

    OptionGroup group = OptionGroup( parser, "Options for glm", "These options will change the behavior of the glm model." );
    group.add_option( "-m", "--model" ).choices( &model_choices[ 0 ], &model_choices[ 2 ] ).metavar( "model" ).help( "The model to use for the phenotype, 'binomial' or 'normal', default = 'binomial'." ).set_default( "binomial" );
    group.add_option( "-f", "--factor" ).choices( &factor_choices[ 0 ], &factor_choices[ 3 ] ).help( "Determines how to code the SNPs, in 'factor' no order of the alleles is assumed, in 'additive' the SNPs are coded as the number of minor alleles, in 'tukey' the coding is the same as factor except that a single parameter for the interaction is used." ).set_default( "factor" );
    group.add_option( "-b", "--box-cox" ).action( "store_true" ).help(  "Runs a box-cox:ish methodology using a family of link functions." );
    group.add_option( "--bc-start" ).set_default( -2.0 ).help(  "Start lambda for box-cox (default=-2.0)." );
    group.add_option( "--bc-end" ).set_default( 3.0 ).help(  "End lambda for box-cox (default=3.0)." );
    group.add_option( "--bc-step" ).set_default( 0.5 ).help(  "Step lambda for box-cox (default=0.5)." );
    group.add_option( "--power-odds" ).action( "store_true" ).help(  "Use the power-odds family for continuous data." );
    parser.add_option_group( group );

    Values options = parser.parse_args( argc, argv );
    if( parser.args( ).size( ) != 2 )
    {
        parser.print_help( );
        exit( 1 );
    }
    shared_ptr<common_options> parsed_data = parse_common_options( options, parser.args( ) );
    parsed_data->data->phenotype.elem( arma::find_nonfinite( parsed_data->data->phenotype ) ).zeros( );

    model_matrix *model_matrix = make_model_matrix( options[ "factor" ], parsed_data->data->covariate_matrix, parsed_data->data->phenotype.n_elem );

    method_type *m = NULL;

    if( !options.is_set( "box_cox" ) )
    {
        if( options[ "model" ] == "binomial" )
        {
            m = new scaleinv_method( parsed_data->data, *model_matrix, false );
        }
        else if( options[ "model" ] == "normal" )
        {
            m = new scaleinv_method( parsed_data->data, *model_matrix, true );
        }
    }
    else
    {
        float lambda_start = (float) options.get( "bc_start" );
        float lambda_end = (float) options.get( "bc_end" );
        float lambda_step = (float) options.get( "bc_step" );

        if( options[ "model" ] == "binomial" )
        {
            m = new boxcox_method( parsed_data->data, *model_matrix, false, lambda_start, lambda_end, lambda_step );
        }
        else if( options[ "model" ] == "normal" )
        {
            m = new boxcox_method( parsed_data->data, *model_matrix, true, lambda_start, lambda_end, lambda_step, options.is_set( "power_odds" ) );
        }
    }
    
    run_method( *m, parsed_data->genotypes, *parsed_data->pairs, *parsed_data->result_file );

    delete m;
    delete model_matrix;

    return 0;
}
