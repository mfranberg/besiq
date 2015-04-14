#include <iostream>

#include <armadillo>

#include <cpp-argparse/OptionParser.h>

#include <bayesic/method/wald_method.hpp>
#include <bayesic/method/wald_lm_method.hpp>
#include <bayesic/method/method.hpp>

#include "common_options.hpp"

using namespace arma;
using namespace optparse;

const std::string USAGE = "bayesic-wald [OPTIONS] pairs genotype_plink_prefix";
const std::string DESCRIPTION = "Fast wald tests for genetic interactions.";

int
main(int argc, char *argv[])
{
    OptionParser parser = create_common_options( USAGE, DESCRIPTION, false );
    
    char const* const model_choices[] = { "binomial", "normal" };
    parser.add_option( "-m", "--model" ).choices( &model_choices[ 0 ], &model_choices[ 2 ] ).metavar( "model" ).help( "The model to use for the phenotype, 'binomial' or 'normal', default = 'binomial'." ).set_default( "binomial" );
    
    Values options = parser.parse_args( argc, argv );
    if( parser.args( ).size( ) != 2 )
    {
        parser.print_help( );
        exit( 1 );
    }
    shared_ptr<common_options> parsed_data = parse_common_options( options, parser.args( ) );

    method_type *m = NULL;
    if( options[ "model" ] == "binomial" )
    {
        m = new wald_method( parsed_data->data );
    }
    else if( options[ "model" ] == "normal" )
    {
        m = new wald_lm_method( parsed_data->data );
    }
    
    run_method( *m, parsed_data->genotypes, *parsed_data->pairs, *parsed_data->result_file );

    delete m;
    
    return 0;
}
