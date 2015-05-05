#include <iostream>

#include <armadillo>

#include <cpp-argparse/OptionParser.h>

#include <bayesic/method/caseonly_method.hpp>
#include <bayesic/method/peer_method.hpp>
#include <bayesic/method/method.hpp>

#include "common_options.hpp"

using namespace arma;
using namespace optparse;

const std::string USAGE = "bayesic-caseonly [OPTIONS] pairs genotype_plink_prefix";
const std::string DESCRIPTION = "A stage-wise case-only filter for genetic interactions.";

int
main(int argc, char *argv[])
{
    char const* const method_choices[] = { "css", "r2", "contrast", "peer" };
    OptionParser parser = create_common_options( USAGE, DESCRIPTION, false );
    parser.add_option( "-m", "--method" ).choices( &method_choices[ 0 ], &method_choices[ 4 ] ).metavar( "method" ).help( "The type of test statistic to compute 'r2' or 'css'." ).set_default( "css" );
    
    Values options = parser.parse_args( argc, argv );
    if( parser.args( ).size( ) != 2 )
    {
        parser.print_help( );
        exit( 1 );
    }
    shared_ptr<common_options> parsed_data = parse_common_options( options, parser.args( ) );

    method_type *m = NULL;
    if( options[ "method" ] != "peer" )
    {
        m = new caseonly_method( parsed_data->data, options[ "method" ] );
    }
    else if( options[ "method" ] == "peer" )
    {
        m = new peer_method( parsed_data->data );
    }
    
    run_method( *m, parsed_data->genotypes, *parsed_data->pairs, *parsed_data->result_file );

    delete m;

    return 0;
}
