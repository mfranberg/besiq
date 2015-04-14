#include <iostream>

#include <armadillo>

#include <cpp-argparse/OptionParser.h>

#include <bayesic/method/loglinear_method.hpp>
#include <bayesic/method/method.hpp>

#include "common_options.hpp"

using namespace arma;
using namespace optparse;

const std::string USAGE = "bayesic-loglinear [OPTIONS] pairs genotype_plink_prefix";
const std::string DESCRIPTION = "A log-linear based test for genetic interactions.";

int
main(int argc, char *argv[])
{
    OptionParser parser = create_common_options( USAGE, DESCRIPTION, false );
    
    Values options = parser.parse_args( argc, argv );
    if( parser.args( ).size( ) != 2 )
    {
        parser.print_help( );
        exit( 1 );
    }
    shared_ptr<common_options> parsed_data = parse_common_options( options, parser.args( ) );

    method_type *m = new loglinear_method( parsed_data->data );
    
    run_method( *m, parsed_data->genotypes, *parsed_data->pairs, *parsed_data->result_file );

    delete m;

    return 0;
}
