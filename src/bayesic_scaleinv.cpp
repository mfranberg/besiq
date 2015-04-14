#include <iostream>

#include <armadillo>

#include <glm/irls.hpp>
#include <glm/models/binomial.hpp>
#include <cpp-argparse/OptionParser.h>

#include <bayesic/method/scaleinv_method.hpp>
#include <bayesic/method/method.hpp>

#include "common_options.hpp"

using namespace arma;
using namespace optparse;

const std::string USAGE = "bayesic-scaleinv [OPTIONS] pairs genotype_plink_prefix";
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
    parser.add_option_group( group );

    Values options = parser.parse_args( argc, argv );
    if( parser.args( ).size( ) != 2 )
    {
        parser.print_help( );
        exit( 1 );
    }
    shared_ptr<common_options> parsed_data = parse_common_options( options, parser.args( ) );

    model_matrix *model_matrix = make_model_matrix( options[ "factor" ], parsed_data->data->covariate_matrix, parsed_data->data->phenotype.n_elem );

    method_type *m = NULL;
    if( options[ "model" ] == "binomial" )
    {
        m = new scaleinv_method( parsed_data->data, *model_matrix, false );
    }
    else if( options[ "model" ] == "normal" )
    {
        m = new scaleinv_method( parsed_data->data, *model_matrix, true );
    }
    
    run_method( *m, parsed_data->genotypes, *parsed_data->pairs, *parsed_data->result_file );

    delete m;
    delete model_matrix;

    return 0;
}
