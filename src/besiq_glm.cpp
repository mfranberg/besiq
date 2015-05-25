#include <iostream>

#include <armadillo>

#include <glm/models/binomial.hpp>
#include <glm/models/normal.hpp>
#include <cpp-argparse/OptionParser.h>

#include <besiq/method/glm_method.hpp>
#include <besiq/method/method.hpp>

#include "common_options.hpp"

using namespace arma;
using namespace optparse;

const std::string USAGE = "besiq-glm [OPTIONS] pairs genotype_plink_prefix";
const std::string DESCRIPTION = "Generalized linear models for genetic interactions.";

int
main(int argc, char *argv[])
{
    OptionParser parser = create_common_options( USAGE, DESCRIPTION, true );
 
    char const* const model_choices[] = { "binomial", "normal" };
    char const* const link_choices[] = { "logit", "logc", "odds", "identity", "log" };
    char const* const factor_choices[] = { "factor", "additive", "tukey" };

    OptionGroup group = OptionGroup( parser, "Options for glm", "These options will change the behavior of the glm model." );
    group.add_option( "-m", "--model" ).choices( &model_choices[ 0 ], &model_choices[ 2 ] ).metavar( "model" ).help( "The model to use for the phenotype, 'binomial' or 'normal', default = 'binomial'." ).set_default( "binomial" );
    group.add_option( "-l", "--link-function" ).choices( &link_choices[ 0 ], &link_choices[ 5 ] ).metavar( "link" ).help( "The link function, or scale, that is used for the penetrance: 'logit' log(p/(1-p)), 'logc' log(1 - p), 'odds' p/(1-p), 'identity' p, 'log' log(p)." );
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
        std::string link = "logit";
        if( options.is_set( "link_function" ) )
        {
            link = options[ "link_function" ];
        }

        binomial *model = new binomial( link );
        m = new glm_method( parsed_data->data, *model, *model_matrix );
    }
    else if( options[ "model" ] == "normal" )
    {
        std::string link = "identity";
        if( options.is_set( "link_function" ) )
        {
            link = options[ "link_function" ];
        }

        normal *model = new normal( link );
        m = new glm_method( parsed_data->data, *model, *model_matrix );
    }

    run_method( *m, parsed_data->genotypes, *parsed_data->pairs, *parsed_data->result_file );

    delete m;
    delete model_matrix;

    return 0;
}
