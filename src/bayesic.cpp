#include <armadillo>

#include <covariates.hpp>
#include <irls.hpp>
#include <models/binomial.hpp>
#include <OptionParser.h>

using namespace arma;
using namespace optparse;

const std::string USAGE = "bayesic [-c covariates.csv] [-m method] genotype_plink_prefix";
const std::string VERSION = "Bayesic 0.0.1";
const std::string DESCRIPTION = "A tool for inferring genetic interactions.";
const std::string EPILOG = "";


int
main(int argc, char *argv[])
{
    OptionParser parser = OptionParser( ).usage( USAGE )
                                         .version( VERSION )
                                         .description( DESCRIPTION )
                                         .epilog( EPILOG );
    
    parser.add_option( "-c" ).action( "store" ).type( "string" ).metavar( "filename" ).help( "Performs the analysis by including the covariates in this file." );
    
    char const* const choices[] = { "bayes", "logistic", "loglinear" };
    parser.add_option( "-m" ).choices( &choices[ 0 ], &choices[ 2 ] ).metavar( "method" ).help( "Which method to use, one of: 'bayes', 'logistic' or 'loglinear'." );

    Values options = parser.parse_args( argc, argv );
    if( options.is_set( "c" ) )
    {
        std::ifstream covariate_file( options[ "c" ].c_str( ) );
        std::cout << "reading covariates" << std::endl;
        mat X = parse_covariate_matrix( covariate_file );
    }

    double A_aux[] = { 1.0, 1.0, 1.0, 1.0,
                       -1.160, -0.655, 0.4156, -1.740 };

    double y_aux[] = { 0.1131, 0.4271, 0.9694, 0.0164 };

    mat A( A_aux, 4, 2 );
    vec y( y_aux, 4 );
    binomial binomial_model;

    vec b = irls( A, y, binomial_model );
    std::cout << "b: " << std::endl << b << std::endl;

    return 0;
}
