#include <iostream>

#include <armadillo>

#include <cpp-argparse/OptionParser.h>

#include <plink/plink_file.hpp>
#include <bayesic/io/covariates.hpp>
#include <bayesic/io/pairfile.hpp>
#include <bayesic/io/resultfile.hpp>
#include <bayesic/method/stepwise_method.hpp>
#include <bayesic/method/lm_stepwise_method.hpp>
#include <bayesic/method/method.hpp>

using namespace arma;
using namespace optparse;

const std::string USAGE = "bayesic-stagewise [OPTIONS] pairs genotype_plink_prefix";
const std::string VERSION = "Bayesic 0.0.1";
const std::string DESCRIPTION = "A stage-wise test for genetic interactions.";
const std::string EPILOG = "";

int
main(int argc, char *argv[])
{
    OptionParser parser = OptionParser( ).usage( USAGE )
                                         .version( VERSION )
                                         .description( DESCRIPTION )
                                         .epilog( EPILOG );
    
    char const* const model_choices[] = { "binomial", "normal" };
    parser.add_option( "-p", "--pheno" ).help( "Read phenotypes from this file instead of a plink file." );
    parser.add_option( "-e", "--mpheno" ).help( "Name of the phenotype that you want to read (if there are more than one in the phenotype file)." );

    parser.add_option( "--print-params" ).action( "store_true" ).set_default( 0 ).help( "Print parameter estimates in factor GLM models." );
    parser.add_option( "-m", "--model" ).choices( &model_choices[ 0 ], &model_choices[ 2 ] ).metavar( "model" ).help( "The model to use for the phenotype, 'binomial' or 'normal', default = 'binomial'." ).set_default( "binomial" );

    Values options = parser.parse_args( argc, argv );
    std::vector<std::string> args = parser.args( );
    if( args.size( ) != 2 )
    {
        std::cerr << "bayesic: error: Pairs or genetypes is missing." << std::endl;
        parser.print_help( );
        exit( 1 );
    }

    /* Read all genotypes */
    plink_file_ptr genotype_file = open_plink_file( args[ 1 ] );
    genotype_matrix_ptr genotypes = create_genotype_matrix( genotype_file );

    /* Create pair iterator */
    const std::vector<std::string> &locus_names = genotype_file->get_locus_names( );
    pairfile *pairs = open_pair_file( args[ 0 ].c_str( ), locus_names );
    if( pairs == NULL || !pairs->open( ) )
    {
        std::cerr << "bayesic: error: Could not open pair file." << std::endl;
        exit( 1 );
    }

    /* Read additional data  */
    method_data_ptr data( new method_data( ) );
    data->print_params = (bool) options.get( "print_params" );
    data->missing = zeros<uvec>( genotype_file->get_samples( ).size( ) );
    std::vector<std::string> order = genotype_file->get_sample_iids( );
    if( options.is_set( "pheno" ) )
    {
        std::ifstream phenotype_file( options[ "pheno" ].c_str( ) );
        data->phenotype = parse_phenotypes( phenotype_file, data->missing, order );
    }
    else
    {
        data->phenotype = create_phenotype_vector( genotype_file->get_samples( ), data->missing );
    }

    /* XXX: Implement proper log file. */
    std::ostream nullstream( 0 );
    arma::set_stream_err2( nullstream );

    method_type *m;
    if( options[ "model" ] == "binomial" )
    {
        m = new stepwise_method( data );
    }
    else if( options[ "model" ] == "normal" )
    {
        m = new lm_stepwise_method( data );
    }
    
    /* Run method */
    resultfile *result = NULL;
    if( options.is_set( "out" ) )
    {
        result = new bresultfile( options[ "out" ], locus_names );
    }
    else
    {
        std::ios_base::sync_with_stdio( false );
        result = new tresultfile( "-", "w", locus_names );
    }
    if( !result->open( ) )
    {
        std::cerr << "bayesic: error: Can not open result file." << std::endl;
        exit( 1 );
    }

    run_method( *m, genotypes, *pairs, *result );

    delete m;
    delete pairs;
    delete result;

    return 0;
}
