#include <iostream>

#include <armadillo>

#include <plinkio/plinkio.h>

#include <gzstream/gzutil.hpp>
#include <bayesic/io/covariates.hpp>
#include <cpp-argparse/OptionParser.h>

#include <plink/plink_file.hpp>
#include <bayesic/env_method/lm_env_stepwise.hpp>
#include <bayesic/method/method.hpp>
#include <bayesic/method/env_method.hpp>

using namespace arma;
using namespace optparse;

const std::string USAGE = "bayesic-env [-c covariates.csv] [-m method] environment_file genotype_plink_prefix";
const std::string VERSION = "Bayesic 0.0.1";
const std::string DESCRIPTION = "A tool for inferring variant-environment interactions.";
const std::string EPILOG = "";

int
main(int argc, char *argv[])
{
    OptionParser parser = OptionParser( ).usage( USAGE )
                                         .version( VERSION )
                                         .description( DESCRIPTION )
                                         .epilog( EPILOG );
    
    
    char const* const choices[] = { "stepwise" };
    char const* const link_choices[] = { "logistic", "log-complement", "odds-additive", "penetrance-additive", "penetrance-multiplicative" };
    char const* const factor_choices[] = { "factor", "additive", "tukey" };

    parser.add_option( "-m", "--method" ).choices( &choices[ 0 ], &choices[ 1 ] ).metavar( "method" ).help( "Which method to use, one of: 'stepwise'." );
    parser.add_option( "-p", "--pheno" ).help( "Read phenotypes from this file instead of a plink file." );
    parser.add_option( "-n", "--mpheno" ).help( "Name of the phenotype to use." );
    parser.add_option( "-e", "--levels" ).type( "int" ).help( "The number of levels of the environmental variable." ).set_default( 1 );
    parser.add_option( "-f", "--factor" ).choices( &factor_choices[ 0 ], &factor_choices[ 3 ] ).help( "Determines how to code the SNPs, in 'factor' no order of the alleles is assumed, in 'additive' the SNPs are coded as the number of minor alleles, in 'tukey' the coding is the same as factor except that a single parameter for the interaction is used." ).set_default( "factor" );
    parser.add_option( "-c", "--cov" ).action( "store" ).type( "string" ).metavar( "filename" ).help( "Performs the analysis by including the covariates in this file." );
    
    Values options = parser.parse_args( argc, argv );
    std::vector<std::string> args = parser.args( );
    if( args.size( ) != 2 )
    {
        std::cerr << "bayesic: error: Genotypes or environment factor is missing." << std::endl;
        parser.print_help( );
        exit( 1 );
    }
    else if( !options.is_set( "method" ) )
    {
        std::cerr << "bayesic:error: No method selected." << std::endl;
        exit( 1 );
    }

    /* Read all genotypes */
    plink_file_ptr genotype_file = open_plink_file( args[ 1 ] );
    genotype_matrix_ptr genotypes = create_genotype_matrix( genotype_file );
    std::vector<std::string> locus_names = genotype_file->get_locus_names( );
 
    /* Read additional data  */
    method_data_ptr data( new method_data( ) );
    data->missing = zeros<uvec>( genotype_file->get_samples( ).size( ) );
    std::vector<std::string> order = genotype_file->get_sample_iids( );
    if( options.is_set( "pheno" ) )
    {
        std::ifstream phenotype_file( options[ "pheno" ].c_str( ) );
        data->phenotype = parse_phenotypes( phenotype_file, data->missing, order, options[ "mpheno" ] );
    }
    else
    {
        data->phenotype = create_phenotype_vector( genotype_file->get_samples( ), data->missing );
    }

    /* Read covariates */
    if( options.is_set( "cov" ) )
    {
        std::ifstream covariate_file( options[ "cov" ].c_str( ) );
        data->covariate_matrix = parse_covariate_matrix( covariate_file, data->missing, order );
    }

    /* Read environment factor */
    std::ifstream environment_file( args[ 0 ].c_str( ) );
    arma::mat E = parse_environment( environment_file, data->missing, order, (unsigned int) options.get( "levels" ) );

    /* XXX: Implement proper log file. */
    std::ostream nullstream( 0 );
    arma::set_stream_err2( nullstream );

    if( options[ "method" ] == "stepwise" )
    {
        if( E.n_cols > 1 )
        {
            /* We just pass the non-redundant levels */
            E = E.cols( 1, E.n_cols - 1 );
        }

        lm_env_stepwise stepwise_env( data, E );
        run_env_method( stepwise_env, genotypes );
    }

    return 0;
}
