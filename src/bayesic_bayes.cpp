#include <iostream>

#include <armadillo>

#include <cpp-argparse/OptionParser.h>

#include <plink/plink_file.hpp>
#include <bayesic/io/covariates.hpp>
#include <bayesic/io/pairfile.hpp>
#include <bayesic/io/resultfile.hpp>
#include <bayesic/prior.hpp>
#include <bayesic/method/bayesic_method.hpp>
#include <bayesic/method/bayesic_fine_method.hpp>

using namespace arma;
using namespace optparse;

const std::string USAGE = "bayesic-bayes [OPTIONS] pairs genotype_plink_prefix";
const std::string VERSION = "Bayesic 0.0.1";
const std::string DESCRIPTION = "Bayesian inference of genetic interactions.";
const std::string EPILOG = "";

int
main(int argc, char *argv[])
{
    OptionParser parser = OptionParser( ).usage( USAGE )
                                         .version( VERSION )
                                         .description( DESCRIPTION )
                                         .epilog( EPILOG );
    
    
    parser.add_option( "-p", "--pheno" ).help( "Read phenotypes from this file instead of a plink file." );
    parser.add_option( "-e", "--mpheno" ).help( "Name of the phenotype that you want to read (if there are more than one in the phenotype file)." );
    parser.add_option( "-o", "--out" ).help( "The output file that will contain the results (binary)." );
    parser.add_option( "--print-params" ).action( "store_true" ).set_default( 0 ).help( "Print parameter estimates in factor GLM models." );
    
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

    /* Read prior parameters */
    arma::vec alpha = arma::ones<arma::vec>( 2 );
    alpha[ 0 ] = (float) options.get( "beta_prior_param1" );
    alpha[ 1 ] = (float) options.get( "beta_prior_param2" );
    if( options.is_set( "estimate_prior_params" ) )
    {
        alpha = estimate_prior_parameters( genotypes, data->phenotype, data->missing, 5000 );
    }

    /* Count the number of interactions to adjust for */
    data->single_prior = (float) options.get( "single_prior" );
    data->num_single = (unsigned int) options.get( "num_single" );
    data->num_interactions = (unsigned int) options.get( "num_interactions" );

    /* XXX: Implement proper log file. */
    std::ostream nullstream( 0 );
    arma::set_stream_err2( nullstream );

    method_type *m = new bayesic_method( data, alpha );
    if( options.is_set( "additive" ) )
    {
        m = new bayesic_method( data, alpha );
    }
    else
    {
        m = new bayesic_fine_method( data, (int) options.get( "mc_iterations" ), alpha );
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
