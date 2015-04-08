#include <iostream>

#include <armadillo>

#include <glm/irls.hpp>
#include <glm/models/binomial.hpp>
#include <cpp-argparse/OptionParser.h>

#include <plink/plink_file.hpp>
#include <bayesic/io/covariates.hpp>
#include <bayesic/io/pairfile.hpp>
#include <bayesic/io/resultfile.hpp>
#include <bayesic/method/lm_method.hpp>
#include <bayesic/method/glm_method.hpp>
#include <bayesic/method/method.hpp>

using namespace arma;
using namespace optparse;

const std::string USAGE = "bayesic-glm [OPTIONS] pairs genotype_plink_prefix";
const std::string VERSION = "Bayesic 0.0.1";
const std::string DESCRIPTION = "Generalized linear models for genetic interactions.";
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
    parser.add_option( "-c", "--cov" ).action( "store" ).type( "string" ).metavar( "filename" ).help( "Performs the analysis by including the covariates in this file." );
    parser.add_option( "--print-params" ).action( "store_true" ).set_default( 0 ).help( "Print parameter estimates in factor GLM models." );
    
    char const* const model_choices[] = { "binomial", "normal" };
    char const* const link_choices[] = { "logit", "logc", "odds", "identity", "log" };
    char const* const factor_choices[] = { "factor", "additive", "tukey" };

    OptionGroup group = OptionGroup( parser, "Options for glm", "These options will change the behavior of the glm model." );
    group.add_option( "-m", "--model" ).choices( &model_choices[ 0 ], &model_choices[ 2 ] ).metavar( "model" ).help( "The model to use for the phenotype, 'binomial' or 'normal', default = 'binomial'." ).set_default( "binomial" );
    group.add_option( "-l", "--link-function" ).choices( &link_choices[ 0 ], &link_choices[ 5 ] ).metavar( "link" ).help( "The link function, or scale, that is used for the penetrance: 'logit' log(p/(1-p)), 'logc' log(1 - p), 'odds' p/(1-p), 'identity' p, 'log' log(p)." ).set_default( "logit" );
    group.add_option( "-f", "--factor" ).choices( &factor_choices[ 0 ], &factor_choices[ 3 ] ).help( "Determines how to code the SNPs, in 'factor' no order of the alleles is assumed, in 'additive' the SNPs are coded as the number of minor alleles, in 'tukey' the coding is the same as factor except that a single parameter for the interaction is used." ).set_default( "factor" );
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
    if( options.is_set( "cov" ) )
    {
        std::ifstream covariate_file( options[ "cov" ].c_str( ) );
        data->covariate_matrix = parse_covariate_matrix( covariate_file, data->missing, order );
    }

    /* XXX: Implement proper log file. */
    std::ostream nullstream( 0 );
    arma::set_stream_err2( nullstream );

    model_matrix *model_matrix = make_model_matrix( options[ "factor" ], data->covariate_matrix, data->phenotype.n_elem );

    method_type *m = NULL;
    if( options[ "model" ] == "binomial" )
    {
        binomial *glm = new binomial( options[ "link_function" ] );
        m = new glm_method( data, *glm, *model_matrix );
    }
    else if( options[ "model" ] == "normal" )
    {
        m = new lm_method( data, *model_matrix );
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
    delete model_matrix;
    delete pairs;
    delete result;

    return 0;
}
