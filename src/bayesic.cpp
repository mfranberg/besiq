#include <iostream>

#include <armadillo>

#include <gzstream/gzutil.hpp>
#include <glm/irls.hpp>
#include <glm/models/binomial.hpp>
#include <cpp-argparse/OptionParser.h>

#include <plink/plink_file.hpp>
#include <bayesic/io/covariates.hpp>
#include <bayesic/io/pairfile.hpp>
#include <bayesic/io/resultfile.hpp>
#include <bayesic/prior.hpp>
#include <bayesic/method/bayesic_method.hpp>
#include <bayesic/method/bayesic_fine_method.hpp>
#include <bayesic/method/lm_method.hpp>
#include <bayesic/method/glm_method.hpp>
#include <bayesic/method/loglinear_method.hpp>
#include <bayesic/method/caseonly_method.hpp>
#include <bayesic/method/stepwise_method.hpp>
#include <bayesic/method/lm_stepwise_method.hpp>
#include <bayesic/method/wald_method.hpp>
#include <bayesic/method/wald_lm_method.hpp>
#include <bayesic/method/method.hpp>

using namespace arma;
using namespace optparse;

const std::string USAGE = "bayesic [-c covariates.csv] [-m method] pairs genotype_plink_prefix";
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
    
    
    char const* const choices[] = { "bayes", "bayes-fine", "lm", "glm", "loglinear", "caseonly", "stepwise", "lm-stepwise", "wald", "wald-lm" };
    char const* const link_choices[] = { "logit", "logc", "odds", "identity", "log" };
    char const* const factor_choices[] = { "factor", "additive", "tukey" };

    parser.add_option( "-m", "--method" ).choices( &choices[ 0 ], &choices[ 10 ] ).metavar( "method" ).help( "Which method to use, one of: 'bayes', 'bayes-fine', 'lm', 'glm', 'loglinear', 'caseonly', 'stepwise', 'lm-stepwise', 'wald', or 'wald-lm'." );
    parser.add_option( "-p", "--pheno" ).help( "Read phenotypes from this file instead of a plink file." );
    parser.add_option( "-e", "--mpheno" ).help( "Name of the phenotype that you want to read (if there are more than one in the phenotype file)." );
    parser.add_option( "-o", "--out" ).help( "The output file that will contain the results (binary)." );
    parser.add_option( "-l", "--link-function" ).choices( &link_choices[ 0 ], &link_choices[ 5 ] ).metavar( "link" ).help( "The link function, or scale, that is used for the penetrance: 'logit' log(p/(1-p)), 'logc' log(1 - p), 'odds' p/(1-p), 'identity' p, 'log' log(p)." ).set_default( "logit" );
    parser.add_option( "-f", "--factor" ).choices( &factor_choices[ 0 ], &factor_choices[ 3 ] ).help( "Determines how to code the SNPs, in 'factor' no order of the alleles is assumed, in 'additive' the SNPs are coded as the number of minor alleles, in 'tukey' the coding is the same as factor except that a single parameter for the interaction is used." ).set_default( "factor" );

    parser.add_option( "-c", "--cov" ).action( "store" ).type( "string" ).metavar( "filename" ).help( "Performs the analysis by including the covariates in this file." );
    parser.add_option( "--print-params" ).action( "store_true" ).set_default( 0 ).help( "Print parameter estimates in factor GLM models." );
    
    OptionGroup group = OptionGroup( parser, "Options for bayes", "These options will change the behaviour of bayes and fine." );
    group.add_option( "-n", "--num-interactions" ).type( "int" ).help( "The number of interactions to correct for, this is used in the model prior (default: all)." );
    group.add_option( "-s", "--num-single" ).type( "int" ).help( "The number of snps to consider when correcting (default: proportional to square of the number of interactions)." );
    group.add_option( "-t", "--single-prior" ).type( "float" ).help( "The probability that a single snp is associated (default: %default)." ).set_default( 0.0 );
    group.add_option( "-i", "--mc-iterations" ).type( "int" ).help( "The number of monte carlo iterations to use in the fine method (default: %default)." ).set_default( 4000000 );
    group.add_option( "-a", "--beta-prior-param1" ).type( "float" ).help( "First shape parameter of beta prior (default: %default)." ).set_default( 2.0 );
    group.add_option( "-b", "--beta-prior-param2" ).type( "float" ).help( "Second shape parameter of beta prior (default: %default)." ).set_default( 2.0 );
    group.add_option( "-e", "--estimate-prior-params" ).action( "store_true" ).help( "Estimate prior parameters from data by permuting phenotype (default: off)." );
    parser.add_option_group( group );

    Values options = parser.parse_args( argc, argv );
    std::vector<std::string> args = parser.args( );
    if( args.size( ) != 2 )
    {
        std::cerr << "bayesic: error: Pairs or genetypes is missing." << std::endl;
        parser.print_help( );
        exit( 1 );
    }
    else if( !options.is_set( "method" ) )
    {
        std::cerr << "bayesic: error: No method selected." << std::endl;
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
        data->phenotype = parse_phenotypes( phenotype_file, data->missing, order, options[ "mpheno" ] );
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

    model_matrix *model_matrix = make_model_matrix( options[ "factor" ], data->covariate_matrix, data->phenotype.n_elem );

    method_type *m = NULL;
    if( options[ "method" ] == "bayes" )
    {
        m = new bayesic_method( data, alpha );
    }
    else if( options[ "method" ] == "bayes-fine" )
    {
        m = new bayesic_fine_method( data, (int) options.get( "mc_iterations" ), alpha );
    }
    else if( options[ "method" ] == "glm" )
    {
        binomial *glm = new binomial( options[ "link_function" ] );
        m = new glm_method( data, *glm, *model_matrix );
    }
    else if( options[ "method" ] == "lm" )
    {
        m = new lm_method( data, *model_matrix );
    }
    else if( options[ "method" ] == "loglinear" )
    {
        m = new loglinear_method( data );
    }
    else if( options[ "method" ] == "caseonly" )
    {
        m = new caseonly_method( data );
    }
    else if( options[ "method" ] == "stepwise" )
    {
        m = new stepwise_method( data );
    }
    else if( options[ "method" ] == "lm-stepwise" )
    {
        m = new lm_stepwise_method( data );
    }
    else if( options[ "method" ] == "wald" )
    {
        m = new wald_method( data );
    }
    else if( options[ "method" ] == "wald-lm" )
    {
        m = new wald_lm_method( data );
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
