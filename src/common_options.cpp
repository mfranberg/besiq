#include <armadillo>

#include "common_options.hpp"

const std::string VERSION = "Bayesic 0.5.9";
const std::string EPILOG = "";

using namespace optparse;
using namespace arma;

OptionParser
create_common_options(const std::string &usage, const std::string &description, bool support_cov)
{
    OptionParser parser = OptionParser( ).usage( usage )
                                         .version( VERSION )
                                         .description( description )
                                         .epilog( EPILOG );
    
    parser.add_option( "-p", "--pheno" ).help( "Read phenotypes from this file instead of a plink file." );
    parser.add_option( "-e", "--mpheno" ).help( "Name of the phenotype that you want to read (if there are more than one in the phenotype file)." );
    parser.add_option( "-o", "--out" ).help( "The output file that will contain the results (binary)." );
    parser.add_option( "-c", "--cov" ).action( "store" ).type( "string" ).metavar( "filename" ).help( "Performs the analysis by including the covariates in this file." );
    parser.add_option( "--split" ).help( "Runs the analysis on a part of the pair file, and this is part X of 1-<num_splits> parts (default = 1)." ).set_default( 1 );
    parser.add_option( "--num-splits" ).help( "Sets the number of parts to split the pair file in (default = 1)." ).set_default( 1 );
    parser.add_option( "--print-params" ).action( "store_true" ).set_default( 0 ).help( "Print parameter estimates in result file." );
    
    return parser;
}

shared_ptr<common_options>
parse_common_options(optparse::Values &options, const std::vector<std::string> &args)
{
    shared_ptr<common_options> result;
    if( args.size( ) != 2 )
    {
        std::cerr << "besiq: error: Pairs or genetypes is missing." << std::endl;
        exit( 1 );
    }
    /* Read all genotypes */
    plink_file_ptr genotype_file = open_plink_file( args[ 1 ] );
    genotype_matrix_ptr genotypes = create_genotype_matrix( genotype_file );
    
    /* Create pair iterator */
    size_t split = (size_t) options.get( "split" );
    size_t num_splits = (size_t) options.get( "num_splits" );
    if( split > num_splits || split == 0 || num_splits == 0 )
    {
        std::cerr << "besiq: error: Num splits and split must be > 0, and split <= num_splits." << std::endl;
        exit( 1 );
    }
    
    pairfile *pairs = open_pair_file( args[ 0 ].c_str( ), genotype_file->get_locus_names( ) );
    if( pairs == NULL || !pairs->open( split, num_splits ) )
    {
        std::cerr << "besiq: error: Could not open pair file." << std::endl;
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

    /* XXX: Implement proper log file. */
    arma::set_stream_err1( std::cerr );
    arma::set_stream_err2( std::cerr );
    
    /* Open results. */
    resultfile *result_file = NULL;
    if( options.is_set( "out" ) )
    {
        result_file = new bresultfile( options[ "out" ], genotype_file->get_locus_names( ) );
    }
    else
    {
        std::ios_base::sync_with_stdio( false );
        result_file = new tresultfile( "-", "w" );
    }
    if( result_file == NULL || !result_file->open( ) )
    {
        std::cerr << "besiq: error: Can not open result file." << std::endl;
        exit( 1 );
    }

    return shared_ptr<common_options>( new common_options( genotype_file, genotypes, data, pairs, result_file ) );
}

