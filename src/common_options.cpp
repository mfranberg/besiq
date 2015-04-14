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
        std::cerr << "bayesic: error: Pairs or genetypes is missing." << std::endl;
        exit( 1 );
    }
    /* Read all genotypes */
    result->genotype_file = open_plink_file( args[ 1 ] );
    result->genotypes = create_genotype_matrix( result->genotype_file );
    
    /* Create pair iterator */
    size_t split = (size_t) options.get( "split" );
    size_t num_splits = (size_t) options.get( "num_splits" );
    if( split > num_splits || split == 0 || num_splits == 0 )
    {
        std::cerr << "bayesic: error: Num splits and split must be > 0, and split <= num_splits." << std::endl;
        exit( 1 );
    }
    
    result->pairs = shared_ptr<pairfile>( open_pair_file( args[ 0 ].c_str( ), result->genotype_file->get_locus_names( ) ) );
    if( result->pairs == NULL || !result->pairs->open( split, num_splits ) )
    {
        std::cerr << "bayesic: error: Could not open pair file." << std::endl;
        exit( 1 );
    }
    
    /* Read additional data  */
    result->data = shared_ptr<method_data>( new method_data( ) );
    result->data->print_params = (bool) options.get( "print_params" );
    result->data->missing = zeros<uvec>( result->genotype_file->get_samples( ).size( ) );
    std::vector<std::string> order = result->genotype_file->get_sample_iids( );
    if( options.is_set( "pheno" ) )
    {
        std::ifstream phenotype_file( options[ "pheno" ].c_str( ) );
        result->data->phenotype = parse_phenotypes( phenotype_file, result->data->missing, order );
    }
    else
    {
        result->data->phenotype = create_phenotype_vector( result->genotype_file->get_samples( ), result->data->missing );
    }
    if( options.is_set( "cov" ) )
    {
        std::ifstream covariate_file( options[ "cov" ].c_str( ) );
        result->data->covariate_matrix = parse_covariate_matrix( covariate_file, result->data->missing, order );
    }

    /* XXX: Implement proper log file. */
    std::ostream nullstream( 0 );
    arma::set_stream_err2( nullstream );
    
    /* Open results. */
    if( options.is_set( "out" ) )
    {
        result->result_file = shared_ptr<resultfile>( new bresultfile( options[ "out" ], result->genotype_file->get_locus_names( ) ) );
    }
    else
    {
        std::ios_base::sync_with_stdio( false );
        result->result_file = shared_ptr<resultfile>( new tresultfile( "-", "w", result->genotype_file->get_locus_names( ) ) );
    }
    if( !result->result_file->open( ) )
    {
        std::cerr << "bayesic: error: Can not open result file." << std::endl;
        exit( 1 );
    }

    return result;
}

