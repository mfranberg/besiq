#include <iostream>
#include <sstream>

#include <armadillo>

#include <cpp-argparse/OptionParser.h>
#include <besiq/io/covariates.hpp>
#include <glm/irls.hpp>
#include <glm/models/normal.hpp>

#include <plink/plink_file.hpp>

#include "gene_environment.hpp"

using namespace arma;
using namespace optparse;

const std::string USAGE = "besiq-gxe [OPTIONS] plink_file variable_file";
const std::string DESCRIPTION = "Calculates the full regression model (only linear).";
const std::string VERSION = "besiq 0.0.1";
const std::string EPILOG = "";

std::vector<std::string>
parse_variables(std::istream &stream)
{
    std::string line;
    std::vector<std::string> variable_list;
    while( std::getline( stream, line ) )
    {
        std::istringstream line_stream( line );
        std::string variable_name;

        line_stream >> variable_name;
            
        variable_list.push_back( variable_name );
    }

    return variable_list;
}

arma::uvec get_included_indices(const std::vector<std::string> variable_order, const std::vector<std::string> &labels)
{
    std::vector<size_t> indices;
    std::map<std::string, int> variable_to_index;
    for(int i = 0; i < variable_order.size( ); i++)
    {
        variable_to_index[ variable_order[ i ] ] = i;
    }

    for(int i = 0; i < labels.size( ); i++)
    {
        if( variable_to_index.count( labels[ i ] ) > 0 )
        {
            indices.push_back( variable_to_index[ labels[ i ] ] );
        }
    }

    arma::uvec result = arma::zeros<arma::uvec>( indices.size( ) );
    for(int i = 0; i < indices.size( ); i++)
    {
        result[ i ] = indices[ i ];
    }

    return result;
}

int
main(int argc, char *argv[])
{
    OptionParser parser = OptionParser( ).usage( USAGE )
                                         .version( VERSION )
                                         .description( DESCRIPTION )
                                         .epilog( EPILOG ); 
 
    parser.add_option( "-p", "--pheno" ).help( "Read phenotypes from this file instead of a plink file." );
    parser.add_option( "-e", "--env" ).action( "store" ).type( "string" ).metavar( "filename" ).help( "Environemental variables." );
    parser.add_option( "-o", "--out" ).help( "The output file that will contain the results (binary)." );
    parser.add_option( "--standardize" ).help( "Use standardized genotypes and covariates." ).action( "store_true" );

    Values options = parser.parse_args( argc, argv );
    if( parser.args( ).size( ) != 2 )
    {
        parser.print_help( );
        exit( 1 );
    }

    std::ios_base::sync_with_stdio( false );
    
    /* Read all genotypes */
    plink_file_ptr genotype_file = open_plink_file( parser.args( )[ 0 ] );
    genotype_matrix_ptr genotypes = create_genotype_matrix( genotype_file );
    std::vector<std::string> order = genotype_file->get_sample_iids( );
    std::vector< std::pair<std::string, std::string> > fid_iid = genotype_file->get_sample_fid_iid( );

    std::ifstream variable_file( parser.args( )[ 1 ].c_str( ) );
    std::vector<std::string> variable_list = parse_variables( variable_file );

    /* Make error streams separate from stdout */
    arma::set_stream_err1( std::cerr );
    arma::set_stream_err2( std::cerr );

    /* Parse phenotypes */
    arma::uvec missing_not_used = arma::zeros<arma::uvec>( genotype_file->get_samples( ).size( ) );
    arma::vec phenotype;
    if( options.is_set( "pheno" ) )
    {
        std::ifstream phenotype_file( options[ "pheno" ].c_str( ) );
        phenotype = parse_phenotypes( phenotype_file, missing_not_used, order, options[ "mpheno" ] );
    }
    else
    {
        phenotype = create_phenotype_vector( genotype_file->get_samples( ), missing_not_used );
    }
    
    /* Parse env */
    arma::uvec env_missing = arma::zeros<arma::uvec>( genotype_file->get_samples( ).size( ) );
    arma::mat env;
    std::vector<std::string> env_names;
    if( options.is_set( "env" ) )
    {
        std::ifstream env_file( options[ "env" ].c_str( ) );
        env = parse_covariate_matrix( env_file, env_missing, order, &env_names );
    }

    /* Open output stream */
    std::ofstream output_file;
    if( options.is_set( "out" ) )
    {
        output_file.open( options[ "out" ].c_str( ) );
    }
    std::ostream &out = options.is_set( "out" ) ? output_file : std::cout;

    gene_environment ge( genotypes, env, phenotype, env_names, false );
    ge.impute_missing( );

    arma::uvec included = get_included_indices( ge.get_names( ), variable_list );
    if( included.n_elem != variable_list.size( ) )
    {
        std::cout << "besiq-gxe: Error could not find all genotype or covariate variables.\n";
        return 1;
    }

    bool standardize = options.is_set( "only_pvalues" );
    
    arma::mat X;
    if( standardize )
    {
        X = ge.get_active( included );
    }
    else
    {
        X = ge.get_active_raw( included );
    }

    glm_info result;
    normal model( "identity" );
    arma::uvec missing = arma::zeros<arma::uvec>( genotype_file->get_samples( ).size( ) );
    arma::vec beta = irls( X, ge.get_centered_phenotype(), missing, model, result);

    out << "snp\tbeta\tse_beta\tpvalue\n";
    if( result.success && result.converged )
    {
        for(int i = 0; i < included.n_elem; i++)
        {
            out << ge.get_name( included[ i ] ) << "\t" <<
                beta[ i ] << "\t" <<
                result.se_beta[ i ] << "\t" <<
                result.p_value[ i ] << "\n";
        }
    }

    return 0;
}
