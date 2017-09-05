#include <iostream>
#include <sstream>

#include <armadillo>

#include <cpp-argparse/OptionParser.h>
#include <besiq/io/covariates.hpp>

#include <plink/plink_file.hpp>

#include "gene_environment.hpp"

using namespace arma;
using namespace optparse;

const std::string USAGE = "besiq-predict [OPTIONS] plink_file";
const std::string DESCRIPTION = "Predicts phenotypes given a list of betas.";
const std::string VERSION = "besiq 0.0.1";
const std::string EPILOG = "";

arma::vec
parse_beta(std::istream &stream, std::vector<std::string> &labels)
{
    std::string line;
    std::vector<double> beta_list;
    while( std::getline( stream, line ) )
    {
        std::istringstream line_stream( line );
        std::string rs_name;
        double beta;

        line_stream >> rs_name >> beta;
            
        labels.push_back( rs_name );
        beta_list.push_back( beta );
    }

    arma::vec b = arma::zeros<arma::vec>( beta_list.size( ) );
    for(int i = 0; i < beta_list.size( ); i++)
    {
        b[ i ] = beta_list[ i ];
    }
    
    return b;
}

arma::uvec get_causal_indices(const std::vector<std::string> variable_order, const std::vector<std::string> &labels)
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
 
    parser.add_option( "-c", "--cov" ).action( "store" ).type( "string" ).metavar( "filename" ).help( "Performs the analysis by including the covariates in this file." );
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

    std::vector<std::string> labels;
    std::ifstream beta_file( parser.args( )[ 1 ].c_str( ) );
    arma::vec beta = parse_beta( beta_file, labels );

    /* Make error streams separate from stdout */
    arma::set_stream_err1( std::cerr );
    arma::set_stream_err2( std::cerr );

    /* Parse phenotypes */
    arma::uvec cov_missing = arma::zeros<arma::uvec>( genotype_file->get_samples( ).size( ) );
    arma::mat cov;
    std::vector<std::string> cov_names;
    if( options.is_set( "cov" ) )
    {
        std::ifstream covariate_file( options[ "cov" ].c_str( ) );
        cov = parse_covariate_matrix( covariate_file, cov_missing, order, &cov_names );
    }

    /* Open output stream */
    std::ofstream output_file;
    if( options.is_set( "out" ) )
    {
        output_file.open( options[ "out" ].c_str( ) );
    }
    std::ostream &out = options.is_set( "out" ) ? output_file : std::cout;

    arma::vec phenotype = arma::zeros<arma::vec>( fid_iid.size( ) );
    gene_environment ge( genotypes, cov, phenotype, cov_names, true );
    ge.impute_missing( );

    arma::uvec causal = get_causal_indices( ge.get_names( ), labels );
    if( causal.n_elem != labels.size( ) )
    {
        std::cout << "besiq-predict: Error could not find all genotype or covariate variables.\n";
        return 1;
    }

    bool standardize = options.is_set( "only_pvalues" );
    
    arma::mat X;
    if( standardize )
    {
        X = ge.get_active( causal );
    }
    else
    {
        X = ge.get_active_raw( causal );
    }
    arma::vec y = X * beta;
    arma::vec y_std = ( y - arma::mean( y ) ) / arma::stddev( y );

    out << "FID\tIID\tPhenotype\n";
    for(int i = 0; i < fid_iid.size( ); i++)
    {
        out << fid_iid[ i ].first << "\t" << fid_iid[ i ].second << "\t" << y_std[ i ] << "\n";
    } 

    return 0;
}
