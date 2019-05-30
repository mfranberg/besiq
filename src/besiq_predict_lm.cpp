#include <iostream>
#include <sstream>

#include <armadillo>

#include <cpp-argparse/OptionParser.h>
#include <besiq/io/covariates.hpp>

#include <plink/plink_file.hpp>
#include <glm/glm.hpp>
#include <glm/models/normal.hpp>

#include "gene_environment.hpp"

using namespace arma;
using namespace optparse;

const std::string USAGE = "besiq-predict-lm [OPTIONS] plink_file vars_file";
const std::string DESCRIPTION = "Models phenotypes given a list of variables.";
const std::string VERSION = "besiq 0.0.1";
const std::string EPILOG = "";

std::vector<std::string>
parse_vars(std::istream &stream)
{
    std::string line;
    std::set<std::string> vars;
    while( std::getline( stream, line ) )
    {
        std::istringstream line_stream( line );
        std::string rs_name;
        line_stream >> rs_name;
        vars.insert( rs_name ); 

        std::istringstream var_stream( rs_name );
        std::string var_name;
        while( std::getline( var_stream, var_name, ':' ) )
        {
            vars.insert( var_name );
        }
    }

    std::vector<std::string> all_vars;
    all_vars.assign( vars.begin( ), vars.end( ) );
 
    return all_vars;
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
        else
        {
            std::cerr << "besiq-predict-lm: Error could not find: '" << labels[ i ] << "'.\n";
            exit( 1 );
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
    parser.add_option( "-p", "--pheno" ).help( "Read phenotypes from this file instead of a plink file." );
    parser.add_option( "-o", "--out" ).help( "The output file that will contain the results (binary)." );
    parser.add_option( "--matrix" ).help( "Write the model matrix instead." ).action( "store_true" );

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

    std::ifstream var_file( parser.args( )[ 1 ].c_str( ) );
    std::vector<std::string> labels = parse_vars( var_file );

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

    arma::mat cov;
    std::vector<std::string> cov_names;
    if( options.is_set( "cov" ) )
    {
        std::ifstream covariate_file( options[ "cov" ].c_str( ) );
        cov = parse_covariate_matrix( covariate_file, missing_not_used, order, &cov_names );
    }

    /* Open output stream */
    std::ofstream output_file;
    if( options.is_set( "out" ) )
    {
        output_file.open( options[ "out" ].c_str( ) );
    }
    std::ostream &out = options.is_set( "out" ) ? output_file : std::cout;

    gene_environment ge( genotypes, cov, phenotype, cov_names, false );
    ge.impute_missing( );

    arma::uvec causal = get_causal_indices( ge.get_names( ), labels );
 
    arma::mat X = ge.get_active( causal );
    arma::vec y = ge.get_centered_phenotype( );
    
    if( options.is_set("matrix") )
    {
        out << "pheno,";
        for(int i = 0; i < causal.size( ); i++)
        {
            out << ge.get_name( causal[ i ] );
            if(i < causal.size( ) - 1)
            {
                out << ",";
            }
        }
        out << "\n";

        X.insert_cols( 0, y );
        X.save( out, arma::csv_ascii );
    }
    else
    {
        X.insert_cols( 0, arma::ones<arma::vec>( X.n_rows ) );
        
        glm_info result;
        normal model( "identity" );
        arma::uvec missing = arma::zeros<arma::uvec>( genotype_file->get_samples( ).size( ) );
        arma::vec beta = glm_fit( X, y, missing, model, result);

        out << "snp\tbeta\tse_beta\tpvalue\n";
        if( result.success && result.converged )
        {
            for(int i = 1; i < beta.n_elem; i++)
            {
                out << ge.get_name( causal[ i - 1 ] ) << "\t" <<
                    beta[ i ] << "\t" <<
                    result.se_beta[ i ] << "\t" <<
                    result.p_value[ i ] << "\n";
            }
        }
        else
        {
            std::cerr << "beisq-predict-lm: Couldn't converge.\n";
        }
    }

    return 0;
}
