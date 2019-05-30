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

const std::string USAGE = "besiq-extract [OPTIONS] plink_file vars_file";
const std::string DESCRIPTION = "Extract genotypes, covariates and phenotype into csv.";
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
    
void
find_geno_and_cov_indices(genotype_matrix_ptr genotypes,
                          const std::vector<std::string> cov_names,
                          const std::vector<std::string> labels,
                          std::vector<unsigned int> &geno_indices,
                          std::vector<unsigned int> &cov_indices)
{
    for(int i = 0; i < labels.size( ); i++)
    {
        int geno_index = genotypes->get_index( labels[ i ] );
        std::vector<std::string>::const_iterator cov_elem = std::find( cov_names.begin( ), cov_names.end( ), labels[ i ] );
        if( geno_index != -1 )
        {
            geno_indices.push_back( geno_index );
        }
        else if( cov_elem != labels.end( ) )
        {
            cov_indices.push_back( cov_elem - cov_names.begin( ) );
        }
        else
        {
            std::cerr << "besiq-extract: Error could not find: '" << labels[ i ] << "'.\n";
            exit( 1 );
        }
    }
}

std::string
geno_to_string(snp_t s)
{
    if( s == 3 )
    {
        return std::string( "NA" );
    }
    else
    {
        return std::to_string( s );
    }
}

std::string
cov_to_string(double cov)
{
    if( cov != cov )
    {
        return std::string( "NA" );
    }

    double intpart;
    if( std::modf(cov, &intpart) == 0.0 )
    {
        return std::to_string( (int) cov );
    }
    else
    {
        return std::to_string( cov );
    }
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


    std::vector<std::string> loci = genotypes->get_snp_names( );
    std::vector<unsigned int> geno_indices;
    std::vector<unsigned int> cov_indices;
    find_geno_and_cov_indices( genotypes, cov_names, labels, geno_indices, cov_indices );

    out << "FID,IID,pheno";
    for(int i = 0; i < geno_indices.size(); i++)
    {
        out << "," << loci[ geno_indices[ i ] ];
    }
    for(int i = 0; i < cov_indices.size( ); i++)
    {
        out << "," << cov_names[ cov_indices[ i ] ];
    }
    out << "\n";

    for(int i = 0; i < phenotype.n_elem; i++)
    {
        out << fid_iid[i].first;
        out << "," << fid_iid[i].second;
        out << "," << phenotype[ i ];
        for(int j = 0; j < geno_indices.size( ); j++)
        {
            out << "," << geno_to_string( genotypes->get_row( geno_indices[ j ] )[ i ] );
        }
        for(int j = 0; j < cov_indices.size( ); j++)
        {
            out << "," << cov_to_string( cov( i, cov_indices[ j ] - 2 ) ); // -2 for FID and IID
        }
        out << "\n";
    }

    return 0;
}
