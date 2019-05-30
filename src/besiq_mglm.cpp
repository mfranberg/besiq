#include <iostream>
#include <sstream>
#include <cfloat>
#include <set>

#include <armadillo>

#include <cpp-argparse/OptionParser.h>
#include <besiq/io/covariates.hpp>
#include <besiq/stats/snp_count.hpp>
#include <glm/irls.hpp>
#include <glm/models/normal.hpp>
#include <glm/models/binomial.hpp>

#include <plink/plink_file.hpp>
#include <dcdflib/libdcdf.hpp>

using namespace arma;
using namespace optparse;

const std::string USAGE = "besiq-mglm [OPTIONS] plink_file";
const std::string DESCRIPTION = "Multivariate additive GLM.";
const std::string VERSION = "besiq 0.0.1";
const std::string EPILOG = "";

std::set<std::string>
parse_variants(std::ifstream &input_file)
{
    std::set<std::string> variant_set;
    std::string variant;
    while( input_file >> variant )
    {
        variant_set.insert( variant );
    }

    return variant_set;
}

std::vector<std::string>
create_design_matrix(genotype_matrix_ptr genotypes, const arma::mat &cov, const std::set<std::string> &variant_set, arma::uvec &missing, arma::mat &X)
{

    bool all_variants = variant_set.size( ) == 0;
    size_t num_variants = variant_set.size( );
    if( all_variants )
    {
        num_variants = genotypes->size( );
    }

    X.resize( missing.n_elem, num_variants + cov.n_cols + 1 );

    std::vector<std::string> names = genotypes->get_snp_names( );
    std::vector<std::string> valid_names;
    size_t cur_col = 0;
    for(int i = 0; i < genotypes->size( ); i++)
    {
        if( !all_variants && variant_set.count( names[ i ] ) <= 0 )
        {
            continue;
        }
        
        snp_row &cur_row = genotypes->get_row( i );
        for(int j = 0; j < cur_row.size( ); j++)
        {
            if( cur_row[ j ] != 3 )
            {
                X( j, cur_col ) = cur_row[ j ];   
            }
            else
            {
                X( j, cur_col ) = 0;
                missing[ j ] = 1;
            }
        }

        if( arma::var( X.col( cur_col ) ) > 1e-4 )
        {
            valid_names.push_back( names[ i ] );
            cur_col++;
        }
    }

    for(int i = 0; i < cov.n_cols; i++)
    {
        X.col( cur_col ) = cov.col( i );
        cur_col++;
    }

    X.col( cur_col ) = arma::ones<arma::vec>( missing.n_elem );
    cur_col++;
    X.resize( missing.n_elem, cur_col );

    return valid_names;
}

int
main(int argc, char *argv[])
{
    char const* const model_choices[] = { "binomial", "normal" };
    char const* const link_choices[] = { "logit", "logc", "odds", "identity", "log" };
    OptionParser parser = OptionParser( ).usage( USAGE )
                                         .version( VERSION )
                                         .description( DESCRIPTION )
                                         .epilog( EPILOG ); 
 
    parser.add_option( "-p", "--pheno" ).help( "Read phenotypes from this file instead of a plink file." );
    parser.add_option( "-e", "--mpheno" ).help( "Name of the phenotype that you want to read (if there are more than one in the phenotype file)." );
    parser.add_option( "-c", "--cov" ).action( "store" ).type( "string" ).metavar( "filename" ).help( "Performs the analysis by including the covariates in this file." );
    parser.add_option( "-s", "--variants" ).action( "store" ).type( "string" ).metavar( "filename" ).help( "Only use these variants." );
    parser.add_option( "-o", "--out" ).help( "The output file that will contain the results (binary)." );
    parser.add_option( "-t", "--threshold" ).help( "Only output p-values less or equal to this value." ).set_default( 1.0 );
    parser.add_option( "-m", "--model" ).choices( &model_choices[ 0 ], &model_choices[ 2 ] ).metavar( "model" ).help( "The model to use for the phenotype, 'binomial' or 'normal', default = 'binomial'." ).set_default( "binomial" );
    parser.add_option( "-l", "--link-function" ).choices( &link_choices[ 0 ], &link_choices[ 5 ] ).metavar( "link" ).help( "The link function, or scale, that is used for the penetrance: 'logit' log(p/(1-p)), 'logc' log(1 - p), 'odds' p/(1-p), 'identity' p, 'log' log(p)." );

    Values options = parser.parse_args( argc, argv );
    if( parser.args( ).size( ) != 1 )
    {
        parser.print_help( );
        exit( 1 );
    }

    std::ios_base::sync_with_stdio( false );
    
    /* Read all genotypes */
    plink_file_ptr genotype_file = open_plink_file( parser.args( )[ 0 ] );
    genotype_matrix_ptr genotypes = create_genotype_matrix( genotype_file );
    std::vector<std::string> order = genotype_file->get_sample_iids( );

    /* Parse phenotypes */
    arma::uvec missing = arma::zeros<arma::uvec>( genotype_file->get_samples( ).size( ) );
    arma::vec phenotype;
    arma::mat cov;
    std::vector<std::string> cov_names;
    if( options.is_set( "pheno" ) )
    {
        std::ifstream phenotype_file( options[ "pheno" ].c_str( ) );
        phenotype = parse_phenotypes( phenotype_file, missing, order, options[ "mpheno" ] );
    }
    else
    {
        phenotype = create_phenotype_vector( genotype_file->get_samples( ), missing );
    }
    if( options.is_set( "cov" ) )
    {
        std::ifstream covariate_file( options[ "cov" ].c_str( ) );
        cov = parse_covariate_matrix( covariate_file, missing, order, &cov_names );
    }

    std::set<std::string> variant_set;
    if( options.is_set( "variants" ) )
    {
        std::ifstream variant_file( options[ "variants" ].c_str( ) );
        variant_set = parse_variants( variant_file );
    }

    /* Open output stream */
    std::ofstream output_file;
    if( options.is_set( "out" ) )
    {
        output_file.open( options[ "out" ].c_str( ) );
    }
    std::ostream &out = options.is_set( "out" ) ? output_file : std::cout;

    arma::mat X;
    std::vector<std::string> variant_names = create_design_matrix( genotypes, cov, variant_set, missing, X );
    
    glm_model *model = NULL;
    if( options[ "model" ] == "binomial" )
    {
        std::string link = "logit";
        if( options.is_set( "link_function" ) )
        {
            link = options[ "link_function" ];
        }

        model= new binomial( link );
    }
    else if( options[ "model" ] == "normal" )
    {
        std::string link = "identity";
        if( options.is_set( "link_function" ) )
        {
            link = options[ "link_function" ];
        }

        model = new normal( link );
    }

    glm_info result;
    arma::vec beta = irls( X, phenotype, missing, *model, result);

    /* The design matrix is constructed such that the variants appear first */
    out << "snp\tbeta\tse_beta\tpvalue\n";
    if( result.success && result.converged )
    {
        for(int i = 0; i < variant_names.size( ); i++)
        {
            out << variant_names[ i ] << "\t" <<
                beta[ i ] << "\t" <<
                result.se_beta[ i ] << "\t" <<
                result.p_value[ i ] << "\n";
        }
    }

    if( model != NULL )
    {
        delete model; 
    }

    return 0;
}
