#include <iostream>

#include <armadillo>

#include <plinkio/plinkio.h>

#include <gzstream/gzutil.hpp>
#include <bayesic/covariates.hpp>
#include <glm/irls.hpp>
#include <glm/models/binomial.hpp>
#include <cpp-argparse/OptionParser.h>

#include <plink/plink_file.hpp>
#include <bayesic/pair_iter.hpp>
#include <bayesic/method/bayesic_method.hpp>
#include <bayesic/method/bayesic_fine_method.hpp>
#include <bayesic/method/logistic_method.hpp>
#include <bayesic/method/logistic_factor_method.hpp>
#include <bayesic/method/loglinear_method.hpp>
#include <bayesic/method/method.hpp>

using namespace arma;
using namespace optparse;

const std::string USAGE = "bayesic [-c covariates.csv] [-m method] pairs genotype_plink_prefix";
const std::string VERSION = "Bayesic 0.0.1";
const std::string DESCRIPTION = "A tool for inferring genetic interactions.";
const std::string EPILOG = "";

std::vector<snp_row>
create_genotype_matrix(plink_file_ptr genotype_file)
{
    std::vector<snp_row> genotype_matrix;
    snp_row row;
    while( genotype_file->next_row( row ) )
    {
        genotype_matrix.push_back( row );
    }

    return genotype_matrix;
}

vec
create_phenotype_vector(plink_file_ptr genotype_file, uvec &missing)
{
    const std::vector<pio_sample_t> &samples = genotype_file->get_samples( );
    vec phenotype( samples.size( ) );
    for(int i = 0; i < samples.size( ); i++)
    {
        phenotype[ i ] = samples[ i ].phenotype;
        if( samples[ i ].affection == PIO_MISSING )
        {
            missing[ i ] = 1;
        }
    }

    return phenotype;
}

unsigned int
count_interactions(const char *pair_file_path, const std::vector<std::string> &loci)
{
    shared_ptr<std::istream> pair_file = open_possible_gz( pair_file_path );
    pair_iter pairs( *pair_file, loci );
    unsigned int num_interactions = 0;
    std::pair<size_t, size_t> pair;
    while( pairs.get_pair( &pair ) )
    {
        num_interactions++;
    }

    return num_interactions;
}

int
main(int argc, char *argv[])
{
    OptionParser parser = OptionParser( ).usage( USAGE )
                                         .version( VERSION )
                                         .description( DESCRIPTION )
                                         .epilog( EPILOG );
    
    parser.add_option( "-c", "--cov" ).action( "store" ).type( "string" ).metavar( "filename" ).help( "Performs the analysis by including the covariates in this file." );
    
    char const* const choices[] = { "bayes", "logistic", "factor", "loglinear", "fine" };
    parser.add_option( "-m", "--method" ).choices( &choices[ 0 ], &choices[ 5 ] ).metavar( "method" ).help( "Which method to use, one of: 'bayes', 'logistic' or 'loglinear'." );

    parser.add_option( "-n" ).type( "int" ).help( "The number of interactions to correct for." ).set_default( 1 );
    parser.add_option( "-p", "--pheno" ).help( "Read phenotypes from this file instead of a plink file." );
    parser.add_option( "-i", "--mc-iterations" ).type( "int" ).help( "The number of monte carlo iterations to use in the fine method." ).set_default( 50000 );

    Values options = parser.parse_args( argc, argv );
    std::vector<std::string> args = parser.args( );
    if( args.size( ) != 2 )
    {
        printf( "bayesic: error: Pairs or genotypes is missing.\n" );
        parser.print_help( );
        exit( 1 );
    }

    /* Read all genotypes */
    plink_file_ptr genotype_file = open_plink_file( args[ 1 ] );
    std::vector<snp_row> genotype_matrix = create_genotype_matrix( genotype_file );
    
    /* Create pair iterator */
    shared_ptr<std::istream> pair_file = open_possible_gz( args[ 0 ].c_str( ) );
    std::vector<std::string> locus_names = genotype_file->get_locus_names( );
    pair_iter pairs( *pair_file, locus_names );
    
    /* Read additional data  */
    method_data_ptr data( new method_data( ) );
    data->missing = zeros<uvec>( genotype_file->get_samples( ).size( ) );
    data->phenotype = create_phenotype_vector( genotype_file, data->missing );
    std::vector<std::string> order = genotype_file->get_sample_iids( );
    if( options.is_set( "pheno" ) )
    {
        std::ifstream phenotype_file( options[ "p" ].c_str( ) );
        data->phenotype = parse_phenotypes( phenotype_file, data->missing, order );
    }
    if( options.is_set( "cov" ) )
    {
        std::ifstream covariate_file( options[ "c" ].c_str( ) );
        data->covariate_matrix = parse_covariate_matrix( covariate_file, data->missing, order );
    }

    /* Count the number of interactions to adjust for */
    data->num_interactions = (unsigned int) options.get( "n" );
    if( !options.is_set( "n" ) )
    {
        data->num_interactions = count_interactions( args[ 0 ].c_str( ), locus_names );
    }

    /* Run method */
    if( options[ "method" ] == "bayes" )
    {
        bayesic_method bayesic( data );
        run_method( bayesic, genotype_matrix, locus_names, pairs );
    }
    else if( options[ "method" ] == "logistic" )
    {
        logistic_method logistic( data );
        run_method( logistic, genotype_matrix, locus_names, pairs );
    }
    else if( options[ "method" ] == "factor" )
    {
        logistic_factor_method logistic( data );
        run_method( logistic, genotype_matrix, locus_names, pairs );
    }
    else if( options[ "method" ] == "loglinear" )
    {
        loglinear_method loglinear( data );
        run_method( loglinear, genotype_matrix, locus_names, pairs );
    }
    if( options[ "method" ] == "fine" )
    {
        bayesic_fine_method bayesic_fine( data, (int) options.get( "mc_iterations" ) );
        run_method( bayesic_fine, genotype_matrix, locus_names, pairs );
    }

    return 0;
}
