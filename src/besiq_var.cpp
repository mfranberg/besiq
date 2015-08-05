#include <iostream>

#include <armadillo>

#include <cpp-argparse/OptionParser.h>
#include <besiq/io/covariates.hpp>

#include <plink/plink_file.hpp>
#include <dcdflib/libdcdf.hpp>

using namespace arma;
using namespace optparse;

const std::string USAGE = "besiq-var [OPTIONS] plink_file";
const std::string DESCRIPTION = "Computes variance hetrogenity for single variants.";
const std::string VERSION = "besiq 0.0.1";
const std::string EPILOG = "";

arma::mat
compute_medians(const snp_row &row, const arma::vec &pheno, const arma::uvec &missing)
{
    std::vector<std::vector<double>> groups( 3 );
    arma::mat medians = arma::zeros<arma::mat>( 3, 2 );
    for(int i = 0; i < row.size( ); i++)
    {
        if( row[ i ] != 3 && missing[ i ] == 0 )
        {
            groups[ row[ i ] ].push_back( pheno[ i ] );
            medians( row[ i ], 1 ) += 1.0;
        }
    }

    for(int i = 0; i < 3; i++)
    {
        std::sort( groups[ i ].begin( ), groups[ i ].end( ) );
        if( groups[ i ].size( ) > 0 )
        {
            medians( i, 0 ) = groups[ i ][ groups[ i ].size( ) / 2 ];
        }
        else
        {
            medians( i, 0 ) = 0.0;
        }
    }

    return medians;
}

bool
compute_brown_forsythe(const snp_row &row, const arma::vec &pheno, const arma::uvec &missing, double *W, double *p, size_t *N)
{
    double k = 3;
    arma::mat medians = compute_medians( row, pheno, missing );

    if( arma::min( medians.col( 1 ) ) <= 20 )
    {
        *N = arma::accu( medians.col( 1 ) );
        return false;
    }
    
    arma::vec z_i = arma::zeros<arma::vec>( 3 );
    double z = 0.0;
    *N = arma::accu( medians.col( 1 ) );
    for(int i = 0; i < row.size( ); i++)
    {
        if( row[ i ] != 3 && missing[ i ] == 0 )
        {
            double z_ij = std::abs( pheno[ i ] - medians( row[ i ], 0 ) );
            z_i[ row[ i ] ] += z_ij / medians( row[ i ], 1 );
            z += z_ij / *N;
        }
    }

    double W_sq = 0.0;
    for(int i = 0; i < row.size( ); i++)
    {
        if( row[ i ] != 3 && missing[ i ] == 0 )
        {
            double z_ij = std::abs( pheno[ i ] - medians( row[ i ], 0 ) );
            W_sq += pow( z_ij - z_i[ row[ i ] ], 2 );
        }
    }

    double numerator = dot( medians.col( 1 ), pow( z_i - z, 2 ) );

    *W = ( (*N - k) * numerator ) / ( ( k - 1 ) * W_sq );
    *p = 1 - f_cdf( *W, k - 1, *N - k );

    return true;
}


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
    arma::vec phenotypes;
    if( options.is_set( "pheno" ) )
    {
        std::ifstream phenotype_file( options[ "pheno" ].c_str( ) );
        phenotypes = parse_phenotypes( phenotype_file, missing, order, options[ "mpheno" ] );
    }
    else
    {
        phenotypes = create_phenotype_vector( genotype_file->get_samples( ), missing );
    }


    std::vector<pio_locus_t> loci = genotype_file->get_loci( );
    std::cout << "chr\tpos\tsnp\tW\tP\tN\n";
    for( int i = 0; i < loci.size( ); i++)
    {
        const snp_row &row = *genotypes->get_row( loci[ i ].name );
        double W;
        double p;
        size_t N;
        if( compute_brown_forsythe( row, phenotypes, missing, &W, &p, &N ) )
        {
            std::cout << (int) loci[ i ].chromosome << "\t" << loci[ i ].bp_position << "\t"  << loci[ i ].name << "\t" << W << "\t" << p << "\t" << N << "\n";
        }
        else
        {
            std::cout << (int) loci[ i ].chromosome << "\t" << loci[ i ].bp_position << "\t" << loci[ i ].name << "\t" << "NA" << "\t" << "NA" << "\t" << N << "\n";
        }
    }

    return 0;
}
