#include <iostream>
#include <string>
#include <vector>

#include <plink_file.hpp>
#include <OptionParser.h>

#include <stats/snp_count.hpp>


using namespace optparse;

const std::string USAGE = "bayesic-pairs genotype_plink_prefix";
const std::string VERSION = "bayesic-pairs 1.0.0";
const std::string DESCRIPTION = "Generates a list of interactions to test with Bayesic.";
const std::string EPILOG = "";

std::vector<float>
compute_maf(plink_file_ptr &genotype_file)
{
    std::vector<float> maf_vec;
    snp_row row;
    while( genotype_file->next_row( row ) )
    {
        double maf = compute_real_maf( row );
        if( maf > 0.5 )
        {
            maf = 1.0 - maf;
        }

        maf_vec.push_back( maf );
    }

    return maf_vec;
}

typedef bool (*maf_validator)(double snp1_maf, double snp2_maf, double maf_threshold);

bool all_ok(double snp1_maf, double snp2_maf, double maf_threshold)
{
    return true;
}

bool maf_ok(double snp1_maf, double snp2_maf, double maf_threshold)
{
    return snp1_maf >= maf_threshold && snp2_maf >= maf_threshold;
}

bool combined_maf_ok(double snp1_maf, double snp2_maf, double maf_threshold)
{
    return snp1_maf * snp2_maf >= maf_threshold;
}

int
main(int argc, char *argv[])
{
    OptionParser parser = OptionParser( ).usage( USAGE )
                                         .version( VERSION )
                                         .description( DESCRIPTION )
                                         .epilog( EPILOG );
    
    parser.add_option( "-m", "--maf" ).type( "float" ).help( "Remove pairs where one of the SNPs have a maf less than this." );
    parser.add_option( "-c", "--combined-maf" ).type( "float" ).help( "Remve pairs where the product of the MAFs is less than this." );

    Values options = parser.parse_args( argc, argv );
    std::vector<std::string> args = parser.args( );
    if( args.size( ) != 1 )
    {
        printf( "bayesic-pairs: error: Pairs or genotypes is missing.\n" );
        parser.print_help( );
        exit( 1 );
    }

    if( options.is_set( "maf" ) && options.is_set( "combined_maf" ) )
    {
        printf( "bayesic-pairs: error: Option --maf and --combined-maf is mututally exclusive.\n" );
        exit( 1 );
    }

    maf_validator validate_maf = all_ok;
    double maf_threshold = 0.0;
    if( options.is_set( "maf" ) )
    {
        validate_maf = maf_ok;
        maf_threshold = (double) options.get( "maf" );
    }
    if( options.is_set( "combined_maf" ) )
    {
        validate_maf = combined_maf_ok;
        maf_threshold = (double) options.get( "combined_maf" );
    }

    plink_file_ptr genotype_file = open_plink_file( args[ 0 ] );
    std::vector<float> maf_vec = compute_maf( genotype_file );
    std::vector<pio_locus_t> loci = genotype_file->get_loci( );
    for(int i = 0; i < loci.size( ); i++)
    {
        for(int j = i + 1; j < loci.size( ); j++)
        {
            if( validate_maf( maf_vec[ i ], maf_vec[ j ],  maf_threshold ) )
            {
                printf( "%s %s\n", loci[ i ].name, loci[ j ].name );
            }
        }
    }

    return 1;
}
