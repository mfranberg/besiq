#include <iostream>

#include <plink_file.hpp>
#include <OptionParser.h>

using namespace optparse;

const std::string USAGE = "generate_pairs genotype_plink_prefix";
const std::string VERSION = "generate_pairs 1.0.0";
const std::string DESCRIPTION = "Generates a list of interactions to test with Bayesic.";
const std::string EPILOG = "";

int
main(int argc, char *argv[])
{
    OptionParser parser = OptionParser( ).usage( USAGE )
                                         .version( VERSION )
                                         .description( DESCRIPTION )
                                         .epilog( EPILOG );
    
    parser.add_option( "-m" ).type( "float" ).help( "Remove pairs where one of the SNPs have a maf less than this." );

    Values options = parser.parse_args( argc, argv );
    std::vector<std::string> args = parser.args( );
    if( args.size( ) != 1 )
    {
        printf( "bayesic: error: Pairs or genotypes is missing.\n" );
        parser.print_help( );
        exit( 1 );
    }

    plink_file_ptr genotype_file = open_plink_file( args[ 0 ] );
    std::vector<pio_locus_t> loci = genotype_file->get_loci( );
    for(int i = 0; i < loci.size( ); i++)
    {
        for(int j = i + 1; j < loci.size( ); j++)
        {
            printf( "%s %s\n", loci[ i ].name, loci[ j ].name );
        }
    }

    return 1;
}
