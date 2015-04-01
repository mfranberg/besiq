#include <iostream>
#include <iomanip>

#include <cpp-argparse/OptionParser.h>

#include <bayesic/io/resultfile.hpp>

using namespace optparse;

const std::string USAGE = "bayesic-view result_file [result_file2 ...]";
const std::string VERSION = "Bayesic 0.0.1";
const std::string DESCRIPTION = "A tool for viewing binary result file.";
const std::string EPILOG = "";

struct comparator
{
    virtual bool operator()(const float& x, const float&y) const = 0;
};

struct return_true : comparator
{
    bool operator()(const float& x, const float&y) const { return true; };
};
struct less : comparator
{
    bool operator()(const float& x, const float&y) const { return x < y; };
};
struct less_equal : comparator
{
    bool operator()(const float& x, const float&y) const { return x <= y; };
};
struct greater : comparator
{
    bool operator()(const float& x, const float&y) const { return x > y; };
};
struct greater_equal : comparator
{
    bool operator()(const float& x, const float&y) const { return x >= y; };
};

int
main(int argc, char *argv[])
{
    OptionParser parser = OptionParser( ).usage( USAGE )
                                         .version( VERSION )
                                         .description( DESCRIPTION )
                                         .epilog( EPILOG ); 

    parser.add_option( "-c", "--count" ).action( "store_true" ).set_default( 0 ).help( "Count the number of pairs and exit." );
    char const * const operations[] = { "none", "lt", "le", "gt", "ge" };
    parser.add_option( "-p", "--operation" ).set_default( "none" ).choices( &operations[ 0 ], &operations[ 5 ] ).help( "The filtering operation to use 'none', 'lt', 'le', 'gt' or 'ge' (default = none)." );
    parser.add_option( "-t", "--threshold" ).set_default( 0.05 ).help( "Filter using this threshold (default = 0.05)." );
    parser.add_option( "-f", "--field" ).set_default( 0 ).help( "The value field to filter on, the field index of the first non snp name is 0." );
    
    Values options = parser.parse_args( argc, argv );
    std::vector<std::string> args = parser.args( );
    if( args.size( ) < 1 )
    {
        std::cerr << "bayesic-view: error: Need at least one result file." << std::endl;
        parser.print_help( );
        exit( 1 );
    }

    std::vector<bresultfile *> result_files;
    for(int i = 0; i < args.size( ); i++)
    {
        bresultfile *result = new bresultfile( args[ i ] );
        result_files.push_back( result );
        if( !result->open( ) )
        {
            std::cerr << "bayesic-view: error: Could not open result file: " << args[ i ] << std::endl;
        }
    }
   
    if( (bool) options.get( "count" ) )
    {
        for(int i = 0; i < result_files.size( ); i++)
        {
            std::cout << result_files[ i ]->num_pairs( ) << std::endl;
        }
        return 0;
    }

    float threshold = (float) options.get( "threshold" );
    size_t field = (size_t) options.get( "field" );
    std::map< std::string, comparator * > op;
    op[ "none" ] = new return_true( );
    op[ "lt" ] = new less( );
    op[ "le" ] = new less_equal( );
    op[ "gt" ] = new greater( );
    op[ "ge" ] = new greater_equal( );

    comparator &compare = *op[ (std::string) options.get( "operation" ) ];

    size_t header_size = result_files[ 0 ]->get_header( ).size( );
    for(int i = 1; i < result_files.size( ); i++)
    {
        if( result_files[ i ]->get_header( ).size( ) != header_size )
        {
            std::cerr << "bayesic-view: error: Different number of columns in result files." << std::endl;
            return 1;
        }
    }

    std::setprecision( 4 );
    std::cout << "snp1 snp2";
    std::vector<std::string> header = result_files[ 0 ]->get_header( );
    for(int i = 0; i < header.size( ); i++)
    {
        std::cout << "\t" << header[ i ];
    }
    std::cout << std::endl;

    for(int i = 0; i < result_files.size( ); i++)
    {
        bresultfile &res = *result_files[ i ];

        float *output = new float[ header_size ];
        std::pair<std::string, std::string> pair;

        while( res.read( &pair, output ) )
        {
            if( !compare( output[ field ], threshold ) )
            {
                continue;
            }

            std::cout << pair.first << " " << pair.second;
            for(int j = 0; j < header_size; j++)
            {
                float value = output[ j ];
                if( output[ j ] != result_get_missing( ) )
                {
                    std::cout << "\t" << output[ j ];
                }
                else
                {
                    std::cout << "\tNA";
                }
            }
            std::cout << std::endl;
        }

        delete output;
    }
    
    return 0;
}
