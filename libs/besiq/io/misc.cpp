#include <iterator>
#include <sstream>

#include <besiq/io/misc.hpp>

std::vector<std::string>
unpack_string(const char *snp_name_str)
{
    std::string snp_names( snp_name_str );
    std::stringstream strstr( snp_names );
    std::istream_iterator<std::string> it( strstr );
    std::istream_iterator<std::string> end;
    std::vector<std::string> results( it, end );

    return results;
}

std::string
pack_string(const std::vector<std::string> &snp_names)
{
    if( snp_names.size( ) == 0 )
    {
        return std::string( );
    }

    std::stringstream output;
    output << snp_names[ 0 ];
    for(size_t i = 1; i < snp_names.size( ); i++)
    {
        output << "\t" << snp_names[ i ];
    }

    return output.str( );
}
