#include <gzutil.hpp>

bool
ends_with(const std::string &str, const std::string &end)
{
    return str.rfind( end ) == ( str.size( ) - end.size( ) );
}

shared_ptr<std::istream>
open_possible_gz(const std::string &path)
{
    if( ends_with( path, ".gz" ) )
    {
        return shared_ptr<std::istream>( new gz::igzstream( path.c_str( ) ) );
    }
    else
    {
        return shared_ptr<std::istream>( new std::ifstream( path.c_str( ) ) );
    }
}

shared_ptr<std::ostream>
create_possible_gz(const std::string &path)
{
    if( ends_with( path, ".gz" ) )
    {
        return shared_ptr<std::ostream>( new gz::ogzstream( path.c_str( ) ) );
    }
    else
    {
        return shared_ptr<std::ostream>( new std::ofstream( path.c_str( ) ) );
    }
}
