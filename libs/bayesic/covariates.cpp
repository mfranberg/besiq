#include <algorithm>
#include <iterator>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <vector>

#include <armadillo>

#include <covariates.hpp>

using namespace arma;

class bad_conversion : public std::runtime_error
{
    public:
        bad_conversion(std::string const& s)
        : std::runtime_error(s)
        {

        }
};

std::vector<std::string>
get_fields(const std::string &header)
{
    std::istringstream header_stream( header );
    std::vector<std::string> fields;
    std::copy( std::istream_iterator<std::string>( header_stream ),
               std::istream_iterator<std::string>( ),
               std::back_inserter< std::vector<std::string> >( fields ) );

    return fields;
}

mat
parse_covariate_matrix(std::istream &stream, const char *missing_string)
{
    std::string header;
    std::getline( stream, header );
    std::vector<std::string> header_fields = get_fields( header );
    std::vector<bool> missing;

    mat X;
    std::string line;
    int row_num = 0;
    rowvec row( header_fields.size( ) );
    while( std::getline( stream, line ) )
    {
        missing.push_back( false );

        std::istringstream line_stream( line );
        for(int i = 0; i < header_fields.size( ); i++)
        {
            std::string field_str;
            line_stream >> field_str;
            if( field_str != missing_string )
            {
                std::istringstream field_stream( field_str );
                float field;
                char c;
                if( ! (field_stream >> field) || field_stream.get( c ) )
                {
                    std::ostringstream error_message;
                    error_message << "Could not parse covariate file, error on line: " << row_num + 2 << " column " << i;
                    throw bad_conversion( error_message.str( ) );
                }

                row[ i ] = field;
            }
            else
            {
                missing[ row_num ] = true;
                row[ i ] = 0.0;
            }
        }

        X.insert_rows( row_num, row );
        
        row_num++;
    }

    std::cout << X << std::endl;
    std::cout << X.n_cols << std::endl;
    std::cout << X.n_rows << std::endl;

    return X;
}
