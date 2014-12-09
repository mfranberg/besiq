#include <algorithm>
#include <iterator>
#include <iostream>
#include <map>
#include <sstream>
#include <stdexcept>
#include <vector>

#include <assert.h>

#include <armadillo>

#include <bayesic/covariates.hpp>

using namespace arma;

/**
 * Tokenizes a header, making sure that it starts with FID and IID.
 *
 * @param header A line that represents the header of a CSV file.
 *
 * @throws std::runtime_error.
 * 
 * @return The tokenized header.
 */
std::vector<std::string>
get_fields(const std::string &header)
{
    std::istringstream header_stream( header );
    std::vector<std::string> fields;
    std::copy( std::istream_iterator<std::string>( header_stream ),
               std::istream_iterator<std::string>( ),
               std::back_inserter< std::vector<std::string> >( fields ) );

    if( fields.size( ) <= 0 || fields[ 0 ] != "FID" || fields[ 1 ] != "IID" )
    {
        throw std::runtime_error( "get_fields: The first two fields must be named FID and IID." );
    }

    return fields;
}

/**
 * Parses a single field expecting it to be a double.
 *
 * @param field_str The field string.
 * @param row_num The line number.
 * @param column The column number.
 *
 * @throws std::runtime_error.
 *
 * @return A parsed double.
 */
double
parse_field(const std::string &field_str, unsigned int row_num, unsigned int column)
{
    std::istringstream field_stream( field_str );
    float field;
    char c;
    if( ! (field_stream >> field) || field_stream.get( c ) )
    {
        std::ostringstream error_message;
        error_message << "Could not parse file, error on line: " << row_num + 2 << " column " << column;
        throw std::runtime_error( error_message.str( ) );
    }
    else
    {
        return field;
    }
}

/**
 * Creates a map from iid to its row in the matrix.
 *
 * @param order A list of iids in the order they should appear
 *              in the parsed covariate matrix.
 *
 * @return A map from iid to row number.
 */
std::map<std::string, size_t>
create_iid_map(const std::vector<std::string> &order)
{
    std::map<std::string, size_t> iid_index;
    for(int i = 0; i < order.size( ); i++)
    {
           iid_index[ order[ i ] ] = i;
    }

    return iid_index;
}

mat
parse_covariate_matrix(std::istream &stream, arma::uvec &missing, const std::vector<std::string> &order, const char *missing_string)
{
    std::string header;
    std::getline( stream, header );
    std::vector<std::string> header_fields = get_fields( header );
    
    std::map<std::string, size_t> iid_index = create_iid_map( order );
    mat X = arma::zeros<arma::mat>( order.size( ), header_fields.size( ) - 2 );

    std::string line;
    int row_num = 0;
    rowvec row( header_fields.size( ) - 2 );
    while( std::getline( stream, line ) )
    {
        std::istringstream line_stream( line );
        std::string fid;
        std::string iid;

        line_stream >> fid;
        line_stream >> iid;
        if( iid_index.count( iid ) <= 0 )
        {
            continue;
        }

        for(int i = 2; i < header_fields.size( ); i++)
        {
            std::string field_str;
            if( ! (line_stream >> field_str) )
            {
                std::ostringstream error_message;
                error_message << "Missing column on line " << row_num + 2;
                throw std::runtime_error( error_message.str( ) );
            }

            if( field_str != missing_string )
            {
                row[ i - 2 ] = parse_field( field_str, row_num, i );
            }
            else
            {
                missing[ iid_index[ iid ] ] = 1;
                row[ i - 2 ] = 0.0;
            }
        }

        X.row( iid_index[ iid ] ) = row;
        iid_index.erase( iid );
        
        row_num++;
    }

    /* Any missing iids marked as missing */
    if( iid_index.size( ) > 0 )
    {
        std::map<std::string, size_t>::iterator it;
        for(it = iid_index.begin( ); it != iid_index.end( ); ++it)
        {
            missing[ it->second ] = 1;
        }
    }

    return X;
}

arma::vec
parse_phenotypes(std::istream &stream, arma::uvec &missing, const std::vector<std::string> &order, const char *missing_string)
{
    mat phenotype_matrix = parse_covariate_matrix( stream, missing, order, missing_string );
    assert( phenotype_matrix.n_cols == 1 );
    return phenotype_matrix.col( 0 );
}

arma::mat
parse_environment(std::istream &stream, arma::uvec &missing, const std::vector<std::string> &order, unsigned int levels = 1, const char *missing_string)
{
    std::string header;
    std::getline( stream, header );
    std::vector<std::string> header_fields = get_fields( header );
    std::map<std::string, int> level_map;
    if( header_fields.size( ) != 3 )
    {
        throw std::runtime_error( "parse_environment: Environment file must contain exactly 3 columns." );
    }
    
    std::map<std::string, size_t> iid_index = create_iid_map( order );
    mat X = arma::zeros<arma::mat>( order.size( ), levels );

    std::string line;
    int row_num = 0;
    rowvec row( levels );
    int current_level = 0;
    while( std::getline( stream, line ) )
    {
        std::istringstream line_stream( line );
        std::string fid;
        std::string iid;
        std::string env;

        line_stream >> fid;
        line_stream >> iid;
        line_stream >> env;
        
        if( iid_index.count( iid ) <= 0 )
        {
            continue;
        }

        if( env == missing_string )
        {
            missing[ iid_index[ iid ] ] = 1;
            row = zeros<rowvec>( levels );
        }

        if( levels == 1 )
        {
            row[ 0 ] = parse_field( env, row_num, 2 );
        }
        else
        {
            row = zeros<rowvec>( levels );
            if( level_map.count( env ) == 0 )
            {
                level_map[ env ] = current_level;
                current_level++;

                if( current_level > levels )
                {
                    throw std::runtime_error( "parse_environment: The environmental factor has more levels than specified." );
                }
            }
            row[ level_map[ env ] ] = 1.0;
        }

        X.row( iid_index[ iid ] ) = row;
        iid_index.erase( iid );
        
        row_num++;
    }

    /* Any missing iids marked as missing */
    if( iid_index.size( ) > 0 )
    {
        std::map<std::string, size_t>::iterator it;
        for(it = iid_index.begin( ); it != iid_index.end( ); ++it)
        {
            missing[ it->second ] = 1;
        }
    }

    return X;
}

