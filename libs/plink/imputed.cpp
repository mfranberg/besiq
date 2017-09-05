#include <iostream>
#include <fstream>
#include <sstream>

#include <plink/snp_row.hpp>
#include <plink/imputed.hpp>

imputed_matrix_ptr
parse_genotypes(const std::string &path, float call_rate, const std::vector<std::string> &samples)
{
    std::ifstream stream( path.c_str( ) );
    std::string line;
    
    imputed_matrix_ptr matrix( new std::vector<snp_row>( ) );
    while( std::getline( stream, line ) )
    {
        std::string tmp;
        std::string name;
        unsigned long long pos;
        std::string a1;
        std::string a2;
        float p_AA;
        float p_Aa;
        float p_aa;

        std::istringstream line_stream( line );

        line_stream >> tmp;
        line_stream >> name;
        line_stream >> pos;
        line_stream >> a1;
        line_stream >> a2;

        snp_row row;
        row.resize( samples.size( ) );
        int i = 0;
        while( (line_stream >> p_AA) && (line_stream >> p_Aa) && (line_stream >> p_aa) && i < samples.size( ) )
        {
            if( p_AA > p_Aa && p_AA > p_aa && p_AA > call_rate )
            {
                row.assign( i,  0 );
            }
            else if( p_Aa > p_AA && p_Aa > p_aa && p_Aa > call_rate )
            {
                row.assign( i,  1 );
            }
            else if( p_aa > p_AA && p_aa > p_Aa && p_aa > call_rate )
            {
                row.assign( i,  2 );
            }
            else
            {
                row.assign( i,  3 );
            }
            i++;
        }

        matrix->push_back( row );
    }

    return matrix;
}

std::vector<std::string>
parse_sample(const std::string &path)
{
    std::ifstream stream( path.c_str( ) );
    std::vector<std::string> sample_list;
    std::string line;
    
    /* Ignore header lines */
    std::getline( stream, line );
    std::getline( stream, line );

    while( std::getline( stream, line ) )
    {
        std::string fid;
        std::string iid;
        std::istringstream line_stream( line );

        line_stream >> fid;
        line_stream >> iid;

        sample_list.push_back( iid );
    }

    return sample_list;
}

std::vector<imputed_info>
parse_info(const std::string &path)
{
    std::ifstream stream( path.c_str( ) );

    std::vector<imputed_info> info_list;
    std::string line;
    
    /* Ignore header */
    std::getline( stream, line );
    while( std::getline( stream, line ) )
    {
        std::string tmp;
        std::istringstream line_stream( line );

        imputed_info info;
        
        line_stream >> tmp;
        line_stream >> info.name;
        line_stream >> info.pos;
        line_stream >> info.maf;
        line_stream >> info.info;

        info_list.push_back( info );
    }

    return info_list;
}

imputed_data
parse_imputed_data(const std::string &prefix, float call_rate)
{
    imputed_data data;

    data.info = parse_info( prefix + "_info" );
    data.samples = parse_sample( prefix + "_samples" );
    data.genotypes = parse_genotypes( prefix, call_rate, data.samples );
        
    return data;
}
