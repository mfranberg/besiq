#include <string>
#include <utility>
#include <iostream>

#include <bayesic/pair_iter.hpp>

pair_iter::pair_iter(std::istream &pair_stream, const std::vector<std::string> &loci)
    : m_pair_stream( pair_stream )
{
    for(int i = 0; i < loci.size( ); i++)
    {
        m_snp_to_index[ loci[ i ] ] = i;
    }
}

bool
pair_iter::get_pair(std::pair<size_t, size_t> *next_pair)
{
    while( m_pair_stream.good( ) )
    {
        std::string name1;
        std::string name2;
        if( !( m_pair_stream >> name1 ) || !( m_pair_stream >> name2 ) )
        {
            return false;
        }

        if( m_snp_to_index.count( name1 ) > 0 && m_snp_to_index.count( name2 ) > 0 )
        {
            next_pair->first = m_snp_to_index[ name1 ];
            next_pair->second = m_snp_to_index[ name2 ];
            return true;
        }
    }

    return false;
}
