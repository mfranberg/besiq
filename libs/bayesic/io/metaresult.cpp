#include <bayesic/io/resultfile.hpp>

#include <bayesic/io/metaresult.hpp>

metaresultfile::metaresultfile(const std::vector<resultfile *> &result_files)
    : m_results( result_files ),
      m_cur_file( 0 )
{
}

bool
metaresultfile::read(std::pair<std::string, std::string> *pair, float *value)
{
    while( m_cur_file < m_results.size( ) )
    {

        if( m_results[ m_cur_file ]->read( pair, value ) )
        {
            return true;
        }
        m_cur_file++;
    }

    return false;
}

uint64_t
metaresultfile::num_pairs()
{
    uint64_t num_pairs = 0;
    for(int i = 0; i < m_results.size( ); i++)
    {
        num_pairs += m_results[ i ]->num_pairs( );
    }

    return num_pairs;
}

std::vector<std::string>
metaresultfile::get_header()
{
    if( m_results.size( ) > 0 )
    {
        return m_results[ 0 ]->get_header( );
    }
    else
    {
        return std::vector<std::string>( );
    }
}

std::vector<std::string>
metaresultfile::get_snp_names()
{
    if( m_results.size( ) > 0 )
    {
        return m_results[ 0 ]->get_snp_names( );
    }
    else
    {
        return std::vector<std::string>( );
    }
}

std::vector<resultfile *> open_result_files(const std::vector<std::string> &paths)
{
    std::vector<resultfile *> result_files;
    for(int i = 0; i < paths.size( ); i++)
    {
        bresultfile *result = new bresultfile( paths[ i ] );
        result_files.push_back( result );
        if( !result->open( ) )
        {
            throw result_open_error( "open_result_files: error: Could not open result file." );
        }
    }
    
    size_t header_size = result_files[ 0 ]->get_header( ).size( );
    for(int i = 1; i < result_files.size( ); i++)
    {
        if( result_files[ i ]->get_header( ).size( ) != header_size )
        {
            throw result_open_error( "open_result_files: error: Different number of columns in result files." );
        }
    }

    return result_files;
}

metaresultfile *open_meta_result_file(const std::vector<std::string> &paths)
{
    std::vector<resultfile *> results = open_result_files( paths );
    return new metaresultfile( results );
}
