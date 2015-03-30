#include <bayesic/io/misc.hpp>
#include <bayesic/io/resultfile.hpp>

bresultfile::bresultfile(const std::string &path)
    : m_mode( "r" ),
      m_path( path ),
      m_fp( NULL )
{
    
}

bresultfile::bresultfile(const std::string &path, const std::vector<std::string> &snp_names)
    : m_mode( "w" ),
      m_path( path ),
      m_fp( NULL ),
      m_snp_names( snp_names )
{
    for(int i = 0; i < snp_names.size( ); i++)
    {
        m_snp_to_index[ snp_names[ i ] ] = i;
    }
}

bresultfile::~bresultfile()
{
    close( );
}

bool
bresultfile::open()
{
    if( m_fp != NULL )
    {
        fseek( m_fp, 0L, SEEK_SET );
    }

    m_header.version = RESULT_CUR_VERSION;
    m_header.format = 0;
    m_header.snp_names_length = 0;
    m_header.col_names_length = 0;
    m_header.num_pairs = 0;
    m_header.num_float_cols = 0;

    m_fp = fopen( m_path.c_str( ), m_mode.c_str( ) );
    if( m_fp == NULL )
    {
        return false;
    }

    if( m_mode == "r" )
    {
        size_t bytes_read = fread( &m_header, sizeof( result_header ), 1, m_fp );
        if( bytes_read != 1 || m_header.version != RESULT_CUR_VERSION )
        {
            fclose( m_fp );
            m_fp = NULL;
            return false;
        }

        char *buffer = (char *) malloc( m_header.snp_names_length );
        bytes_read = fread( buffer, 1, m_header.snp_names_length, m_fp );
        if( bytes_read != m_header.snp_names_length )
        {
            free( buffer );
            fclose( m_fp );
            m_fp = NULL;
            return false;
        }

        m_snp_names = unpack_string( buffer );
        free( buffer );

        buffer = (char *) malloc( m_header.col_names_length );
        bytes_read = fread( buffer, 1, m_header.col_names_length, m_fp );
        if( bytes_read != m_header.col_names_length )
        {
            free( buffer );
            fclose( m_fp );
            m_fp = NULL;
            return false;
        }

        m_col_names = unpack_string( buffer );
    }
    else
    {
        return set_header( m_col_names );
    }

    return m_fp != NULL;
}

bool
bresultfile::read(std::pair<std::string, std::string> *pair, float *values)
{
    if( m_fp == NULL || m_mode != "r" )
    {
        return false;
    }

    uint32_t snps[ 2 ];
    size_t bytes_read = fread( snps, sizeof( uint32_t ), 2, m_fp );
    if( bytes_read != 2 )
    {
        return false;
    }
    
    pair->first = m_snp_names[ snps[ 0 ] ];
    pair->second = m_snp_names[ snps[ 1 ] ];

    bytes_read = fread( values, sizeof( float ), m_header.num_float_cols, m_fp );
    if( bytes_read != m_header.num_float_cols )
    {
        return false;
    }

    return true;
}

bool
bresultfile::write(const std::pair<std::string, std::string> &pair, float *values)
{
    if( m_mode != "w" || m_fp == NULL )
    {
        return false;
    }

    std::map<std::string, size_t>::const_iterator snp1 = m_snp_to_index.find( pair.first );
    std::map<std::string, size_t>::const_iterator snp2 = m_snp_to_index.find( pair.second );

    if( ( snp1 == m_snp_to_index.end( ) ) || ( snp2 == m_snp_to_index.end( ) ) )
    {
        return false;
    }

    uint32_t write_pair[] = { snp1->second, snp2->second };
    size_t n_snp = fwrite( write_pair, sizeof( uint32_t ), 2, m_fp );
    size_t n_cols = fwrite( values, sizeof( float ), m_header.num_float_cols, m_fp );
    if( n_snp == 2 && n_cols == m_header.num_float_cols )
    {
        m_header.num_pairs++;
    }
    else
    {
        return false;
    }

    return true;
}

void
bresultfile::close()
{
    if( m_fp != NULL )
    {
        if( m_mode == "w" )
        {
            fseek( m_fp, 0L, SEEK_SET );
            fwrite( &m_header, sizeof( result_header ), 1, m_fp );
        }

        fclose( m_fp );
        m_fp = NULL;
    }
}

uint64_t
bresultfile::num_pairs()
{
    if( m_fp != NULL )
    {
        return m_header.num_pairs;
    }
    else
    {
        return 0;
    }
}

const std::vector<std::string> &
bresultfile::get_header()
{
    return m_col_names;
}
const std::vector<std::string> &
bresultfile::get_snp_names()
{
    return m_snp_names;
}

bool
bresultfile::set_header(const std::vector<std::string> &col_names)
{
    if( m_fp == NULL )
    {
        return false;
    }

    fseek( m_fp, 0L, SEEK_SET );
    m_col_names = col_names;
    
    std::string packed_snp_names = pack_string( m_snp_names );
    std::string packed_col_names = pack_string( m_col_names );

    m_header.snp_names_length = packed_snp_names.size( ) + 1;
    m_header.col_names_length = packed_col_names.size( ) + 1;
    m_header.num_float_cols = m_col_names.size( );

    size_t bytes_written = fwrite( &m_header, sizeof( result_header ), 1, m_fp );
    if( bytes_written != 1 )
    {
        return false;
    }

    bytes_written = fwrite( packed_snp_names.c_str( ), 1, m_header.snp_names_length, m_fp );
    if( bytes_written != m_header.snp_names_length )
    {
        return false;
    }

    bytes_written = fwrite( packed_col_names.c_str( ), 1, m_header.col_names_length, m_fp );
    if( bytes_written != m_header.col_names_length )
    {
        return false;
    }

    return true;
}

resultfile *
open_result_file(const std::string &path, const std::vector<std::string> &snp_names)
{
    FILE *fp = fopen( path.c_str( ), "r" );
    if( fp == NULL )
    {
        return NULL;
    }

    result_header header;
    size_t bytes_read = fread( &header, sizeof( result_header ), 1, fp );
    if( bytes_read != 1 )
    {
        return NULL;
    }

    if( header.version == RESULT_CUR_VERSION )
    {
        fclose( fp );
        return new bresultfile( path );
    }
    else
    {
        return NULL;
    }
}

float
result_get_missing()
{
    return -9.0f;
}
