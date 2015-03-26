#include <iterator>
#include <sstream>
#include <iostream>
#include <fstream>

#include <bayesic/pairfile.hpp>

bpairfile::bpairfile(const std::string &path)
    : m_path( path ),
      m_mode( "r" ),
      m_fp( NULL )
{
}

bpairfile::bpairfile(const std::string &path, const std::vector<std::string> &snp_names)
    : m_path( path ),
      m_mode( "w" ),
      m_fp( NULL ),
      m_snp_names( snp_names )
{
}

bpairfile::~bpairfile()
{
    close( );
}

bool
bpairfile::open()
{
    if( m_fp == NULL )
    {
        m_header.version = PAIR_CUR_VERSION;
        m_header.format = 0;
        m_header.num_pairs = 0;
        m_header.header_length = 0;

        if( m_mode == "r" )
        {
            m_fp = fopen( m_path.c_str( ), "r" );
        }
        else if( m_mode == "w" )
        {
            m_fp = fopen( m_path.c_str( ), "w" );
        }
    }
    else
    {
        fseek( m_fp, 0L, SEEK_SET );
    }

    if( m_mode == "r" )
    {
        size_t bytes_read = fread( &m_header, sizeof( bpair_header ), 1, m_fp );
        if( bytes_read != 1 || m_header.version != PAIR_CUR_VERSION )
        {
            fclose( m_fp );
            m_fp = NULL;
            return false;
        }

        char *buffer = (char *) malloc( m_header.header_length );
        bytes_read = fread( buffer, 1, m_header.header_length, m_fp );
        if( bytes_read != m_header.header_length )
        {
            free( buffer );
            fclose( m_fp );
            m_fp = NULL;
            return false;
        }

        m_snp_names = parse_snp_names( buffer );
        free( buffer );
    }
    else
    {
        std::string snp_names = make_snp_names( m_snp_names );
        m_header.header_length = snp_names.size( ) + 1;
        size_t bytes_written = fwrite( &m_header, sizeof( bpair_header ), 1, m_fp );
        if( bytes_written != 1 )
        {
            return false;
        }

        bytes_written = fwrite( snp_names.c_str( ), 1, m_header.header_length, m_fp );
        if( bytes_written != m_header.header_length )
        {
            return false;
        }
    }

    return m_fp != NULL;
}

void
bpairfile::close()
{
    if( m_fp != NULL )
    {
        if( m_mode == "w" )
        {
            fseek( m_fp, 0L, SEEK_SET );
            fwrite( &m_header, sizeof( bpair_header ), 1, m_fp );
        }

        fclose( m_fp );
        m_fp = NULL;
    }
}

const std::vector<std::string> &
bpairfile::get_snp_names()
{
    return m_snp_names;
}

bool
bpairfile::read(std::pair<std::string, std::string> &pair)
{
    if( m_mode != "r" || m_fp == NULL )
    {
        return false;
    }

    uint32_t read_pair[ 2 ];
    size_t bytes_read = fread( read_pair, sizeof( uint32_t ), 2, m_fp );
    if( bytes_read != 2 )
    {
        return false;
    }

    pair.first = m_snp_names[ read_pair[ 0 ] ];
    pair.second = m_snp_names[ read_pair[ 1 ] ];

    return true;
}

bool
bpairfile::write(size_t snp_id1, size_t snp_id2)
{
    if( m_mode != "w" || m_fp == NULL )
    {
        return false;
    }

    uint32_t pair[] = { snp_id1, snp_id2 };
    size_t bytes_written = fwrite( pair, sizeof( uint32_t ), 2, m_fp );

    if( bytes_written == sizeof( uint32_t ) * 2 )
    {
        m_header.num_pairs++;
        return true;
    }
    else
    {
        return false;
    }
}

std::vector<std::string>
bpairfile::parse_snp_names(const char *snp_name_str)
{
    std::string snp_names( snp_name_str );
    std::stringstream strstr( snp_names );
    std::istream_iterator<std::string> it( strstr );
    std::istream_iterator<std::string> end;
    std::vector<std::string> results( it, end );

    return results;
}

std::string
bpairfile::make_snp_names(const std::vector<std::string> &snp_names)
{
    std::stringstream output;
    output << snp_names[ 0 ];
    for(size_t i = 1; i < snp_names.size( ); i++)
    {
        output << "\t" << snp_names[ i ];
    }

    return output.str( );
}

size_t
bpairfile::num_pairs()
{
    return m_header.num_pairs;
}

tpairfile::tpairfile(const std::string &path, std::vector<std::string> snp_names, const char *mode)
    : m_path( path ),
      m_mode( mode ),
      m_input( NULL ),
      m_output( NULL ),
      m_snp_names( snp_names )
{
    for(int i = 0; i < m_snp_names.size( ); i++)
    {
        m_snp_to_index[ snp_names[ i ] ] = i;
    }
}

tpairfile::~tpairfile()
{
    close( );
}

bool
tpairfile::open()
{
    if( m_mode == "r" )
    {
        if( m_input != NULL )
        {
            delete m_input;
            m_input = NULL;
        }

        if( m_path == "-" )
        {
            m_input = &std::cin;
        }
        else
        {
            m_input = new std::ifstream( m_path.c_str( ) );
        }

        return m_input->good( );
    }
    else
    {
        if( m_output != NULL )
        {
            delete m_output;
            m_output = NULL;
        }

        if( m_path == "-" )
        {
            m_output = &std::cout;
        }
        else
        {
            m_output = new std::ofstream( m_path.c_str( ) );
        }

        return m_output->good( );
    }
}

void
tpairfile::close()
{
    if( m_input != NULL && m_input != &std::cin)
    {
        delete m_input;
    }
    if( m_output != NULL && m_output != &std::cout)
    {
        delete m_output;
    }
}

bool tpairfile::read(std::pair<std::string, std::string> &pair)
{
    std::string snp1;
    std::string snp2;

    if( !( *m_input >> snp1 ) || !( *m_input >> snp2 ) )
    {
        return false;
    }

    pair.first = snp1;
    pair.second = snp2;
    
    return true;
}

bool tpairfile::write(size_t snp1_id, size_t snp2_id)
{
    *m_output << m_snp_names[ snp1_id ] << " " << m_snp_names[ snp2_id ] << "\n";

    return true;
}

size_t
tpairfile::num_pairs()
{
    return 0;
}

pairfile *
open_pair_file(const std::string &path, const std::vector<std::string> &snp_names)
{
    FILE *fp = fopen( path.c_str( ), "r" );
    if( fp == NULL )
    {
        return NULL;
    }

    bpair_header header;
    size_t bytes_read = fread( &header, sizeof( bpair_header ), 1, fp );
    if( bytes_read != 1 )
    {
        return NULL;
    }

    if( header.version == PAIR_CUR_VERSION )
    {
        fclose( fp );
        return new bpairfile( path );
    }
    else
    {
        return new tpairfile( path, snp_names, "r" );
    }
}
