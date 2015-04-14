#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>

#include <bayesic/io/misc.hpp>
#include <bayesic/io/pairfile.hpp>

#define BLOCK_SIZE 4194304ULL

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
bpairfile::open(size_t split, size_t num_splits)
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

        m_snp_names = unpack_string( buffer );
        free( buffer );

        /* Only read a part of the pair file */
        uint64_t pairs_per_split = ( m_header.num_pairs + num_splits - 1 ) / num_splits;
        uint64_t seek_length = sizeof( uint32_t ) * 2 * pairs_per_split * (split - 1);
        m_pairs_left = pairs_per_split;

        if( seek_length > 0 )
        {
            fseek( m_fp, seek_length, SEEK_CUR );
        }
    }
    else
    {
        std::string snp_names = pack_string( m_snp_names );
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
    if( m_mode != "r" || m_fp == NULL || m_pairs_left <= 0 )
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
    m_pairs_left--;

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

    if( bytes_written == 2 )
    {
        m_header.num_pairs++;
        return true;
    }
    else
    {
        return false;
    }
}

size_t
bpairfile::num_pairs()
{
    return m_header.num_pairs;
}

tpairfile::tpairfile(const std::string &path, std::vector<std::string> snp_names, const char *mode)
    : m_path( path ),
      m_mode( mode ),
      m_num_pairs( 0 ),
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
tpairfile::open(size_t split, size_t num_splits)
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
            if( num_splits > 1 )
            {
                uint64_t npairs = num_pairs( );
                uint64_t pairs_per_split = ( npairs + num_splits - 1 ) / num_splits;
                uint64_t pairs_to_skip = pairs_per_split * (split - 1);

                std::pair<std::string, std::string> pair;
                while( pairs_to_skip-- > 0 )
                {
                    read( pair );
                }
                m_pairs_left = pairs_per_split;
            }
            else
            {
                m_pairs_left = -1;
            }
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
    if( m_path != "-" && m_pairs_left <= 0 )
    {
        return false;
    }

    std::string snp1;
    std::string snp2;

    if( !( *m_input >> snp1 ) || !( *m_input >> snp2 ) )
    {
        return false;
    }

    pair.first = snp1;
    pair.second = snp2;
    m_pairs_left--;
    
    return true;
}

bool tpairfile::write(size_t snp1_id, size_t snp2_id)
{
    *m_output << m_snp_names[ snp1_id ] << " " << m_snp_names[ snp2_id ] << "\n";
    m_num_pairs++;

    return true;
}

size_t
tpairfile::num_pairs()
{
    if( m_path != "-" && m_mode == "r" && m_num_pairs == 0 )
    {
        std::ifstream *input = static_cast<std::ifstream *>( m_input );
        std::streampos pos = input->tellg( );
        input->seekg( 0, input->beg );
        
        std::pair<std::string, std::string> pair;
        while( read( pair ) )
        {
            m_num_pairs++;
        }
       
        input->clear( );
        input->seekg( pos, input->beg );
    }

    return m_num_pairs;
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

    fclose( fp );
    if( header.version == PAIR_CUR_VERSION )
    {
        return new bpairfile( path );
    }
    else
    {
        return new tpairfile( path, snp_names, "r" );
    }
}

bool write_pairs_in_block(uint32_t *snp_buffer, uint64_t num_pairs, FILE *in_fp, FILE *out_fp)
{
    for(int64_t block_pairs_left = num_pairs; block_pairs_left > 0; block_pairs_left -= BLOCK_SIZE )
    {
        uint64_t cur_block_size = 2 * BLOCK_SIZE;
        if( block_pairs_left < BLOCK_SIZE )
        {
            cur_block_size = 2 * block_pairs_left;
        }

        size_t bytes;
        if( (bytes = fread( snp_buffer, sizeof( uint32_t ), cur_block_size, in_fp ) ) != cur_block_size )
        {
            return false;
        }
        if( fwrite( snp_buffer, sizeof( uint32_t ), cur_block_size, out_fp ) != cur_block_size )
        {
            return false;
        }
    }
    return true;
}

char * parse_header(FILE *fp, bpair_header *header)
{
    /* Parse header info */
    size_t bytes_read = fread( header, sizeof( bpair_header ), 1, fp );
    if( bytes_read != 1 || header->version != PAIR_CUR_VERSION )
    {
        fclose( fp );
        return NULL;
    }
    
    char *buffer = (char *) malloc( header->header_length );
    bytes_read = fread( buffer, 1, header->header_length, fp );
    if( bytes_read != header->header_length )
    {
        free( buffer );
        fclose( fp );
        return NULL;
    }

    return buffer;
}

bool split_pair_file(const std::string &all_pairs, size_t num_splits, const std::string &output_path)
{
    FILE *fp = fopen( all_pairs.c_str( ), "r" );
    if( fp == NULL || num_splits <= 1 )
    {
        return false;
    }

    bpair_header header;
    char *buffer = parse_header( fp, &header );
    if( buffer == NULL )
    {
        fclose( fp );
        return false;
    }

    /* Create the split files and copy snps */
    uint64_t num_pairs_in_each_split = (header.num_pairs + num_splits - 1) / num_splits;
    if( num_pairs_in_each_split <= 0 )
    {
        fclose( fp );
        free( buffer );
        return false;
    }

    bool error = false;
    uint64_t total_pairs = header.num_pairs;
    int split = 1;
    FILE *cur_fp;
    uint32_t *snp_buffer = (uint32_t *) malloc( sizeof( uint32_t ) * BLOCK_SIZE * 2 );
    for(int64_t pairs_left = total_pairs; pairs_left > 0; pairs_left -= num_pairs_in_each_split)
    {
        std::stringstream ss;
        ss << output_path << ".split" << split;
        std::string filename = ss.str( );

        cur_fp = fopen( filename.c_str( ), "w" );
        if( cur_fp == NULL )
        {
            error = true;
            goto split_error;
        }

        bpair_header cur_header( header );
        cur_header.num_pairs = num_pairs_in_each_split;
        if( pairs_left < num_pairs_in_each_split )
        {
            cur_header.num_pairs = pairs_left;
        }

        if( fwrite( &cur_header, sizeof( bpair_header ), 1, cur_fp ) != 1 )
        {
            error = true;
            fclose( cur_fp );
            goto split_error;
        }
        if( fwrite( buffer, 1, cur_header.header_length, cur_fp ) != cur_header.header_length )
        {
            error = true;
            fclose( cur_fp );
            goto split_error;
        }
        if( !write_pairs_in_block( snp_buffer, cur_header.num_pairs, fp, cur_fp ) )
        {
            error = true;
            fclose( cur_fp );
            goto split_error;
        }

        fclose( cur_fp );

        split++;
    }

split_error: 
    fclose( fp );
    free( buffer );
    free( snp_buffer );

    return !error;
}


