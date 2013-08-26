#include <plink_file.hpp>

plink_file::plink_file(const pio_file_t &file, const std::vector<pio_sample_t> &samples, const std::vector<pio_locus_t> &loci)
    : m_file( file ),
      m_samples( samples ),
      m_loci( loci )
{
    m_row_buffer = (snp_t *) malloc( sizeof( snp_t ) * samples.size( ) );
}

const std::vector<pio_sample_t> &
plink_file::get_samples() const
{
    return m_samples;
}

const std::vector<pio_locus_t> &
plink_file::get_loci() const
{
    return m_loci;
}

bool
plink_file::next_row(snp_row &row)
{
    if( pio_next_row( &m_file, m_row_buffer ) == PIO_OK )
    {
        row.resize( pio_row_size( &m_file ) );
        for(int i = 0; i < row.size( ); i++)
        {
            row[ i ] = m_row_buffer[ i ];
        }
        return true;
    }
    else
    {
        return false;
    }
}

plink_file::~plink_file()
{
    free( m_row_buffer );
    pio_close( &m_file );
}

snp_row::snp_row()
{
}


void
snp_row::resize(size_t new_size)
{
    m_genotypes.resize( new_size );
}

size_t
snp_row::size() const
{
    return m_genotypes.size( );
}

const snp_t &
snp_row::operator[](size_t index) const
{
    return m_genotypes[ index ];
}

snp_t &
snp_row::operator[](size_t index)
{
    return m_genotypes[ index ];
}

plink_file_ptr
open_plink_file(const std::string &plink_prefix)
{
    pio_file_t file;
    if( pio_open( &file, plink_prefix.c_str( ) ) != PIO_OK )
    {
        throw plink_error( "Couldn't not open file " + plink_prefix );
    }

    std::vector<pio_sample_t> samples;
    for(int i = 0; i < pio_num_samples( &file ); i++)
    {
        samples.push_back( *pio_get_sample( &file, i ) );
    }

    std::vector<pio_locus_t> loci;
    for(int i = 0; i < pio_num_loci( &file ); i++)
    {
        loci.push_back( *pio_get_locus( &file, i ) );
    }

    return plink_file_ptr( new plink_file( file, samples, loci ) );
}

bool
get_snp_row(plink_file *file, snp_row &row)
{
    return file->next_row( row );
}
