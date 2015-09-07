#include <plink/plink_file.hpp>

plink_file::plink_file(const pio_file_t &file, const std::vector<pio_sample_t> &samples, const std::vector<pio_locus_t> &loci, bool mafflip)
    : m_file( file ),
      m_samples( samples ),
      m_loci( loci ),
      m_mafflip( mafflip )
{
    m_row_buffer = (snp_t *) malloc( sizeof( snp_t ) * samples.size( ) );
}

const std::vector<pio_sample_t> &
plink_file::get_samples() const
{
    return m_samples;
}

std::vector<std::string>
plink_file::get_sample_iids() const
{
    std::vector<std::string> iids;
    for(int i = 0; i < m_samples.size( ); i++)
    {
        iids.push_back( m_samples[ i ].iid );
    }

    return iids;
}

const std::vector<pio_locus_t> &
plink_file::get_loci() const
{
    return m_loci;
}

std::vector<std::string>
plink_file::get_locus_names() const
{
    std::vector<std::string> loci_names;
    for(int i = 0; i < m_loci.size( ); i++)
    {
        loci_names.push_back( m_loci[ i ].name );
    }

    return loci_names;
}

float compute_maf(snp_t *row, size_t length)
{
    int mac = 0;
    int total = 0;
    for(int i = 0; i < length; i++)
    {
        if( row[ i ] != 3 )
        {
            mac += row[ i ];
            total++;
        }
    }

    return ((float) mac) / ( 2 * total );
}

bool
plink_file::next_row(snp_row &row)
{
    if( pio_next_row( &m_file, m_row_buffer ) == PIO_OK )
    {
        row.resize( pio_row_size( &m_file ) );
        float maf = compute_maf( m_row_buffer, row.size( ) );
        bool flip = m_mafflip && maf > 0.5;
        for(int i = 0; i < row.size( ); i++)
        {
            if( flip && m_row_buffer[ i ] != 3 )
            {
                row.assign( i, 2 - m_row_buffer[ i ] );
            }
            else
            {
                row.assign( i, m_row_buffer[ i ] );
            }
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

plink_file_ptr
open_plink_file(const std::string &plink_prefix, bool mafflip)
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

    return plink_file_ptr( new plink_file( file, samples, loci, mafflip ) );
}

bool
get_snp_row(plink_file *file, snp_row &row)
{
    return file->next_row( row );
}

genotype_matrix_ptr
create_genotype_matrix(plink_file_ptr genotype_file)
{
    shared_ptr< std::vector<snp_row> > genotypes( new std::vector<snp_row>( ) );
    snp_row row;
    while( genotype_file->next_row( row ) )
    {
        genotypes->push_back( row );
    }

    return genotype_matrix_ptr( new genotype_matrix( genotypes, genotype_file->get_locus_names( ) ) );
}

genotype_matrix::genotype_matrix(shared_ptr< std::vector<snp_row> > matrix, const std::vector<std::string> &snp_names) 
    : m_matrix( matrix ),
    m_snp_names( snp_names )
{
    for(int i = 0; i < snp_names.size( ); i++)
    {
        m_snp_to_index[ snp_names[ i ] ] = i;
    }
}

snp_row const *
genotype_matrix::get_row(const std::string &name) const
{
    std::map<std::string, size_t>::const_iterator it = m_snp_to_index.find( name );
    if( it != m_snp_to_index.end( ) )
    {
        return &(*m_matrix)[ it->second ];
    }
    else
    {
        return NULL;
    }
}
snp_row &
genotype_matrix::get_row(size_t index) const
{
    return (*m_matrix)[ index ];
}

const std::vector<std::string> &
genotype_matrix::get_snp_names() const
{
    return m_snp_names;
}

size_t 
genotype_matrix::size() const
{
    return m_matrix->size( );
}
