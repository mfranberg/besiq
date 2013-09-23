#include <plink/snp_row.hpp>

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

const unsigned char &
snp_row::operator[](size_t index) const
{
    return m_genotypes[ index ];
}

unsigned char &
snp_row::operator[](size_t index)
{
    return m_genotypes[ index ];
}
