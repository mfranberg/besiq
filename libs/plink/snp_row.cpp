#include <plink/snp_row.hpp>

snp_row::snp_row()
{

}

void
snp_row::resize(size_t new_size)
{
    unsigned int total_number_of_bits = new_size * 2;
    unsigned int bits_per_element = 8 * sizeof( int );
    unsigned int elements_required = (total_number_of_bits + bits_per_element - 1) / bits_per_element;

    m_size = new_size;
    m_genotypes.resize( elements_required );
}

size_t
snp_row::size() const
{
    return m_size;
}

unsigned char
snp_row::operator[](size_t index) const
{
    size_t element = ( index * 2 ) / ( 8 * sizeof( unsigned int ) );
    unsigned int element_index = ( index * 2 ) - element * ( 8 * sizeof( unsigned int ) );
    unsigned int element_mask = 0x3 << element_index;
    
    return ( m_genotypes[ element ] & element_mask ) >> element_index;
}

void
snp_row::assign(size_t index, unsigned char value)
{
    size_t element = ( index * 2 ) / ( 8 * sizeof( unsigned int ) );
    unsigned int element_index = ( index * 2 ) - element * ( 8 * sizeof( int ) );
    unsigned int element_mask = ~( 0x3 << element_index );
    unsigned int positioned_value = ( value & 0x3 ) << element_index;

    m_genotypes[ element ] = ( m_genotypes[ element ] & element_mask ) | positioned_value;
}
