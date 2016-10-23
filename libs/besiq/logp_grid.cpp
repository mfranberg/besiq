#include <besiq/logp_grid.hpp>

#include <cmath>

bool locus_cmp(const pio_locus_t &a, const pio_locus_t &b)
{
    if( a.chromosome == b.chromosome )
    {
        return a.bp_position < b.bp_position;
    }
    else
    {
        return a.chromosome < b.chromosome;
    }
}

logp_grid::logp_grid(std::vector<pio_locus_t> loci, size_t max_grid_size, size_t bp_window)
{
    std::sort( loci.begin( ), loci.end( ), locus_cmp );


    unsigned int prev_chrom = 0;
    long long grid_start_pos = -1;
    int grid_point = -1;
    for(int i = 0; i < loci.size( ); i++)
    {
        if( loci[ i ].chromosome != prev_chrom )
        {
            grid_start_pos = loci[ i ].bp_position;
            grid_point++;
        }

        if( loci[ i ].bp_position - grid_start_pos > bp_window )
        {
            grid_start_pos = loci[ i ].bp_position;
            grid_point++;
        }

        m_variant_map[ std::string( loci[ i ].name ) ] = grid_point;
        prev_chrom = loci[ i ].chromosome;
    }

    if( grid_point + 1 > max_grid_size )
    {
        std::cerr << "Too many grid points: " << grid_point + 1 << std::endl;
        exit( 1 );
    }

    for(int i = 0; i < grid_point + 1; i++)
    {
        m_grid.push_back( std::vector<grid_data>( grid_point + 1 ) );
    }
}

bool
logp_grid::add_pvalue(const std::string &snp1, const std::string &snp2, double p)
{
    std::map<std::string, int>::iterator s1 = m_variant_map.find( snp1 );
    std::map<std::string, int>::iterator s2 = m_variant_map.find( snp2 );

    if( s1 == m_variant_map.end( ) || s2 == m_variant_map.end( ) )
    {
        return false;
    }

    grid_data &data1 = m_grid[ s1->second ][ s2->second ];
    grid_data &data2 = m_grid[ s2->second ][ s1->second ];

    double mlogp = (p != 1.0) ? -std::log( p ) : 0.0;

    data1.max_value = std::max( data1.max_value, mlogp );
    data1.sum += mlogp;
    data1.sum_sq += mlogp * mlogp;
    data1.n += 1;
    
    data2.max_value = data1.max_value;
    data2.sum = data1.sum;
    data2.sum_sq = data1.sum_sq;
    data2.n = data1.n;

    return true;
}

void
logp_grid::write_grid(std::ostream &stream)
{
    stream << "x\ty\tmax\tsum\tsum_sq\tn\n";
    for(int i = 0; i < m_grid.size( ); i++)
    {
        std::vector<grid_data> &row = m_grid[ i ];

        for(int j = 0; j < row.size( ); j++)
        {
            stream << i << "\t" << j << "\t" << row[ j ].max_value << "\t" << row[ j ].sum << "\t" << row[ j ].sum_sq << "\t" << row[ j ].n << "\n";
        }
    }
}
