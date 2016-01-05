#ifndef __LOGP_GRID_H__
#define __LOGP_GRID_H__

#include <iostream>
#include <vector>
#include <string>

#include <plink/plink_file.hpp>

/**
 * Helper class to store summary stats for each
 * grid point.
 */
struct grid_data
{
    grid_data()
        : max_value( 0.0 ),
          sum( 0.0 ),
          sum_sq( 0.0 ),
          n( 0 )
    {
    }

    /**
     * Maximum statistic.
     */
    double max_value;

    /**
     * Sum of statistics.
     */ 
    double sum;

    /**
     * Squared sum of statistics.
     */
    double sum_sq;

    /**
     * Number of observed statistics.
     */
    size_t n;
};

/**
 * Manages summary statistics of p-values in a grid, to enable
 * high level plots of pair-wise interactions.
 */
class logp_grid
{
public:
    /**
     * Constructor.
     *
     * @param loci A list of loci. It is important that the chromosome
     *             and position of these loci are correct. A copy is made
     *             to sort the variants first.
     * @param max_grid_size Maximum grid size.
     * @param bp_window Size in bp of each grid point.
     */
    logp_grid(std::vector<pio_locus_t> loci, size_t max_grid_size = 7000, size_t bp_window = 500000);

    /**
     * Adds a p-value of a variant pair to the grid.
     *
     * @param snp1 Name of the fist variant (must exist in the list of variants
     *                 supplied to the constructor).
     * @param snp2 Name of the second variant (must exist in the list of variants
     *                 supplied to the constructor).
     * @param p The p-value
     *
     * @return True if successful, False otherwise.
     */
    bool add_pvalue(const std::string &snp1, const std::string &snp2, double p);

    /**
     * Writes the grid to a file on the csv format
     *
     * grid_x grid_y sum max sum_sq n
     *
     * @param stream The stream to write the data to.
     */
    void write_grid(std::ostream &stream);

private:
    /**
     * Maps a variant name to its grid point.
     */
    std::map<std::string, int> m_variant_map;

    /**
     * A grid of summary stats.
     */
    std::vector< std::vector<grid_data> > m_grid;
};

#endif /* End of __LOGP_GRID_H__ */
