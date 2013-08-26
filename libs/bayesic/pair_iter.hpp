#ifndef __PAIR_ITER_H__
#define __PAIR_ITER_H__

#include <iostream>
#include <map>
#include <utility>
#include <vector>

#include <plinkio/plinkio.h>

/**
 * The purpose of this class is to iterate over pairs, but
 * not keep them all in memory, as this list can be very large.
 */
class pair_iter
{
public:
    /**
     * Constructor.
     *
     * @param path Path to the file containing pairs.
     * @param loci The list of loci.
     */
    pair_iter(std::istream &pair_stream, const std::vector<pio_locus_t> &loci);

    /**
     * Retrives a pair and writes it.
     *
     * @param next_pair The pair will be written here.
     *
     * @return Returns true if a pair could be fetched, false otherwise, if
     *         it returns false there are no more pairs to be fetched.
     */
    bool get_pair(std::pair<size_t, size_t> *next_pair);

private:
    /**
     * The stream to fetch pairs from.
     */
    std::istream &m_pair_stream;

    /**
     * A map from locus name to index.
     */
    std::map<std::string, size_t> m_snp_to_index;
};

#endif /* End of __PAIR_ITER_H__ */
