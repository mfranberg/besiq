#ifndef __PLINK_FILE_H__
#define __PLINK_FILE_H__

#include <stdexcept>
#include <string>
#include <vector>
#include <shared_ptr.hpp>

extern "C"
{
    #include <plinkio/plinkio.h>
}

/**
 * General exception class when dealing with plink files.
 */
class plink_error
: public std::runtime_error
{
public:
    plink_error(std::string const& s)
    : std::runtime_error(s)
    {
    }
};

class snp_row
{
public:
    /**
     * Constructor.
     */
    snp_row();

    /**
     * Resizes the row to be able to hold the given size.
     *
     * @param new_size The new size of the row.
     */
    void resize(size_t new_size);

    /**
     * Returns the length of the row.
     *
     * @return The length of the row.
     */
    size_t size() const;

    /**
     * Access operator, not checked for bounds.
     *
     * @param index Index of the SNP to retrive.
     * 
     * @return Return the SNP at the given index.
     */
    const snp_t & operator[](size_t index) const;

    /**
     * Access operator for assignment.
     *
     * @param index Index of the SNP to modify.
     *
     * @return A reference to the SNP to modify.
     */
    snp_t & operator[](size_t index);

private:
    /**
     * Internal data structure, as a vector.
     */
    std::vector<snp_t> m_genotypes;
};

/**
 * Class that represents an opened plink file.
 */
class plink_file
{
public:
    /**
     * Constructor.
     *
     *
     */
    plink_file(const pio_file_t &file, const std::vector<pio_sample_t> &samples, const std::vector<pio_locus_t> &loci);

    /**
     * Returns a vector that contains information about the
     * samples.
     *
     * @return a vector that contains information about the
     * samples.
     */
    const std::vector<pio_sample_t> & get_samples() const;

    /**
     * Reads a row from the underlying plink file, and writes
     * it in the given snp_row.
     *
     * @param row Row to write genotypes in.
     *
     * @return True if the row could be read, false otherwise.
     */
    bool next_row(snp_row &row);

    /**
     * Returns a vector that contains information about the
     * SNPs.
     *
     * @return a vector that contains information about the
     * SNPs.
     */
    const std::vector<pio_locus_t> & get_loci() const;

    /**
     * Destructor.
     *
     * Deallocates read buffer and closes plink file.
     */
    ~plink_file();

private:
    /**
     * Underlying file object.
     */
    pio_file_t m_file;

    /**
     * List of samples.
     */
    std::vector<pio_sample_t> m_samples;

    /**
     * List of SNPs.
     */
    std::vector<pio_locus_t> m_loci;

    /**
     * Buffer for reading rows.
     */
    snp_t *m_row_buffer;
    
};

typedef shared_ptr<plink_file> plink_file_ptr;

/**
 * Opens the given plink file and returns it.
 *
 * @param plink_prefix The path to the plink file.
 *
 * @return A pointer to the plink file.
 */
plink_file_ptr open_plink_file(const std::string &plink_prefix);

/**
 * Retrieves the current row from the opened plink file.
 *
 * @param file An opened plink file.
 * @param row A row object to write to.
 *
 * @return True if the row could be read, false otherwise.
 */
bool get_snp_row(plink_file *file, snp_row &row);

#endif /* End of __PLINK_FILE_H__ */
