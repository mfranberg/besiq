#ifndef __PLINK_FILE_H__
#define __PLINK_FILE_H__

#include <stdexcept>
#include <string>
#include <vector>
#include <map>
#include <shared_ptr/shared_ptr.hpp>

#include <plink/snp_row.hpp>
#include <plinkio/plinkio.h>

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

/**
 * Class that represents an opened plink file.
 */
class plink_file
{
public:
    /**
     * Constructor.
     *
     * @param file Opened file.
     * @param samples List of samples in opened file.
     * @param loci List of loci in opened file.
     * @param mafflip If true, all snps will be flipped so that minor allele is 2.
     */
    plink_file(const pio_file_t &file, const std::vector<pio_sample_t> &samples, const std::vector<pio_locus_t> &loci, bool mafflip = false);

    /**
     * Returns a vector that contains information about the
     * samples.
     *
     * @return a vector that contains information about the
     * samples.
     */
    const std::vector<pio_sample_t> & get_samples() const;

    /**
     * Returns a vector that contains the iids of all samples,
     * in the order they appear.
     *
     * @return A vector that contains the names of all samples.
     */
    std::vector<std::string> get_sample_iids() const;

    /**
     * Returns a vector that contains the fids and iids of all
     * samples in the order they appear.
     *
     * @return A vector of pairs that contain (fid, iid) of all samples.
     */
    std::vector< std::pair<std::string, std::string> > get_sample_fid_iid() const;

    /**
     * Returns a vector that contains the names of all loci,
     * in the order they appear.
     *
     * @return a vector that contains the names of all loci.
     */
    std::vector<std::string> get_locus_names() const;

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

    /**
     * Indiciates whether alleles should be coded according to
     * the minor allele.
     */
    bool m_mafflip;
};

class genotype_matrix
{
public:
    /**
     * Constructor.
     *
     * @param matrix The genotypes. This class now takes responsibility
     *               of the matrix.
     */
    genotype_matrix(shared_ptr< std::vector<snp_row> > matrix, const std::vector<std::string> &snp_names);

    /**
     * Returns the genotypes for the given name.
     *
     * @param name Name of the variant.
     *
     * @return the genotypes for the given index, or NULL
     *         if no genotypes were found.
     */
    snp_row const *get_row(const std::string &name) const;
    
    /**
     * Returns the genotypes for the given index.
     *
     * @param index Name of the variant.
     *
     * @return the genotypes for the given index.
     */
    snp_row &get_row(size_t index) const;

    /**
     * Returns a list of snp names.
     *
     * @return a list of snp names.
     */
    const std::vector<std::string> &get_snp_names() const;

    /**
     * Get the index for a given name.
     *
     * @param name The snp name.
     *
     * @return The index of the snp, -1 if not found.
     */
    int get_index(const std::string &name) const;

    /**
     * Returns the size of the matrix.
     */
    size_t size() const;

private:
    /**
     * The underlying matrix.
     */
    shared_ptr< std::vector<snp_row> > m_matrix;

    /**
     * List of snp names for each row.
     */
    std::vector<std::string> m_snp_names;

    /**
     * An index mapping snp names to indices.
     */
    std::map<std::string, size_t> m_snp_to_index;
};

typedef shared_ptr<genotype_matrix> genotype_matrix_ptr;
typedef shared_ptr<plink_file> plink_file_ptr;

/**
 * Opens the given plink file and returns it.
 *
 * @param plink_prefix The path to the plink file.
 * @param mafflip If true the alleles are coded according to the minor allele.
 *
 * @return A pointer to the plink file.
 */
plink_file_ptr open_plink_file(const std::string &plink_prefix, bool mafflip = false);

/**
 * Retrieves the current row from the opened plink file.
 *
 * @param file An opened plink file.
 * @param row A row object to write to.
 *
 * @return True if the row could be read, false otherwise.
 */
bool get_snp_row(plink_file *file, snp_row &row);

/**
 * Creates a matrix of genotypes by reading the genotypes
 * from the given plink file.
 *
 * @param genotype_file A plink file.
 *
 * @return A matrix of genotypes.
 */
genotype_matrix_ptr create_genotype_matrix(plink_file_ptr genotype_file);

/**
 * Creates a matrix of genotypes by reading the genotypes
 * from the given plink file but filtering on maf.
 *
 * @param genotype_file A plink file.
 * @param maf Minor allele frequency threshold.
 *
 * @return A matrix of genotypes.
 */
genotype_matrix_ptr create_filtered_genotype_matrix(plink_file_ptr genotype_file, float maf);

#endif /* End of __PLINK_FILE_H__ */
