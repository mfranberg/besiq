#ifndef __METHOD_H__
#define __METHOD_H__

#include <vector>
#include <string>

#include <armadillo>

#include <plink/snp_row.hpp>
#include <shared_ptr/shared_ptr.hpp>

class pairfile;
class resultfile;
class genotype_matrix;
typedef shared_ptr<genotype_matrix> genotype_matrix_ptr;

/**
 * Smallest number of samples with a specific genotype for case/control data,
 * for the normal statistical tests to be reasonably valid.
 */
const unsigned int METHOD_SMALLEST_CELL_SIZE_BINOMIAL = 5;

/**
 * Smallest number of samples with a specific genotype for quantative data,
 * for the normal statistical tests to be reasonably valid.
 */
const unsigned int METHOD_SMALLEST_CELL_SIZE_NORMAL = 10;

/**
 * Represents additional data that is required by the method.
 */
struct method_data
{
    /**
     * The phenotype or outcome.
     */
    arma::vec phenotype;

    /**
     * The covariates.
     */
    arma::mat covariate_matrix;

    /**
     * Missing samples are indicated by non-zero values.
     */
    arma::uvec missing;

    /**
    * The number of interactions to correct for.
    */
    unsigned int num_interactions;

    /**
     * The number of single snps.
     */
    unsigned int num_single;

    /**
     * If true print all parameter estimates.
     */
    bool print_params;

    /**
     * The probability that any one snp is associated.
     */
    float single_prior;

    /**
     * Threshold for filtering variant pairs.
     */
    double threshold;

    /**
     * Use fast less robust matrix inversion.
     */
    bool fast_inversion;
};

/**
 * Pointer to method data.
 */
typedef shared_ptr<method_data> method_data_ptr;

class method_type
{
public:
    /**
     * Constructor.
     */
    method_type(method_data_ptr data)
        : m_data( data ),
          m_num_ok_samples( 0 )
    {
    }

    virtual ~method_type()
    {
    }

    /**
     * Returns the additional data.
     */
    method_data_ptr get_data()
    {
        return m_data;
    }

    virtual void set_num_ok_samples(size_t ok_samples)
    {
        m_num_ok_samples = ok_samples;
    }
    
    /**
     * Return the number of samples that could be used in
     * the last call to run.
     *
     * @param row1 The first variant.
     * @param row2 The second variant.
     *
     * @return The number of samples that could be used.
     */
    virtual size_t num_ok_samples(const snp_row &row1, const snp_row &row2)
    {
        return m_num_ok_samples;
    }

    /**
     * Return the column names that will be written by this method.
     */
    virtual std::vector<std::string> init() = 0;

    /**
     * 
     * @param row1 The first genotype.
     * @param row2 The second genotype.
     * @param output The results for this method, the size will be the same
     *               as the length of the header, -9 is interpreted as missing.
     *               You can assume that all values are initialized to -9.
     *
     * @return The value of the test statistic, should return -9 if not computed or missing.
     */
    virtual double run(const snp_row &row1, const snp_row &row2, float *output) = 0;

private:
    /**
     * Additional data required by the method.
     */
    method_data_ptr m_data;

    /**
     * The number of samples
     */
    size_t m_num_ok_samples;
};

/**
 * Runs the given method on the genotype file, traversing
 * the given list of SNPs.
 *
 * @param method A method to run.
 * @param genotype_matix Genotypes for all SNPs.
 * @param pairs The pairs to test.
 * @param result The result file.
 */
void run_method(method_type &method, genotype_matrix_ptr genotype_matrix, pairfile &pairs, resultfile &result);

#endif /* End of __METHOD_H__ */
