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
        : m_data( data )
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
    
    virtual size_t num_usable_samples(const snp_row &row1, const snp_row &row2)
    {
        unsigned int n = 0;
        for(int i = 0; i < row1.size( ); i++)
        {
            if( row1[ i ] != 3 && row2[ i ] != 3 && m_data->missing[ i ] == 0 )
            {
                n++;
            }
        }

        return n;
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
     */
    virtual void run(const snp_row &row1, const snp_row &row2, float *output) = 0;

private:
    /**
     * Additional data required by the method.
     */
    method_data_ptr m_data;
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
