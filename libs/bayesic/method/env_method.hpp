#ifndef __ENV_METHOD_H__
#define __ENV_METHOD_H__

#include <vector>
#include <string>

#include <armadillo>

#include <plink/snp_row.hpp>
#include <shared_ptr/shared_ptr.hpp>
#include <bayesic/pair_iter.hpp>
#include <bayesic/method/method.hpp> 

class method_env_type
{
public:
    /**
     * Constructor.
     */
    method_env_type(method_data_ptr data)
        : m_data( data )
    {
    }

    /**
     * Returns the additional data.
     */
    method_data_ptr get_data()
    {
        return m_data;
    }
    
    virtual size_t num_usable_samples(const snp_row &row1)
    {
        unsigned int n = 0;
        for(int i = 0; i < row1.size( ); i++)
        {
            if( row1[ i ] != 3 && m_data->missing[ i ] == 0 )
            {
                n++;
            }
        }

        return n;
    }

    /**
     * Outputs the column names that will be
     * outputted by this method separated by '\t'.
     */
    virtual void init(std::ostream &output) = 0;


    /**
     * @param row The first genotype.
     * @param output The output stream
     */
    virtual void run(const snp_row &row, std::ostream &output) = 0;

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
 * @param loci The list of loci.
 */
void run_env_method(method_env_type &method, const std::vector<snp_row> &genotype_matrix, const std::vector<std::string> &loci);

#endif /* End of __ENV_METHOD_H__ */
