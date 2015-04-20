#ifndef __LM_ENV_STEPWISE_H__
#define __LM_ENV_STEPWISE_H__

#include <string>
#include <vector>

#include <armadillo>

#include <glm/glm.hpp>
#include <bayesic/method/env_method.hpp>
#include <bayesic/stats/log_scale.hpp>

/**
 * This class is responsible for initializing and repeatedly
 * executing the lm regression on a single variant and an
 * environmental factor. The variant is treated as a factor
 * variable.
 */
class lm_env_stepwise
: public method_env_type
{
public:
    /**
     * Constructor.
     *
     * @param data Additional data required by all methods.
     * @param E Environment matrix that only contains the
     *          non-redundant levels, i.e. the intercept has
     *          been removed.
     */
    lm_env_stepwise(method_data_ptr data, const arma::mat &E);
    
    /**
     * @see method_env_type::init.
     */
    virtual void init(std::ostream &output);

    /**
     * @see method_env_type::run.
     */
    virtual void run(const snp_row &row, std::ostream &output);

private:
    /**
     * Design matrix corresponding to the 4 hypotheses.
     *
     * null: No snp or environment
     * snp: Only snp effect
     * env: Only env effect
     * add: Additive snp and env
     * alt: Full snp and env model
     */
    arma::mat m_null_matrix;
    arma::mat m_snp_matrix;
    arma::mat m_env_matrix;
    arma::mat m_add_matrix;
    arma::mat m_alt_matrix;

    /**
     * Contains all design matrices.
     */
    std::vector<arma::mat *> m_model;

    /**
     * Only contains non-redundant levels, i.e. the
     * "reference level" or "intercept" is removed 
     *
     * Environment matrix, need to save this to compute
     * the variant interaction features. 
     */
    arma::mat m_E;

    /**
     * Initializes all design matrices with the genotypes.
     *
     * @param row Variant.
     * @param missing The missing samples.
     * @param valid Determines whether there were enough samples in each cell.
     */
    void init_matrix_with_snp(const snp_row &row, arma::uvec &missing, bool *valid);    
};

#endif /* End of __LM_ENV_STEPWISE_H__ */
