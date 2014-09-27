#ifndef __STEPWISE_MODELS_H__
#define __STEPWISE_MODELS_H__

#include <armadillo>

#include <plink/snp_row.hpp>
#include <bayesic/stats/log_scale.hpp>
#include <bayesic/stats/snp_count.hpp>

/**
 * This class represents a general model that can compute a
 * model likelihood for two snps.
 */
class stepwise_model
{
public:
    /**
     * Constructor.
     *
     * @param num_params The number of parameters.
     * @param df The degrees of freedom left compared to the full model.
     */
    stepwise_model(unsigned int df)
    : m_df( df )
    {
    }
    
    /**
     * Returns the degrees of freedom left when compared to a
     * full model.
     *
     * @return The degrees of freedom left when compared to a
     * full model.
     */
    unsigned int df()
    {
        return m_df;
    }

    /**
     * Computes the likelihood of the snps and phenotypes under this
     * model. 
     *
     * @param row1 The first snp.
     * @param row2 The second snp.
     * @param phenotype The phenotype, discrete 0.0 and 1.0.
     * @param weight A weight for each sample, this will be used instead of 1.0 as a count.
     * @param is_valid True if model could be estimated, false otherwise.
     *
     * @return The likelihood of the snps.
     */
    virtual log_double prob(const snp_row &row1, const snp_row &row2, const arma::vec &phenotype, const arma::vec &weight, bool *is_valid) = 0;
    
private:
    /**
     * Number of parameters.
     */
    unsigned int m_num_params;
    
    /**
     * Degrees of freedom left when compared to a full model.
     */
    unsigned int m_df;
};

/**
 * This model represents a full model in the sense
 * that each cell in the table has it is own parameter.
 */
class lm_full
: public stepwise_model
{
public:
    /**
     * Constructor.
     */
    lm_full();

    /**
     * @see stepwise_model::prob.
     */
    virtual log_double prob(const snp_row &row1, const snp_row &row2, const arma::vec &phenotype, const arma::vec &weight, bool *is_valid);
};

/**
 * This class represents a model where no snp is associated with
 * the phenotype. The snps may however be in ld.
 */
class intercept
: public stepwise_model
{
public:
    /**
     * Constructor.
     */
    intercept();

    /**
     * @see stepwise_model::prob.
     */
    virtual log_double prob(const snp_row &row1, const snp_row &row2, const arma::vec &phenotype, const arma::vec &weight, bool *is_valid);
};

/**
 * This class represents a model where one snp is associated with
 * the phenotype, and the other possibly in ld with the first.
 */
class single
: public stepwise_model
{
public:
    /**
     * Constructor.
     *
     * @param is_first If true the first snp is associated, otherwise the second.
     */
    single(bool is_first);

    /**
     * @see stepwise_model::prob.
     */
    virtual log_double prob(const snp_row &row1, const snp_row &row2, const arma::vec &phenotype, const arma::vec &weight, bool *is_valid);

private:
    /**
     * Indicates which snp is associated with the phenotype, if true
     * the first, false the second.
     */
    bool m_is_first;
};

#endif /* End of __STEPWISE_MODELS_H__ */
