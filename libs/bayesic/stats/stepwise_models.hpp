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
     * @param count A 9x3 matrix, where the 3 columns contains: the
     *              number of samples, the phenotype sum and 
     *              squared phenotype sum.
     *
     * @return The likelihood of the snps.
     */
    virtual log_double prob(const arma::mat &count) = 0;
    
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
    virtual log_double prob(const arma::mat &count);
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
    virtual log_double prob(const arma::mat &count);
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
    virtual log_double prob(const arma::mat &count);

private:
    /**
     * Indicates which snp is associated with the phenotype, if true
     * the first, false the second.
     */
    bool m_is_first;
};

#endif /* End of __STEPWISE_MODELS_H__ */
