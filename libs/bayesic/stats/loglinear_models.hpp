#ifndef __LOGLINEAR_MODELS_H__
#define __LOGLINEAR_MODELS_H__

#include <armadillo>

#include <plink_file.hpp>
#include <stats/dirichlet.hpp>
#include <stats/log_scale.hpp>
#include <stats/snp_count.hpp>

/**
 * This class represents a general model that can compute a
 * model likelihood for two snps.
 */
class loglinear_model
{
public:
    /**
     * Constructor.
     *
     * @param num_params The number of parameters.
     * @param df The degrees of freedom left compared to the full model.
     */
    loglinear_model(unsigned int num_params, unsigned int df)
    : m_num_params( num_params ),
      m_df( df )
    {
    }
    
    /**
     * Returns the number of parameters in this model.
     *
     * @return the number of parameters in this model.
     */
    unsigned int num_params( )
    {
        return m_num_params;
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
     *
     * @return The likelihood of the snps.
     */
    virtual log_double prob(const snp_row &row1, const snp_row &row2, const arma::vec &phenotype, const arma::vec &weight) = 0;
    
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
class full
: public loglinear_model
{
public:
    /**
     * Constructor.
     */
    full();

    /**
     * @see loglinear_model::prob.
     */
    virtual log_double prob(const snp_row &row1, const snp_row &row2, const arma::vec &phenotype, const arma::vec &weight);
};

/**
 * This class represents a model where no snp is associated with
 * the phenotype. The snps may however be in ld.
 */
class block
: public loglinear_model
{
public:
    /**
     * Constructor.
     */
    block();

    /**
     * @see loglinear_model::prob.
     */
    virtual log_double prob(const snp_row &row1, const snp_row &row2, const arma::vec &phenotype, const arma::vec &weight);
};

/**
 * This class represents a model where one snp is associated with
 * the phenotype, and the other possibly in ld with the first.
 */
class partial
: public loglinear_model
{
public:
    /**
     * Constructor.
     *
     * @param is_first If true the first snp is associated, otherwise the second.
     */
    partial(bool is_first);

    /**
     * @see loglinear_model::prob.
     */
    virtual log_double prob(const snp_row &row1, const snp_row &row2, const arma::vec &phenotype, const arma::vec &weight);

private:
    /**
     * Indicates which snp is associated with the phenotype, if true
     * the first, false the second.
     */
    bool m_is_first;
};

#endif /* End of __LOGLINEAR_MODELS_H__ */
