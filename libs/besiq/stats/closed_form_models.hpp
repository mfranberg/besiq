#ifndef __CLOSED_FORM_MODELS_H__
#define __CLOSED_FORM_MODELS_H__

#include <armadillo>

#include <besiq/stats/log_scale.hpp>
#include <besiq/stats/snp_count.hpp>

/**
 * This class represents a general model that can compute a
 * model likelihood for two snps.
 */
class closed_form_model
{
public:
    /**
     * Constructor.
     *
     * @param df The degrees of freedom.
     */
    closed_form_model(unsigned int df)
      : m_df( df )
    {
    }

    /**
     * Returns the degrees of freedom of this model.
     *
     * @return The degrees of freedom of this model.
     */
    unsigned int df()
    {
        return m_df;
    }

    /**
     * Computes the likelihood of the snps and phenotypes under this
     * model. 
     *
     * @param count The number of counts for each genotype.
     *
     * @return The likelihood of the snps.
     */
    virtual log_double prob(const arma::mat &count) = 0;
    
private:
    
    /**
     * Degrees of freedom.
     */
    unsigned int m_df;
};

#endif /* End of __CLOSED_FORM_MODELS_H__ */
