#ifndef __LM_STEPWISE_METHOD_H__
#define __LM_STEPWISE_METHOD_H__

#include <string>
#include <vector>

#include <armadillo>

#include <bayesic/method/method.hpp>
#include <bayesic/stats/log_scale.hpp>
#include <bayesic/stats/stepwise_models.hpp>

/**
 * This class is responsible for intializing and repeatedly
 * executing a stepwise method where p-values for a series
 * of more complex models are computed.
 */
class lm_stepwise_method
: public method_type
{
public:
    /**
     * Constructor.
     *
     * @param data Additional data required by all methods, such as
     *             covariates.
     */
    lm_stepwise_method(method_data_ptr data);
    
    /**
     * @see method_type::init.
     */
    virtual void init(std::ostream &output);

    /**
     * Returns the number of usable samples.
     *
     * @param row1 The first snp.
     * @param row2 The second snp.
     * @param phenotype The phenotype.
     * 
     * @return The number of usable samples.
     */    
    unsigned int num_ok_samples(const snp_row &row1, const snp_row &row2, const arma::vec &phenotype);

    /**
     * @see method_type::run.
     */
    virtual void run(const snp_row &row1, const snp_row &row2, std::ostream &output);

private: 
    /**
     * A weight > 0 associated with each sample, that allows for
     * covariate adjustment.
     */
    arma::vec m_weight;

    /**
     * The models used.
     */
    std::vector<stepwise_model *> m_models;
};

#endif /* End of __LM_STEPWISE_METHOD_H__ */
