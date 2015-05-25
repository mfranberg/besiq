#ifndef __LOGLINEAR_METHOD_H__
#define __LOGLINEAR_METHOD_H__

#include <string>
#include <vector>

#include <armadillo>

#include <besiq/method/method.hpp>
#include <besiq/stats/log_scale.hpp>
#include <besiq/stats/binomial_models.hpp>

/**
 * This class is responsible for intializing and repeatedly
 * executing the log-linear method on pairs of snps.
 */
class loglinear_method
: public method_type
{
public:
    /**
     * Constructor.
     *
     * @param data Additional data required by all methods, such as
     *             covariates.
     */
    loglinear_method(method_data_ptr data);
    
    /**
     * @see method_type::init.
     */
    virtual std::vector<std::string> init();

    /**
     * @see method_type::run.
     */
    virtual void run(const snp_row &row1, const snp_row &row2, float *output);

private: 
    /**
     * A weight > 0 associated with each sample, that allows for
     * covariate adjustment.
     */
    arma::vec m_weight;

    /**
     * The models used.
     */
    std::vector<closed_form_model *> m_models;
};

#endif /* End of __LOGLINEAR_METHOD_H__ */
