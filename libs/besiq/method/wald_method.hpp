#ifndef __WALD_METHOD_H__
#define __WALD_METHOD_H__

#include <string>
#include <vector>

#include <armadillo>

#include <besiq/method/method.hpp>
#include <besiq/stats/log_scale.hpp>

/**
 * This class is responsible for executing the closed form
 * wald test for a logistic regression model.
 */
class wald_method
: public method_type
{
public:
    /**
     * Constructor.
     *
     * @param data Additional data required by all methods, such as
     *             covariates.
     */
    wald_method(method_data_ptr data);
    
    /**
     * @see method_type::init.
     */
    virtual std::vector<std::string> init();
    
    /**
     * @see method_type::run.
     */
    virtual double run(const snp_row &row1, const snp_row &row2, float *output);
private:
    /**
     * Weight for each sample.
     */
    arma::vec m_weight;
};

#endif /* End of __WALD_METHOD_H__ */
