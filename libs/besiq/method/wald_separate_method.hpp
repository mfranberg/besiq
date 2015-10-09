#ifndef __WALD_SEPARATE_METHOD_H__
#define __WALD_SEPARATE_METHOD_H__

#include <string>
#include <vector>

#include <armadillo>

#include <besiq/method/method.hpp>
#include <besiq/stats/log_scale.hpp>

/**
 * This class is responsible for executing the closed form
 * wald test for a logistic regression model.
 */
class wald_separate_method
: public method_type
{
public:
    /**
     * Constructor.
     *
     * @param data Additional data required by all methods, such as
     *             covariates.
     */
    wald_separate_method(method_data_ptr data, bool is_lm);
    
    /**
     * @see method_type::init.
     */
    virtual std::vector<std::string> init();
    
    /**
     * @see method_type::run.
     */
    virtual double run(const snp_row &row1, const snp_row &row2, float *output);
private:
    void compute_lm(const snp_row &row1, const snp_row &row2, float *output);
    void compute_binomial(const snp_row &row1, const snp_row &row2, float *output);
    /**
     * Weight for each sample.
     */
    arma::vec m_weight;

    /**
     * Indicates whether this is a linear model or not.
     */
    bool m_is_lm;
};

#endif /* End of __WALD_SEPARATE_METHOD_H__ */
