#ifndef __CASEONLY_METHOD_H__
#define __CASEONLY_METHOD_H__

#include <string>
#include <vector>

#include <armadillo>

#include <besiq/method/method.hpp>
#include <besiq/method/wald_method.hpp>
#include <besiq/stats/log_scale.hpp>

/**
 * This class is responsible for intializing and repeatedly
 * executing the case-only test by Lewinger et al.
 */
class caseonly_method
: public method_type
{
public:
    /**
     * Constructor.
     *
     * @param data Additional data required by all methods, such as
     *             covariates.
     * @param method The method to use, 'r2' or 'css'.
     */
    caseonly_method(method_data_ptr data, const std::string &method);
    
    /**
     * @see method_type::init.
     */
    virtual std::vector<std::string> init();

    /**
     * @see method_type::run.
     */
    virtual void run(const snp_row &row1, const snp_row &row2, float *output);

private: 
    virtual void compute_r2(const snp_row &row1, const snp_row &row2, float *output);
    virtual void compute_css(const snp_row &row1, const snp_row &row2, float *output);
    virtual void compute_contrast(const snp_row &row1, const snp_row &row2, float *output);
    /**
     * A weight > 0 associated with each sample, that allows for
     * covariate adjustment.
     */
    arma::vec m_weight;

    /**
     * What type of method to use 'r2' or 'css'.
     */
    std::string m_method;

    /**
     * To run the LR test.
     */
    wald_method m_wald;
};

#endif /* End of __CASEONLY_METHOD_H__ */
