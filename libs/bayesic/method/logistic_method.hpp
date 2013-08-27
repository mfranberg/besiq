#ifndef __LOGISTIC_METHOD_H__
#define __LOGISTIC_METHOD_H__

#include <string>
#include <vector>

#include <armadillo>

#include <models/binomial.hpp>
#include <irls.hpp>
#include <method/method.hpp>
#include <stats/log_scale.hpp>

/**
 * This class is responsible for initializing and repeatedly
 * executing the logistic regression on pairs of snps.
 */
class logistic_method
: public method_type
{
public:
    /**
     * Constructor.
     *
     * @param data Additional data required by all methods.
     */
    logistic_method(method_data_ptr data);
    
    /**
     * @see method_type::run.
     */
    virtual void run(const snp_row &row1, const snp_row &row2, const std::string &name1, const std::string &name2);

private:
    /**
     * The regression design matrix, the first 4 columns are:
     * snp1, snp2, snp1 x snp2 and intercept. The rest are covariates.
     */
    arma::mat m_design_matrix;

    /**
     * The glm model used, in this case a binomial model with logit link.
     */
    binomial m_model;
};

#endif /* End of __LOGISTIC_METHOD_H__ */
